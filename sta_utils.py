import subprocess
from scipy import stats, signal
from scipy import interpolate
import numpy as np
import time, math
import sys
import os.path
import multiprocessing as mp
from matplotlib import animation
from matplotlib.path import Path
import glob #Wildcard searching in dir listing
import skimage
from skimage import data
from skimage.transform import rotate
            
from pylab import *
from scipy.interpolate import interp1d
import struct, array, csv
import scipy.optimize 
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline

import matplotlib.mlab as mlab
from scipy import ndimage

import PIL
from PIL import ImageFont
from PIL import Image
from PIL import ImageDraw

from libtfr import *


def Load_images(file_dir, file_name):
    
    if (os.path.exists(file_dir + file_name + '/' + file_name + '_images.npy')==False):

        print file_dir + file_name + '/' + file_name+'.tif'
        img = Image.open(file_dir + file_name + '/' + file_name+'.tif')

        counter=0
        if True:
            while True:
                try:
                    img.seek(counter)
                except EOFError:
                    break
                counter+=1
                print counter
        img.seek(0)

        #Default pic sizes
        n_pixels = img.size[0]

        #Initialize 3D image array
        n_frames = counter
        #n_frames = 69093

        images_raw = np.zeros((n_frames, 256, 256), dtype = np.float16)

        print "n_frames: ", n_frames
        for i in range(0, n_frames,1): 
            try:
                img.seek(i)
                print "Loading frame: ", i
                #images_raw [i] = np.flipud(np.fliplr(np.float16(img))) #FLIP IMAGES FOR Experiments Nov and Dec 2015
                if n_pixels==128:
                    images_raw[i] = scipy.misc.imresize(np.float16(img),2., interp='bicubic', mode=None)
                else:
                    images_raw[i] = np.float16(img) #2016-1-11 2016-1-14 experiment no flipping needed

            except EOFError:
                break

        images_start= 0
        images_end = n_frames
        
        images_raw=images_raw[images_start:images_end]

        print "Saving imaging array..."

        np.save(file_dir + file_name + '/' + file_name + '_images', images_raw)
        np.savetxt(file_dir + file_name + '/' + file_name + '_images_start_'+str(images_start)+'_end_'+
        str(images_end), [images_start, images_end])

        images_raw = np.load(file_dir + file_name + '/' + file_name + '_images.npy')
        images_raw = np.float16(images_raw)
        
        return images_raw

    else:
        images_raw = np.load(file_dir + file_name + '/' + file_name + '_images.npy')
        images_raw = np.float16(images_raw)



    if len(images_raw[0])==256:
        return images_raw

    else:
        #temp_img = np.float64(images_raw) #load stim data
        temp=[]
        for k in range(len(images_raw)):
            print k, len(images_raw)
            temp.append(scipy.misc.imresize(np.float16(images_raw[k]),2., interp='bicubic', mode=None)) #Interpolate
        
        images_raw = np.float16(temp)
        np.save(file_dir + file_name + '/' + file_name + '_images', images_raw)
        
        return np.array(temp)
                    
                    
        
def Rotate_images(images_raw, file_dir, file_name, overwrite_shifted):

    temp_name = file_dir + file_name + '/' + file_name + '_images_rotated'
    
    if overwrite_shifted or (os.path.exists(temp_name + '.npy')==False):
        #Recenter/shift images:
        print "Rotating images"
        with open(file_dir + "image_shift.txt", "r") as f: #Text file contains image shift and rotation angle info
            data = csv.reader(f)
            temp = []
            for row in data:
                temp.append(int(row[0]))
        x_shift = temp[0]
        y_shift = temp[1]
        angle = temp[2]
        
        n_pixels = len(images_raw[0])
        
        shift_img = np.float64(images_raw)

        #Rotate image
        print "Rotate angle: ", angle
        if angle != 0:

            for i in range(len(shift_img)):
                print "Rotating img: ", i
                shift_img[i] = skimage.transform.rotate(shift_img[i], angle)#, mode='constant', cval=100)
                
        print "...done rotating."
        np.save(temp_name,np.array(shift_img,dtype=np.float16))

    else:
        shift_img = np.load(temp_name+'.npy')

    n_pixels = len(shift_img[0])

    return shift_img
            

def Reduce_images_128pixels(images_raw):
    print images_raw.shape
    
    from skimage.measure import block_reduce
    reduced_frames = np.zeros((len(images_raw),128,128),dtype=np.float32)
    print reduced_frames.shape
    for i in range(len(images_raw)):
        print "Block reducing frame: ", i
        reduced_frames[i] = block_reduce(images_raw[i], block_size=(2,2), func=np.mean)

    return reduced_frames

def Luminance_contrast(frame_in, frame_out):
    for i in range(0,len(frame_in),1):
        for j in range(0,len(frame_in[0]),1):
            frame_out[i,j] = frame_in[i][j]/np.mean(frame_in)   #Definition of luminance-contrast from Ringach 2002
    return frame_out
    
def Baseline_regression(images_raw):        
    print "Removing baseline over all data"
    baseline = np.mean(images_temp, axis=0)
    images_temp = (images_temp - baseline)/baseline


def Load_images_start_end(file_dir, file_name, images_raw):
    ''' This method aligns imaging record to electrophysiology by identifying the on and off 
        trigger times from the original Multi-Channel-Systems .mcd 17th channel 
        
        Other electrophysiology files will require different code for identifying the trigger time
    '''
    
    mcd_file = file_dir + file_name+ '/' + file_name + '.mcd'
    if (os.path.exists(mcd_file)==True):
        data = MCD_read_imagingtimes_old(mcd_file)

        print "Finding beginning and end of imaging trigger on channel 17th .mcd file"
        temp_data = []
        for i in range(data['extra'].item_count):
            temp_data.append(data['extra'].get_data(i)) 
        temp_data = np.array(temp_data)[:,0]    #Select time column only from 17th channel

        start_array = []
        end_array = []
        start_array.append(temp_data[0])
        for i in range(1,len(temp_data)-1,1):
            if temp_data[i+1]-temp_data[i]>1.0:
                end_array.append(temp_data[i])
                start_array.append(temp_data[i+1])
        end_array.append(temp_data[-1])

    else:
        #Load epoch.txt file
        epoch_file = file_dir + file_name+ '/epochs.txt'
        data = np.loadtxt(epoch_file)
        start_array = data[:,0]
        end_array = data[:,1]


    #Load index of recording for series recordings
    rec_index_file_name = file_dir + file_name + '/rec_index.txt'
    if (os.path.exists(rec_index_file_name)==True):
        rec_index = np.loadtxt(rec_index_file_name)
        print "Recording index: ", rec_index

    #Find star/end for multi-part recordings
    img_start = start_array[int(rec_index)-1]
    img_end = end_array[int(rec_index)-1] #temp_data[len(temp_data)-1][0]
    print "img_start: ", img_start
    print "img_end: ", img_end

   
    n_pixels = len(images_raw[0])
    n_frames = len(images_raw)

    #Compute correct img_start and img_end
    len_frame = float(img_end - img_start) / n_frames
    img_rate = 1. / len_frame
    img_rate_file = file_dir + file_name+ '/img_rate.txt'
    np.savetxt(img_rate_file, [img_rate])

    print "shifted img_start, img_end ", img_start, img_end
    
    img_times = []
    for i in range(n_frames):
        img_times.append(img_start+i*len_frame)
    img_times = np.array(img_times)
    
    return img_start, img_end, len_frame, img_rate, n_pixels, n_frames, img_times

def Spike_averages((args)):
    global images_temp
    temp3 = args
    return images_temp[temp3]

def Spike_averages_parallel_prespike_2sec((args)):
    global images_temp
    temp4, index = args
    sum_images = images_temp[temp4[0]]
    for i in range(1, len(temp4),1):
        temp_img = images_temp[temp4[i]] - np.mean(images_temp[temp4[i][0:len(temp4[i])/3]],axis=0) #Remove average of 1st 3rd of images
        sum_images+= temp_img
    return sum_images

def Spike_averages_parallel_prespike_3sec((args)):
    global images_temp, temp_n_pixels
      
    temp3 = args

    sum_images = np.zeros((len(temp3[0]),temp_n_pixels,temp_n_pixels), dtype=np.float32)
    for i in range(0,len(temp3),1):
        baseline = np.nanmean(images_temp[temp3[i][0:len(temp3[i])/2]], axis=0) 
        temp_img = images_temp[temp3[i]] #Remove avg of 1st half of images
        temp_frame = (temp_img - baseline)/baseline
        sum_images += temp_frame
    
    sum_images = sum_images/len(temp3)
    
    return sum_images
    

def Spike_averages_parallel_prespike_3sec_1D((args)):
    ''' Saves 64 x 64 size motifs for all spikes (i.e. not just averages) to be used for clustering 1D signals later'''

    global images_temp, temp_n_pixels

    temp3 = args

    #These numpy arrays store 3D motifs (time x height x width) as 2D motifs (height x (time x width)); TODO: Automate detection of the array sizes as these vary with imaging_rates
    #vectors = np.zeros((len(temp3), 64, 11392), dtype=np.float32) #Later recordings sampling rates have 178 frames per 6 seconds
    #vectors = np.zeros((len(temp3), 64, 11520), dtype=np.float32) #older recordings
    vectors = np.zeros((len(temp3), 64, 19072), dtype=np.float32) #50Hz recs
    
    for i in range(0, len(temp3), 1):
        baseline = np.average(images_temp[temp3[i][0:len(temp3[i])/2]], axis=0) 
        temp_img = images_temp[temp3[i]] #Remove avg of 1st half of images
        temp_frame = (temp_img - baseline)/baseline
        
        #sum_images += temp_frame
        #temp_stack = np.ma.hstack(temp_frame)
        temp_stack = np.hstack(temp_frame)
        indexes = np.isnan(temp_stack)
        temp_stack[indexes]=0

        #vectors[i] = scipy.misc.imresize(temp_stack,.25)
        vectors[i] = scipy.ndimage.interpolation.zoom(temp_stack,.25)

    
    #sum_images = sum_images/len(temp3)
    
    return vectors

def MCD_read_imagingtimes_old(MCDFilePath):
 
    #import necessary libraries
 
    import neuroshare as ns
    import numpy as np
 
    #open file using the neuroshare bindings
    fd = ns.File(MCDFilePath)
 
    #create index
    indx = 0
 
    #create empty dictionary
    data = dict()
 
    #loop through data and find all analog entities
 
    for entity in fd.list_entities():
        #print "looping over entities: ", entity

        if entity.entity_type == 1:
            data["extra"] = fd.entities[indx]
    return data


def Spike_averages_parallel_globalsignalregression_1D((args)):
    ''' Saves 64 x 64 size motifs for all spikes (i.e. not just averages) to be used for clustering 1D signals later'''
    global images_temp, temp_n_pixels

    temp3 = args

    vectors = []
    for i in range(len(temp3)): 
        temp_img = images_temp[temp3[i]] #select the image stack around     #temp3 already contains a stack of 178-180 frames centred on spike
        
        temp_stack = np.ma.hstack(temp_img)
        indexes_nan = np.isnan(temp_stack)
        temp_stack[indexes_nan]=0

        vectors.append(scipy.ndimage.interpolation.zoom(temp_stack,.25))
    
    return vectors

def Spike_averages_parallel((args)):
    global images_temp
    temp4, index = args
    sum_images = images_temp[temp4[0]]
    for i in range(1, len(temp4),1):
        sum_images+= images_temp[temp4[i]]
    return sum_images

def Sum_list((temp_list)):
    global temp_window, temp_img_rate, temp_n_pixels
    temp_sum = np.zeros((int(temp_window*temp_img_rate)*2, temp_n_pixels, temp_n_pixels), dtype=np.float16)
    for i in range(len(temp_list)):
        temp_sum += temp_list[i]
    return temp_sum
    
def Compute_spike_triggered_average(unit, channel, all_spikes, window, img_rate, img_times, n_pixels, images_aligned, file_dir, file_name, n_procs, overwrite, stm_types, random_flag, spiking_modes):
    '''Computes average frame values from t=-window .. +window (usually 180 to 270 frames) '''
    
    global images_temp, temp_window, temp_img_rate, temp_n_pixels
    
    temp_window = window
    temp_img_rate = img_rate
    temp_n_pixels = n_pixels

    print "No. of processors: ", n_procs

    #Remove spikes outside of imaging frames
    print "Total no. spikes: ", len(all_spikes)
    temp0 = np.where(np.logical_and(all_spikes>=img_times[0]+window, all_spikes<=img_times[-1]-window))[0]
    all_spikes = all_spikes[temp0]
    print "No. spikes within imaging period: ", len(all_spikes), " firing rate: ", float(len(all_spikes))/(img_times[-1]-img_times[0]), " Hz."

    #Save only spikes that are within imaging window:
    filename_spikes = file_dir+file_name+'/unit_'+str(unit).zfill(2)+ '_channel_' + str(channel).zfill(2) + '*.csv'
    file_out = glob.glob(filename_spikes)[0]
    np.savetxt(file_out[:-4]+"_imagingspikes.txt", all_spikes)       #Save only the spikes within imaging window

    #If no spikes in window return
    if len(all_spikes)==0:
        return images_aligned, spikes, stm_type

    #Loop over spiking modes
    for stm_type in stm_types: 
        spikes = all_spikes
        
        print "\n... processing stm_type: ", stm_type
        
        #Check to see if images already loaded and saved as .npy
        npy_file_name = file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+stm_type+'_'+str(window)+'sec_window_'+str(len(spikes)).zfill(5)+"_spikes"
        stack_file_name = file_dir + file_name + '/stack1D_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+stm_type+'_'+str(window)+'sec_window_'+str(len(spikes)).zfill(5)+"_spikes"

        #if (overwrite) or (os.path.exists(stack_file_name+'.npy')==False):
        if (overwrite) or (os.path.exists(npy_file_name+'.npy')==False):
            images_temp = np.array(images_aligned.copy(), dtype=np.float32)

            #Use allspikes
            if stm_type =='all':
                #Use all spikes
                print "# Spikes ", len(spikes)

                temp3 = []
                for i in spikes:        #temp3 contains all the frame indexes from img_times for each spike in raster; e.g. 180 frames for each spike automatically aligned
                    temp3.append(np.where(np.logical_and(img_times>=i-window, img_times<=i+window))[0][0:2*int(window*img_rate)]) #Fixed this value or could be off +/-1 frame
            
            #Use only spikes surrounded by 500ms silence on both sides...
            elif stm_type=='1sec':
                temp5 = []
                temp5.append(spikes[0])
                counter = 0
                for i in range(1,len(spikes)-1,1):
                    if ((spikes[i]-temp5[counter])>=.5) and ((spikes[i+1]-spikes[i])>=.5):
                        temp5.append(spikes[i])
                        counter+=1
                print "# Spikes isolated in 1sec windows: ", counter
                spikes = temp5

                temp3 = []
                for i in spikes:
                    temp3.append(np.where(np.logical_and(img_times>=i-window, img_times<=i+window))[0][0:2*int(window*img_rate)])

            #Use only spikes > 500ms apart
            elif stm_type=='burst':
                temp5 = []
                temp5.append(spikes[0])
                counter = 0
                for i in range(1,len(spikes),1):
                    if ((spikes[i]-temp5[counter])>=.5):
                        temp5.append(spikes[i])
                        counter+=1
                print "# Spikes beginning of bursts: ", counter
                spikes = temp5
                
                temp3 = []
                for i in spikes:
                    temp3.append(np.where(np.logical_and(img_times>=i-window, img_times<=i+window))[0][0:2*int(window*img_rate)])
                       
            
            #IF COMPUTING RANDOM/NOISE Flag
            if random_flag:
                print "... computing all motifs..."
                temp3 = []
                clipped_img_times = img_times[180:-90]  #For 30Hz rates make sure 3sec+3sec at beginning and 3sec at end
                for i in clipped_img_times:        #temp3 contains all the frame indexes from img_times for each spike in raster; e.g. 180 frames for each spike automatically aligned
                    temp3.append(np.where(np.logical_and(img_times>=i-window, img_times<=i+window))[0][0:2*int(window*img_rate)]) #Fixed this value or could be off +/-1 frame
            
                print "...# random motifs: ", len(temp3)
                
                np.savetxt(file_out[:-4]+"_imagingspikes.txt", clipped_img_times)       #Save only the spikes within imaging window
                #quit()

            
            #***************************************************************
            #GENERATE ONLY AVERAGE MOTIFS
            if True:
                
                #Multiple spiking modes
                if stm_type =='modes':      
                    mode_spikes = []
                    print file_out[:-4]+"_imagingspikes_grouped_spikes.txt"
                    
                    #Here loading the spikes for each mode; Each row contains the spikes for each mode in txt format to be human readable
                    if os.path.exists(file_out[:-4]+"_imagingspikes_grouped_spikes.txt")==False:
                        print "... spiking mode file '*_grouped_spikes.txt' missing..."
                        print "... can only generate 'all' spiking modes..."
                        print "See sta_maps.py code for details on how to generate mode spiking file..."
                        quit()
                        
                    with open(file_out[:-4]+"_imagingspikes_grouped_spikes.txt", 'rt') as inputfile:
                        reader = csv.reader(inputfile)
                        for row in reader:
                            mode_spikes.append(np.float32(row))
                            
                    #Load spiking modes from disk
                    for m in range(len(mode_spikes)): 
                        npy_file_name = file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+ \
                                        stm_type+'_'+str(window)+'sec_window_'+str(len(spikes)).zfill(5)+"_spikes_"+spiking_modes[m]
                        
                        if os.path.exists(npy_file_name+'.npy')==False: 
                                        
                            temp3 = []
                            for i in mode_spikes[m]:        #temp3 contains all the frame indexes from img_times for each spike in raster; e.g. 180 frames for each spike automatically aligned
                                temp3.append(np.where(np.logical_and(img_times>=i-window, img_times<=i+window))[0][0:2*int(window*img_rate)]) #Fixed this value or could be off +/-1 frame
                    
                            #Initialize list to capture data and store frame sequences;  Make MP Pool for spiking mode
                            print "... computing spiking mode: ", spiking_modes[m], "  #spikes: ", len(mode_spikes[m])
                            images_triggered_temp=[]
                            temp4 = []  #Contains only first spikes from bursts in cell
                            
                            pool = mp.Pool(n_procs)
                            chunks = int(len(mode_spikes[m])/n_procs) #Break up the temp3 array into n_procs that are "chunk" long each
                            for i in range(n_procs):
                                temp4.append(temp3[i*chunks: (i+1)*chunks])

                            print "Removing average of all pre spike frames - (time: -", window, "sec .. 0sec)"
                            images_triggered_temp.extend(pool.map(Spike_averages_parallel_prespike_3sec, temp4))
                            
                            pool.close()
                            print "... done "

                            #Computing averages spikes
                            print "Summing Number of chunks: ", len(images_triggered_temp)

                            temp_images = np.zeros((int(window*img_rate)*2, n_pixels, n_pixels), dtype=np.float16)
                            for i in range(len(images_triggered_temp)):
                                temp_images += images_triggered_temp[i]
                            
                            #DIVIDE BY NUMBER OF CHUNKS; Note used to be divided by number of spikes; also residue is being thrown out...
                            images_processed = temp_images/float(len(images_triggered_temp))

                            np.save(npy_file_name, images_processed)
                        
                        else:
                            images_processed = np.load(npy_file_name+'.npy')
                
                #Use all spikes for STMs.
                else: 
                    #Compute all frames based on image index;
                    print "Computing images from spike averages in parallel for window: ", window, " secs ..."
                    images_triggered_temp=[]
                    temp4 = []  #Contains only first spikes from bursts in cell

                    if len(spikes) < 30:
                        pool = mp.Pool(1)
                        temp4.append(temp3)            
                    else:
                        pool = mp.Pool(n_procs)
                        chunks = int(len(spikes)/n_procs) #Break up the temp3 array into n_procs that are "chunk" long each
                        for i in range(n_procs):
                            temp4.append(temp3[i*chunks: (i+1)*chunks])
                    
                    if True:
                        print "Removing average of all pre spike frames - (time: -", window, "sec .. 0sec)"
                        images_triggered_temp.extend(pool.map(Spike_averages_parallel_prespike_3sec, temp4))

                    if False:
                        indices = np.arange(len(spikes))
                        print "Removing average of 2sec - (time: -3sec .. -1sec)"
                        images_triggered_temp.extend(pool.map(Spike_averages_parallel_prespike_2sec, zip(temp4,indices)))
                    
                    pool.close()
                    print "... done "

                    #Sum over all spikes
                    print "Summing Number of chunks: ", len(images_triggered_temp)

                    temp_images = np.zeros((int(window*img_rate)*2, n_pixels, n_pixels), dtype=np.float16)
                    for i in range(len(images_triggered_temp)):
                        temp_images += images_triggered_temp[i]
                    
                    #Divide by the number of chunks
                    images_processed = temp_images/float(len(images_triggered_temp))

                    npy_file_name = file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+stm_type+'_'+str(window)+'sec_window_'+str(len(spikes)).zfill(5)+"_spikes"

                    np.save(npy_file_name, images_processed)


            #**********************************************************
            #SAVE ALL SINGLE SPIKE MOTIFS 
            if False:
                #SLIDING WINDOW METHOD
                if True: 
                    print "Removing average of all pre spike frames - (time: -", window, "sec .. 0sec)"
                    images_triggered_temp.extend(pool.map(Spike_averages_parallel_prespike_3sec_1D, temp4))
                    
                    print "... done... making 1D stack..."
                    stack_1D = vstack((images_triggered_temp))
                    print stack_1D.shape

                    np.save(stack_file_name, stack_1D)
                    pool.close()
        
                    return 0, 0, 0 #Dummy return variables
                    
                else:    #Global average
                    print "Removing average of all frames..."
                    baseline = np.mean(images_temp, axis=0)

                    vectors = []
                    images_temp = (images_temp - baseline)/baseline #Normalize all data to DF/F using GSR 
                    
                    images_triggered_temp.extend(pool.map(Spike_averages_parallel_globalsignalregression_1D, temp4))

                    print "... done... making 1D stack..."
                    stack_1D = vstack((images_triggered_temp))
                    print stack_1D.shape

                    np.save(stack_file_name+"_globalsignalregression", stack_1D)
                    pool.close()
        
                    return 0, 0, 0 #Dummy return variables
                  
        else: 
            print "Skipping processing of images ... loading from file"
            images_processed = np.load(npy_file_name+'.npy')
   
    return images_processed, spikes, stm_type

def view_static_stm(unit, main_dir, file_dir, file_name, stm_types, img_rate, spiking_modes):
    
    print "... viewing STMs..."
    
    #Select # of frames to average for display purposes only
    if img_rate < 32:       block = 10  #Number of frames to average for each frame plotted below
    else:                   block = 15  #~50Hz recordings have more frames, so average more per block
    
    midline_mask = 1        #Number of pixels to mask at the midline
    
    #Convert firing modes to a list of strings
    plot_strings=[]
    for mode_ in stm_types:
        if mode_=="modes":
            plot_strings.extend(spiking_modes)      #spiking modes are 'first', 'last', 'burst', 'tonic'
        else:
            plot_strings.extend([mode_])            #other spiking modes are 'all'; Others can be added

    fig = plt.figure()


    #Load STMs for each firing mode and plot
    for ctr, plot_string in enumerate(plot_strings):
        
        ax = plt.subplot(len(plot_strings),1,ctr+1)
        
        filenames = glob.glob(file_dir+file_name+'/img_avg_'+file_name+'_unit'+str(unit).zfill(2)+ "*"+plot_string+"*.npy")
        for fname in filenames:
            if 'maps' in fname: pass        #Hacky way of selecting the correct motif/montage map data (There are also single frame STM files containgin terms "minmaps"/'maxmaps')
            else: 
                filename= fname
                break
                
        data = np.load(filename)
        
        print filename

        img_stack = mask_data(data, main_dir, midline_mask)
        
        temp_stack = []
        for p in range(0, len(img_stack), block):
            temp_stack.append(np.mean(img_stack[p:p+block], axis=0))
        
        img_out = np.ma.hstack((temp_stack))
        v_max = np.nanmax(np.abs(img_out)); v_min = -v_max

        img_out[:, img_out.shape[1]/2-2:img_out.shape[1]/2+2] = -v_max
        

        im = plt.imshow(img_out, vmin = v_min, vmax=v_max)

        plt.title(plot_string+ '  DF/F: '+str(round(v_max*100,1))+"%", fontsize=15)
        
        if ctr == (len(plot_strings)-1): 
            new_xlabel = np.arange(-3.0, 3.1, 1)
            old_xlabel = np.linspace(0, img_out.shape[1], 7)
            plt.xticks(old_xlabel, new_xlabel, fontsize=25)
            
            plt.xlabel("Time from spike (sec)", fontsize = 30)
            ax.get_yaxis().set_visible(False)

        else: 
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

    
        #Colorbar
        #from mpl_toolkits.axes_grid1 import make_axes_locatable #Used for colorbar; allocates a bit of space for it

        #divider = make_axes_locatable(plt.gca())
        #cax = divider.append_axes("right", "2%", pad="1%")
        #cbar = plt.colorbar(im, cax=cax, ticks=[v_min, 0, v_max], orientation='vertical')
        
        #cbar.ax.set_yticklabels([str(round(v_min*100,1))+'%', '0%', str(round(v_max*100,1))+'%'])  # horizontal colorbar

        cbar = fig.colorbar(im, ticks = [v_min, 0, v_max], ax=ax, fraction=0.02, pad=0.05, aspect=3)
        cbar.ax.set_yticklabels([str(round(v_min*100,1))+"%", '0'+"%", str(round(v_max*100,1))+"%"])  # vertically oriented colorbar
        cbar.ax.tick_params(labelsize=15) 
        
    plt.suptitle("Experiment: " + file_name + "  unit: "+str(unit), fontsize=15)
    plt.show()

        
def remove_bregma(event):
    global bregma_coords, images_temp, n_pix
    
    if event.inaxes is not None:
        bregma_coords.append((event.ydata, event.xdata))
        for j in range(len(bregma_coords)):
            for k in range(n_pix):
                for l in range(7):
                    images_temp[100][k][min(n_pix-1,int(bregma_coords[j][1])-3+l)]=0
        
        plt.imshow(images_temp[100])
        fig.canvas.draw()
    else:
        plt.close()
        fig.canvas.mpl_disconnect(cid)


def define_area(event):
    
    global area_coords, images_temp, n_pix, win, img_r
    
    if event.inaxes is not None:
        area_coords.append((int(event.ydata), int(event.xdata)))
        #print int(event.ydata), int(event.xdata)
        for j in range(len(area_coords)):
            for k in range(3):
                for l in range(3):
                    images_temp[100][min(n_pix-1,int(area_coords[j][0])-1+k)][min(n_pix-1,int(area_coords[j][1])-1+l)]=0

        ax.imshow(images_temp[100])
        plt.show()

    else:
        #print 'Exiting'
        plt.close()
        fig.canvas.mpl_disconnect(cid)
        
def PointsInCircum(r,n=100):
    return [(math.cos(2*pi/n*x)*r,math.sin(2*pi/n*x)*r) for x in xrange(0,n+1)]
            
def on_click(event):
    
    global coords, images_temp, ax, fig, cid
    
    n_pix = len(images_temp[0])
    
    if event.inaxes is not None:
        coords.append((event.ydata, event.xdata))
        for j in range(len(coords)):
            for k in range(3):
                for l in range(3):
                    images_temp[100][min(n_pix,int(coords[j][0])-1+k)][min(n_pix,int(coords[j][1])-1+l)]=0

        ax.imshow(images_temp[100])
        #plt.show()
        fig.canvas.draw()
                    #figManager = plt.get_current_fig_manager()
                    #figManager.window.showMaximized()
    else:
        print 'Exiting'
        plt.close()
        fig.canvas.mpl_disconnect(cid)

def bregma_point(event):
    
    global coords, images_temp, ax, fig, cid, n_pix, win
    
    if event.inaxes is not None:
        coords.append((event.ydata, event.xdata))
        for j in range(len(coords)):
            for k in range(3):
                for l in range(3):
                    images_temp[100][min(n_pix,int(coords[j][0])-1+k)][min(n_pix,int(coords[j][1])-1+l)]=0

        #ax.imshow(images_temp[100], cmap = cm.Greys_r)

        #fig.canvas.draw()
        
        plt.close()
        #fig.close()
        fig.canvas.mpl_disconnect(cid)
        return

def Define_generic_mask(images_processed, main_dir):

    global coords, images_temp, ax, fig, cid
    
    images_temp = images_processed.copy()

    fig, ax = plt.subplots()

    if (os.path.exists(main_dir + 'genericmask.txt')==False):
        coords=[]

        ax.imshow(images_processed[100])#, vmin=0.0, vmax=0.02)
        ax.set_title("Compute generic (outside the brain) mask")
        #figManager = plt.get_current_fig_manager()
        #figManager.window.showMaximized()
        cid = fig.canvas.mpl_connect('button_press_event', on_click)
        plt.show()

        #******* MASK AND DISPLAY AREAS OUTSIDE GENERAL MASK 
        #Search points outside and black them out:
        all_points = []
        for i in range(len(images_processed[0][0])):
            for j in range(len(images_processed[0][0])):
                all_points.append([i,j])

        all_points = np.array(all_points)
        vertixes = np.array(coords) 
        vertixes_path = Path(vertixes)
        
        mask = vertixes_path.contains_points(all_points)
        counter=0
        coords_save=[]
        for i in range(len(images_processed[0][0])):
            for j in range(len(images_processed[0][0])):
                if mask[counter] == False:
                    images_processed[100][i][j]=0
                    coords_save.append([i,j])
                counter+=1

        fig, ax = plt.subplots()
        ax.imshow(images_processed[100])
        plt.show()
       
        genericmask_file = main_dir + 'genericmask.txt'
        np.savetxt(genericmask_file, coords_save)

        print "Finished Making General Mask"

    #else:
    #    print "Loading saved general mask"
        
        
    if (os.path.exists(main_dir + 'bregmamask.txt')==False):
        bregma_coords = []
        print "Making Bregma mask"
        ax.imshow(images_processed[100])#, vmin=0.0, vmax=0.02)
        ax.set_title("Compute bregma mask")
        #figManager = plt.get_current_fig_manager()
        #figManager.window.showMaximized()
        cid = fig.canvas.mpl_connect('button_press_event', remove_bregma)
        plt.show()

       
        bregmamask_file = main_dir + 'bregmamask.txt'
        np.savetxt(bregmamask_file, bregma_coords)

        print "Finished Bregma Mask"

    #else:
    #    print "Loading saved bregma mask"
        
    return generic_coords
    
def mouse_coords(event):
    
    global xycoords, fig, cid
    
    if event.inaxes is not None:
        print "Mouse coords: ", event.ydata, event.xdata
        xycoords = [event.ydata, event.xdata]
        plt.close()
        fig.canvas.mpl_disconnect(cid)
        
        #return coords
    #else:
    #    print 'Exiting'

def define_area_circle(event):
    
    global area_coords, fig, cid, circle_size
    
    if event.inaxes is not None:
        #Define N points on a circle centred at mouse click; shift circle to location
        points = np.vstack(PointsInCircum(circle_size,n=20))
        points = points + [int(event.ydata), int(event.xdata)]

        area_coords = []
        for i in range(len(points)):
            area_coords.append((points[i][0], points[i][1]))
            
        #Plot recent area
        temp_coords = np.array(area_coords)
        plt.plot(temp_coords[:,1],temp_coords[:,0],color='white',linewidth=4)
        
        fig.canvas.draw()
        
    else:
        plt.close()
        fig.canvas.mpl_disconnect(cid)


def define_area_rectangle(event):

    global area_coords, fig, cid, circle_size
    
    if event.inaxes is not None:
        #Define N points on a circle centred at mouse click; shift circle to location
        points = np.vstack(([0,0],[60,0], [60,12],[0,12], [0,0]))
        points = points + [int(event.ydata), int(event.xdata)]

        area_coords = []
        for i in range(len(points)):
            area_coords.append((points[i][0], points[i][1]))
            
        #Plot recent area
        temp_coords = np.array(area_coords)
        plt.plot(temp_coords[:,1],temp_coords[:,0],color='white',linewidth=4)
        
        fig.canvas.draw()
        
    else:
        plt.close()
        fig.canvas.mpl_disconnect(cid)


def define_area_circle_old(event):
    
    global area_coords, fig, cid, area_coords_left, area_coords_right
    
    if event.inaxes is not None:
        #Define N points on a circle centred at mouse click; shift circle to location
        points = np.vstack(PointsInCircum(3*(256/100),n=20))
        points = points + [int(event.ydata), int(event.xdata)]

        area_coords_left = []
        area_coords_right = []
        for i in range(len(points)):
            area_coords_left.append((points[i][0], points[i][1]))
            area_coords_right.append((points[i][0], 256-points[i][1]))
            
        #Plot recent area
        temp_coords = np.array(area_coords_left)
        plt.plot(temp_coords[:,1],temp_coords[:,0],color='white',linewidth=4)
        
        temp_coords = np.array(area_coords_right)
        plt.plot(temp_coords[:,1],temp_coords[:,0],color='white',linewidth=4)

        fig.canvas.draw()
        
    else:
        plt.close()
        fig.canvas.mpl_disconnect(cid)


def define_area_manual(event):
    
    global area_coords, images_temp, n_pix, win, img_r,fig
    
    if event.inaxes is not None:
        area_coords.append((int(event.ydata), int(event.xdata)))
        #print int(event.ydata), int(event.xdata)
        for j in range(len(area_coords)):
            for k in range(3):
                for l in range(3):
                    images_temp[int(area_coords[j][0])-1+k][int(area_coords[j][1])-1+l]=0
        
        ax.imshow(images_temp)
        fig.canvas.draw()

    else:
        plt.close()
        fig.canvas.mpl_disconnect(cid)


def Load_max_map(file_dir, area_name):
    
    global xycoords, fig, cid, n_pix
    n_pix = 256
    xycoords = []

    #Load normalized maps and stack them for display in proper order
    depths = ['cortex','subcortical']
    states = ['anesthetized', 'awake']
    maptypes = ['max', 'min']
    counter = 0
    for depth in depths:
        for state in states:
            for maptype in maptypes:
                temp = np.load(file_dir+maptype+'_maps_'+depth+'_'+state+'.npy')
                if counter>0:
                    img = np.vstack((img,temp))
                else:
                    img=temp
                counter+=1
    
    plt.close()
    #Display all maps and use mouseclick to select particular map
    fig, ax = plt.subplots()
    ax.imshow(img)
    ax.set_title(file_dir + "\nSelect map to search for: " + area_name, fontsize = 30)
    cid = fig.canvas.mpl_connect('button_press_event', mouse_coords)
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.show()

    #print "Click coords: ", xycoords
    #print "Size of image: ", 8, len(img[0])/n_pix
    #Select correct map from mouse click
    height = 8
    width = len(img[0])/n_pix
    for h in range(height):
        if (xycoords[0]>(h*n_pix)) and (xycoords[0]<((h+1)*n_pix)):
            y_panel = h
            break
    for w in range(width):
        if (xycoords[1]>(w*n_pix)) and (xycoords[1]<((w+1)*n_pix)):
            x_panel = w
            break

    counter = 0
    exit2=False
    for depth in depths:
        for state in states:
            for maptype in maptypes:
                if (counter == y_panel): 
                    exit2 = True
                    break
                counter+=1
            if exit2: break
        if exit2: break

    #Load correct row map and correct column
    images_temp = np.load(file_dir+maptype+'_maps_' + depth+'_'+state+'.npy')  
    max_maps = []
    for i in range(len(images_temp[0])/len(images_temp)):
        max_maps.append(images_temp[:,i*n_pix:(i+1)*n_pix])
    images_temp = max_maps[x_panel]
    
    return images_temp #Return only selected map for annotation


def Define_cortical_areas(file_dir, file_name, area_names, sides):
    print "Defining cortical areas"

    #Define maps from max-min super amax maps
    if True:
        global area_coords, ax, fig, cid, circle_size #Not sure need all these vars
        n_pixels = 256
      
        depth = 'cortex'#, 'subcortical']
        
        circle_sizes = [10,10,15,15,15,15,10,8]
        
        #for depth in depths:
        counter=0
        for area in area_names:
            circle_size = circle_sizes[counter]
            
            for side in sides:
                save_file = file_dir + depth+"_"+area+"_"+ side
                if (os.path.exists(save_file+'.npy')==False):
                    
                    #print "Which map to load for marking hindlimb?"
                    
                    images_temp = Load_max_map(file_dir, depth+ " " + area+" " + side)

                    area_coords = []
                    fig, ax = plt.subplots()
                    ax.imshow(images_temp)
                    ax.set_title(file_dir+"\nDefine Location of "+depth+" " + area+" "+side+' in ', fontsize=30)
                    #cid = fig.canvas.mpl_connect('button_press_event', define_area_manual)
                    
                    if counter==3:
                        cid = fig.canvas.mpl_connect('button_press_event', define_area_rectangle)
                    else:
                        cid = fig.canvas.mpl_connect('button_press_event', define_area_circle)
                    
                    #fig.canvas.update()
                    figManager = plt.get_current_fig_manager()
                    figManager.window.showMaximized()
                    plt.ylim(n_pixels,0)
                    plt.xlim(0,n_pixels)
                    plt.show()

                    #Convert coords into x and y data for plotting; select each column, then tack on first coordinate to complete circle
                    area_coords.append(area_coords[0])
                    area_coords = np.array(area_coords)
                    
                    from matplotlib import path
                    #Compute pixel locations inside cropped area
                    p = path.Path(area_coords)
                    all_pts = []
                    for i in range(256):
                        for j in range(256):
                            all_pts.append([i,j])
                    pts_inside = p.contains_points(all_pts)
                    
                    #Generate mask for saving 
                    mask_save = np.zeros((256,256),dtype=int8)+1
                    for i in range(256):
                        for j in range(256):
                            if pts_inside[i*256+j]==True:
                                mask_save[i,j]=False
                    
                    #Save mask
                    np.save(save_file, mask_save)
                    np.save(file_dir + "subcortical_"+area+"_"+ side, mask_save)
                    #Save contour
                    np.save(save_file+'_contour', area_coords)
                    np.save(file_dir + "subcortical_"+area+"_"+ side+'_contour', mask_save, area_coords)
                    
                    
            counter+=1

    else:
        print "... not defining areas.. loading stim based areas"



def Parallel_find_previous_frame_index((targets)):
    global din_times
    indexes = []
    for i in range(len(targets)):
        index = np.argmin(np.abs(din_times - targets[i]))
        
        if (targets[i] < din_times[index]):
            indexes.append(index-1)
        else:    
            indexes.append(index)
    return indexes


def remove_artifact(event):
    
    global coords, images_temp, n_pix, win, img_r
    
    if event.inaxes is not None:
        coords.append((event.ydata, event.xdata))
        for j in range(len(coords)):
            for k in range(7):
                for l in range(7):
                    xx = min(n_pix-1,int(coords[j][0])-1+k)
                    yy = min(n_pix-1,int(coords[j][1])-1+l)
                    images_temp[100][xx][yy]=0
                    if (np.array(coords) == [xx,yy]).all(-1).any(): #Check to see if pair already in list
                        pass
                    else:
                        coords.append([xx,yy])

        ax.imshow(images_temp[100])
        ax.set_title("Remove Artifacts")
        plt.show()

    else:
        print 'Exiting'
        plt.close()
        fig.canvas.mpl_disconnect(cid)

def load_tif(work_dir, file_name):
    
    if (os.path.exists(work_dir + file_name +'.npy')==False):
        print "Opening: ", work_dir + file_name+'.tif'
        img = Image.open(work_dir + file_name+'.tif')

        counter=0
        if True:
            while True:
                try:
                    img.seek(counter)
                except EOFError:
                    break
                counter+=1

        n_pixels = 128
        images_raw = np.zeros((counter, n_pixels, n_pixels), dtype = np.float32)

        for i in range(0, counter,1): 
            try:
                img.seek(i)
                images_raw [i] = img 
            except EOFError:
                break

        print "Saving imaging array..."

        #n_frames = min(counter,n_frames) #make sure you don't go past end of data; NB: technically it should be n_frames-offset_frames
        #images_raw = images_raw[offset_frames:offset_frames+n_frames]

        np.save(work_dir + file_name, images_raw)
        
    else:
        print "Loading .npy file from disk"
        images_raw = np.load(work_dir+file_name+'.npy')

    return images_raw



def filter_data():
    lowcut = 1
    highcut=70
    img_rate = 150        #Frame rate of imaging

    img_rate = len(images_raw)/pull_times[-1]
    window = int(2 * img_rate) #Number of seconds pre and post level pull
    
    data_array = np.zeros((window*2,128,128), dtype=np.float32)
    print "Computing pull triggered vids..."
    for trigger in frame_triggers:
        print "pulltime frame#: ", trigger
        data_chunk = images_raw[trigger-window:trigger+window]
        
        data_chunk = np.array(butter_bandpass_filter(data_chunk, lowcut, highcut, img_rate, order = 2))
        
        data_array+=data_chunk
    
    data = data_array/len(frame_triggers)
    

def remove_3sec_baseline(images_raw, frame_triggers, pull_times, work_dir, file_name):
    
    temp_fname = work_dir+file_name+"_3sec"
    if (os.path.exists(temp_fname+'.npy')==False):

        img_rate = len(images_raw)/pull_times[-1]
        window = int(2 * img_rate) #Number of seconds pre and post level pull
        
        data_array = np.zeros((window*2,128,128), dtype=np.float32)
        
        print "Computing pull triggered vids..."
        for trigger in frame_triggers:
            print "pulltime frame#: ", trigger
            if trigger <window: continue
            data_chunk = images_raw[trigger-window:trigger+window]
            baseline = np.average(data_chunk[0:window], axis=0)
            data_chunk = (data_chunk-baseline)/baseline
            
            #sigma_value = 1
            #data_array += ndimage.gaussian_filter(data_chunk, sigma=sigma_value) 
        
            data_array += data_chunk
            
        data = data_array/len(frame_triggers)
        
        #plt.imshow(data[0])
        #plt.show()
        
        np.save(temp_fname, data)
    
    else:
        data = np.load(temp_fname+'.npy')
        
    return data

def Define_bregma_lambda(images_rotated, main_dir, file_dir, file_name):
    
    global coords, images_temp, ax, fig, cid, n_pix, bregma_coords
    
    images_temp = images_rotated.copy()
    
    plt.close()
    n_pixels = 256
    n_pix = n_pixels
    
    if (os.path.exists(file_dir + 'bregma.txt')==False):

        coords = []
        fig, ax = plt.subplots()

        ax.imshow(images_temp[100], cmap = cm.Greys_r)
        ax.set_title("Define Bregma")
        plt.ylim(255,0)
        plt.xlim(0,255)
        for k in range(n_pixels/10):
            plt.plot([k*10,k*10],[1,n_pixels-2], color='red')
            plt.plot([1,n_pixels-2],[k*10,k*10], color='red')
        cid = fig.canvas.mpl_connect('button_press_event', bregma_point)
        plt.show()

        #Search points and black them out:
        images_temp[100][coords[0][0]][coords[0][1]]=1.0

        fig, ax = plt.subplots()
        ax.imshow(images_temp[100], cmap = cm.Greys_r)
        plt.show()
        
        bregma_file = file_dir + 'bregma.txt'
        np.savetxt(bregma_file, coords)
        #if '7-22' in file_dir:  np.savetxt(main_dir, coords)
           

    file_name = file_dir + file_name+'/'+file_name+'_images_aligned.npy'
    if (os.path.exists(file_name)==False):
        bregma_loc = np.loadtxt(file_dir+'/bregma.txt')
        #lambda_loc = np.loadtxt(file_dir+'/lambda.txt')
        print "exp bregma: ", bregma_loc 
        
        #Realign images:
        master_bregma = np.loadtxt(main_dir+'bregma.txt')
        print "master bregma: ", master_bregma
        
        y_shift = int(-master_bregma[0]+bregma_loc[0])
        x_shift = int(-master_bregma[1]+bregma_loc[1])
        
        #Shift image - x and y directions
        for k in range(len(images_temp)):
            temp_array = np.hstack((images_temp[k][:,x_shift:n_pixels], images_temp[k][:,0:x_shift]))
            images_temp[k] = np.vstack((temp_array[y_shift:n_pixels,:], temp_array[0:y_shift,:]))
        
        images_temp = np.float16(images_temp)
        np.save(file_name, images_temp)
        
        #plt.imshow(images_temp[100])
        #plt.title("Realigned images")
        #plt.show()
        
    else:
        images_temp = np.load(file_name)
    
        #plt.imshow(images_temp[100])
        #plt.title("Realigned images")
        #plt.show()
    
    return images_temp
   
    
def Define_artifact(images_processed, file_dir, file_name, window, n_pixels,img_rate):
    ''' Tool to manually ablate small regions of imaging area that may be 
        artifacts '''
    
    global coords, images_temp, ax, fig, cid, n_pix, win
    n_pix = n_pixels
    win = window

    print "Manual Mask Mode"
    images_temp = np.array(images_processed).copy()
    
    #Load Generic Mask
    if (os.path.exists(file_dir + 'genericmask.txt')==True):
        print "Loading existing generic mask"
        generic_mask_file = file_dir +'genericmask.txt'
        coords_generic = np.loadtxt(generic_mask_file)
        #Ablate generic map
        for i in range(len(images_temp)):
            for j in range(len(coords_generic)):
                images_temp[i][min(n_pixels-1,int(coords_generic[j][0]))][min(n_pixels-1,int(coords_generic[j][1]))]=0
                
    bregma_mask_file = file_dir + 'bregmamask.txt'
    if (os.path.exists(bregma_mask_file)==True):
        bregma_coords = np.loadtxt(bregma_mask_file)
        print "Loading bregma mask"
        #Remove centreline artifacts
        for i in range(len(images_temp)):
            for j in range(len(bregma_coords)):
                for k in range(n_pix):
                    for l in range(7):
                        images_temp[i][k][min(n_pix-1,int(bregma_coords[1])-3+l)]=0  #DON"T HARDWIRE MID SLICE

    #Load existing artiact  mask file
    coords=[]
    specific_mask_file = file_dir +'artifactmask.txt'
    if (os.path.exists(specific_mask_file)==True):
        temp_data= np.loadtxt(specific_mask_file)
        for i in range(len(temp_data)):
            coords.append(temp_data[i])
        #update_length=len(coords)

        #Ablate specific map
        for i in range(len(images_temp)):
            for j in range(len(coords)):
                for k in range(7):
                    for l in range(7):
                        images_temp[i][min(n_pixels-1,int(coords[j][0])-3+k)][min(n_pixels-1,int(coords[j][1])-3+l)]=0
    #else:
        #update_length=0
        
        fig, ax = plt.subplots()
        ax.imshow(images_temp[100])
        ax.set_title("Remove Artifacts")
        cid = fig.canvas.mpl_connect('button_press_event', remove_artifact)
        plt.show()
    
    #Save total map containing specific coords
    if len(coords)>0: np.savetxt(specific_mask_file, np.array(coords))


def Ablate_outside_area(n_pixels, contour_coords):
    ''' Function takes points from contour and returns list of coordinates
        lying outside of the contour - to be used to ablate image    '''
    
    #Search points outside and black them out:
    all_points = []
    for i in range(n_pixels):
        for j in range(n_pixels):
            all_points.append([i,j])

    all_points = np.array(all_points)
    vertixes = np.array(contour_coords) 
    vertixes_path = Path(vertixes)
    
    mask = vertixes_path.contains_points(all_points)
    ablate_coords=[]
    counter=0
    for i in range(n_pixels):
        for j in range(n_pixels):
            if mask[counter] == False:
                #images_processed[100][i][j]=0
                ablate_coords.append([i,j])
            counter+=1
            
    return ablate_coords

    
def Add_triple_arrays(coords1, coords2, coords3):
    
    coords = []
    if len(coords1)>0: coords.extend(coords1)
    if len(coords2)>0: coords.extend(coords2)
    if len(coords3)>0: coords.extend(coords3)
    
    return coords

def Add_double_arrays(bregma_coords_temp, generic_coords):
    coords = np.array(np.vstack((bregma_coords_temp, generic_coords)), dtype=np.int16)
    
    return coords

def Load_areas_and_mask(depth, unit, channel, n_pixels, main_dir, file_dir, file_name, images_aligned, area_names, sides):
    
    #Create set of images, saved coordinates and borders for output
    images_areas = []   #Make list of lists to hold images for each area

    with open(file_dir+file_name+'/depth.txt', "r") as f:
        data = csv.reader(f)
        for row in data:
            depth = row[0]
    
    #Load stimulus evoked areas:
    print "Loading ROIs..."
    for area in area_names:
        for side in sides:
            #print area+" "+side
            area_file = file_dir+ depth+'_' + area+'_'+side
            #print area_file
            if (os.path.exists(area_file+'.npy')==True):
                area_mask = np.load(area_file+'.npy')     #Load coordinates; 2D vectors stacked; Mask out everything else
               
                a = []
                for i in range(len(images_aligned)): #Mask every processed frame in +/- window period
                    a.append(np.ma.array(images_aligned[i], mask=area_mask, fill_value = 0., hard_mask = True))
                                                    
                images_areas.append(a)
    
    return images_areas



def Compute_static_maps_max(img_rate, window, n_procs, main_dir, file_dir, file_name, n_pixels, unit, channel, n_spikes, ptp, stm_types, plotting, spiking_modes):
    ''' Average pre- and post- spike intervals to obtain static map; See also the "max" version of this function'''
    
    print "Computing static maps"
    
    plot_strings=[]
    stm_filenames=[]

    for mode_ in stm_types:
        if mode_=="modes":
            plot_strings.extend(spiking_modes)    #Need these new modes for the revision analysis
            
            for spiking_mode in spiking_modes:
                print file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_modes_*'+spiking_mode+'.npy'
                stm_filenames.append(glob.glob(file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_modes_*'+spiking_mode+'.npy')[0])
        else:
            plot_strings.extend(mode_)
            stm_filenames.append(glob.glob(file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_all_*.npy')[0])
        
    for stm_filename in stm_filenames:   #Loop over spiking modes

        if os.path.exists(stm_filename)==False: 
            print "... cell STM not computed...exiting..."
            return 
        
        images_areas = np.load(stm_filename)

        #Load General mask (removes background)
        generic_mask_file = []
        generic_mask_file = main_dir + 'genericmask.txt'
        if (os.path.exists(generic_mask_file)==True):
            generic_coords = np.int16(np.loadtxt(generic_mask_file))
        
        generic_mask_indexes=np.zeros((n_pixels,n_pixels))
        for i in range(len(generic_coords)):
            generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True

        #Compute and Save Maxmaps
        temp_array = []
        for i in range(-int(img_rate), int(img_rate),1):
            temp = np.ma.array(images_areas[int(window*img_rate)+i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
            temp_array.append(temp)

        temp_array = np.float32(temp_array)
        images_out2 = np.amax(temp_array,axis=0)
        images_out2 = np.ma.array(images_out2, mask=generic_mask_indexes, fill_value = 0)
        images_out2 = np.ma.filled(images_out2, 0.0)
        images_out2 = np.nan_to_num(images_out2)

        np.save(stm_filename[:-4]+"_maxmaps", images_out2)

        #SAVE MIN MAPS ALSO
        if True:
            temp_array = []
            for i in range(-int(img_rate), int(img_rate),1):
                temp = np.ma.array(images_areas[int(window*img_rate)+i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
                temp_array.append(temp)

            temp_array = np.float32(temp_array)
            images_out2 = np.amin(temp_array,axis=0)
            images_out2 = np.ma.array(images_out2, mask=generic_mask_indexes, fill_value = 0.)
            
            images_out2 = np.ma.filled(images_out2, np.min(images_out2))
            
            np.save(stm_filename[:-4]+"_minmaps", images_out2)

        if plotting: 
            
            #Display max pixel maps
            midline_mask = 3        #Number of pixels to mask at the midline

            fig, ax = plt.subplots()
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

            image_loaded = np.load(stm_filename[:-4]+"_maxmaps.npy")
            image_loaded = mask_data_single_frame(image_loaded, main_dir, midline_mask)

            #Process image loaded
            v_max = np.nanmax(np.abs(image_loaded)); v_min = -v_max
            print "...v_max: ", v_max
            im = plt.imshow(image_loaded, vmin=v_min, vmax=v_max)

            plt.title("STM:    " + stm_filename, fontsize=20)

            cbar = fig.colorbar(im, ticks = [v_min, 0, v_max], ax=ax, shrink=.5, pad=.1, aspect=5)
            cbar.ax.set_yticklabels([str(round(v_min*100,1))+"%", '0'+"%", str(round(v_max*100,1))+"%"])  # vertically oriented colorbar
            cbar.ax.tick_params(labelsize=20) 

            plt.show()
    
def Load_max_maps(ptps, file_dir, file_name):

    with open(file_dir+file_name+'/depth.txt', "r") as f:
        data = csv.reader(f)
        for row in data:
            depth = row[0]
    with open(file_dir+file_name+'/state.txt', "r") as f:
        data = csv.reader(f)
        for row in data:
            state = row[0]

    max_maps = []
    min_maps = []
    max_maps_nspikes = []
    min_maps_nspikes = []
    max_maps_ptp = []
    min_maps_ptp = []
    max_maps_depth = []
    min_maps_depth = []
    max_maps_state = []
    min_maps_state = []
    max_maps_channel = []
    min_maps_channel = []
    max_maps_unit = []
    counter = 0
    files = os.listdir(file_dir+file_name+'/')
    for file_ in files: #Load all cells from list of experiments
        #print file_
        if ('maxmap' in file_):
        #if ('avemap' in file_):

            n_spikes = int(file_[-9:-4])
            ptp = int(file_[-20:-17])
            channel = int(file_[-27:-25])
            unit = int(file_[-38:-36])
            #if qc[unit]==0: continue
            
            maps = np.load(file_dir+file_name+'/'+file_)
            max_maps.append(maps)
            max_maps_nspikes.append(n_spikes)
            max_maps_ptp.append(ptp)
            max_maps_depth.append(depth)
            max_maps_state.append(state)
            max_maps_channel.append(channel)
            max_maps_unit.append(unit)
            #print channel
            #print "Maxmap spikes: ", n_spikes
            
        if ('minmap' in file_):
            n_spikes = int(file_[-9:-4])
            ptp = int(file_[-20:-17])
            channel = int(file_[-27:-25])

            maps = np.load(file_dir+file_name+'/'+file_)
            min_maps.append(maps)
            min_maps_nspikes.append(n_spikes)
            min_maps_ptp.append(ptp)
            min_maps_depth.append(depth)
            min_maps_state.append(state)
            min_maps_channel.append(channel)
            #print "Minmap spikes: ", n_spikes
    #max_maps = max_maps/counter/255
    
    return max_maps, min_maps, max_maps_nspikes, min_maps_nspikes, max_maps_ptp, min_maps_ptp, max_maps_depth,min_maps_depth,max_maps_state,min_maps_state,max_maps_channel,min_maps_channel,max_maps_unit


def Mp_max((temp_array1)):
    global coords_temp
    
    for j in range(len(coords_temp)):
        temp_array1[coords_temp[j][0]][coords_temp[j][1]] = -1000

    return temp_array1

    
def Mp_min((temp_array2)):
    global coords_temp

    for j in range(len(coords_temp)):
        temp_array2[coords_temp[j][0]][coords_temp[j][1]] = 1000

    return temp_array2

def Zero_images((images_areas)):
    global coords_temp
    
    for i in range(len(images_areas)):
        print "Processing time: ", i
        for j in range(len(coords_temp)):
            print "Processing coordinate: ", j
            images_areas[i][min(255,int(coords_temp[j][0]))][min(255,int(coords_temp[j][1]))]=0
            
    return images_areas

def Search_max_min(unit, channel, spikes, file_dir, file_name, img_rate, window, n_procs, area_names, depth, sides, stm_types, spiking_modes):

    global coords_temp

    #Convert firing modes to a list of strings
    plot_strings=[]
    stm_filenames=[]
    for mode_ in stm_types:
        if mode_=="modes":
            plot_strings.extend(spiking_modes)    #Need these new modes for the revision analysis
            for spiking_mode in spiking_modes:
                stm_filenames.append(glob.glob(file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_modes_*'+spiking_mode+'.npy')[0])
        else:
            plot_strings.extend(mode_)
            stm_filenames.append(glob.glob(file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_all_*.npy')[0])
                    
                        
    for plot_string, stm_filename in zip(plot_strings, stm_filenames):
    
        images_processed = np.load(stm_filename)
        n_pixels = images_processed.shape[1]
        
        Max_plot = []
        Min_plot = []
        Max_index = []
        Min_index = []
        Max_pixel_value = []
        Min_pixel_value = []
        Images_areas = []

        #Search max/min for each ROI defined area
        counter=0
        for area in area_names:
            for side in sides:
                print "Searching max/min for area # ", area, " ", side
                
                area_file = file_dir+ '/roi/'+depth+'_' + area+'_'+side

                if (os.path.exists(area_file+'.npy')==True):
                    area_mask = np.load(area_file+'.npy')     #Load coordinates; 2D vectors stacked; Mask out everything else
                else:
                    print "...file not found: ", area_file+'.npy'
                
                coords_temp = area_mask

                temp_array = np.ma.array(np.zeros(images_processed.shape), mask=True)
                for i in range(len(images_processed)):
                    temp_array[i] = np.ma.array(images_processed[i], mask=coords_temp, fill_value = 0., hard_mask = True)
                
                temp_max_array = []
                temp_min_array = []
                temp_array1 = np.ma.array(np.zeros((int(img_rate*2),n_pixels,n_pixels)), mask=True)
                temp_array2 = np.ma.array(np.zeros((int(img_rate*2),n_pixels,n_pixels)), mask=True)
                #temp_array1 = np.zeros((int(img_rate*2),n_pixels,n_pixels), dtype=np.float32)
                #temp_array2 = np.zeros((int(img_rate*2),n_pixels,n_pixels), dtype=np.float32)
                for i in range(int(img_rate*2)):     #This searches +/- 1 sec from time = 0 sec; OR Maximum to begining of data;
                    temp_array1[i] = np.ma.array(temp_array[int(window*img_rate)-int(img_rate)+i]).copy()
                    temp_array2[i] = np.ma.array(temp_array[int(window*img_rate)-int(img_rate)+i]).copy()
                    #temp_array1[i] = temp_array[int(window*img_rate)-int(img_rate)+i].copy()
                    #temp_array2[i] = temp_array[int(window*img_rate)-int(img_rate)+i].copy()

                temp_array1._sharedmask=False
                temp_array2._sharedmask=False


                #Set masked background values to very large/very low values to not confuse search for max/min values
                pool = mp.Pool(n_procs)
                temp_max_array.extend(pool.map(Mp_max, temp_array1))
                pool.close()

                pool = mp.Pool(n_procs)
                temp_min_array.extend(pool.map(Mp_min, temp_array2))
                pool.close()

                #Search for max value using 1D unravel of temp_arrays; assign location of max/min index; detect max/min values overall;
                temp1_array = np.ma.array(np.zeros((int(img_rate*2),n_pixels,n_pixels)), mask=True)
                for k in range(len(temp_max_array)):
                    temp1_array[k] = np.ma.array(temp_max_array[k], mask=area_mask) #np.ma.array(temp_max_array)

                temp_max_array = temp_array1
                max_index = np.unravel_index(np.nanargmax(temp_max_array), temp_max_array.shape)
                max_pixel_value = temp_array[int(window*img_rate)-int(img_rate)+max_index[0]][max_index[1]][max_index[2]]

                #temp1_array = np.ma.array(np.zeros((img_rate*2,256,256)), mask=True)
                #for k in range(len(temp_max_array)):
                    #temp1_array[k] = np.ma.array(temp_min_array[k], mask=area_mask) #np.ma.array(temp_max_array)
                temp_min_array = temp_array2
                min_index = np.unravel_index(np.nanargmin(temp_min_array), temp_min_array.shape)
                min_pixel_value = temp_array[int(window*img_rate)-int(img_rate)+min_index[0]][min_index[1]][min_index[2]]
                
                #max/min_index at index [0] contain time slice info; needed for max/min_pixel_value above, but not after;
                max_index = max_index[1:] #Reduce back to 2D arrays
                min_index = min_index[1:]
               
                max_plot = []
                min_plot = []
                
                for i in range(len(temp_array)):
                    max_plot.append(temp_array[i][max_index[0]][max_index[1]])
                    min_plot.append(temp_array[i][min_index[0]][min_index[1]])
                  
                Max_plot.append(max_plot)
                Min_plot.append(min_plot)
                Max_index.append(max_index)
                Min_index.append(min_index)
                Max_pixel_value.append(max_pixel_value)
                Min_pixel_value.append(min_pixel_value)
                Images_areas.append(temp_array)
                
                counter+=1


        #Make save array; 
        temp_array = []
        temp_array.append([len(spikes)])  #Save # spikes within imaging period
        
        #Add all time courses to 
        print "...Max_plot: ", Max_plot
        print Max_plot[0]
        for i in range(len(Max_plot)):
            max_plot=np.array(Max_plot[i])
            min_plot=np.array(Min_plot[i])

            temp_array.append(min_plot)
            temp_array.append(max_plot)


        #Add max_index and min_index information to the end of the time_course_data*.txt file
        for i in range(len(Max_index)):
            temp_array.append(Max_index[i])
            temp_array.append(Min_index[i])


        #Save time_course_data
        with open(file_dir + file_name + '/time_course_data_' + file_name+'_'+plot_string+ '_'+str(window)+'sec_window_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+"_"+area_names[0], "w") as f:
            writer = csv.writer(f)
            writer.writerows(temp_array)
            
            
        #Add max_index and min_index information to the end of the time_course_data*.txt file
        for i in range(len(Max_index)):
            temp_array.append(Max_index[i])
            temp_array.append(Min_index[i])


        #Save time_course_data
        with open(file_dir + file_name + '/time_course_data_' + file_name+'_'+plot_string+ '_'+str(window)+'sec_window_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+"_"+area_names[0], "w") as f:
            writer = csv.writer(f)
            writer.writerows(temp_array)
        

    return Max_plot, Min_plot, Max_pixel_value, Min_pixel_value, Max_index, Min_index, area_names, Images_areas

def Average_roi(images_areas, img_rate, window, n_procs, area_names, sides):
    
    global coords_temp

    ave_plot = []
    print "Averaging areas..."
    
    #Search max/min for each ROI defined area
    counter=0
    for area in area_names:
        for side in sides:
            
            max_plot = []
            for i in range(len(images_areas[counter])):     #This searches +/- 1 sec from time = 0 sec; OR Maximum to begining of data;
                value = np.nanmean(images_areas[counter][i])
                max_plot.append(value)
            
            ave_plot.append(max_plot)

            counter+=1

    return ave_plot
    


def Save_time_course(unit, channel, spikes, Max_plot, Min_plot, Max_index, Min_index, window, len_frame, file_dir, file_name, area_names, sides, stm_types):
    
    
    temp_array = []
    temp_array.append([len(spikes)])  #Save # spikes within imaging period
    temp_array2 = []
    for area in area_names:
        for side in sides:
            temp_array2.append(area+"_"+side)   #Save names of areas recorded
    
    for i in range(len(Max_plot)):
        ax = plt.subplot(8,2,i+1)
        max_plot=np.array(Max_plot[i])
        min_plot=np.array(Min_plot[i])

        temp_array.append(min_plot)
        temp_array.append(max_plot)

        xx = np.arange(-window, window, len_frame)
        xx = xx[0: len(max_plot)]

        ax.plot(xx, max_plot, color='black', linewidth=2)
        ax.plot(xx, min_plot, color='blue', linewidth=2)
        ax.plot([-3,3],[0,0], color='black')
        ax.plot([0,0],[min(min_plot),max(max_plot)], color='black')
        ax.set_ylim(-.04, .041)
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.yaxis.set_ticks(np.arange(-0.04,0.041,0.04))
        
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        if i==0: 
            plt.title("Left")
            ax.get_yaxis().set_visible(True)
            plt.ylabel(area_names[int(i/2)], fontsize=8)
        elif (i%2)==0 and i<14:
            ax.get_yaxis().set_visible(True)
            plt.ylabel(area_names[int(i/2)], fontsize=8)
        elif i==1:
            #ax.get_xaxis().set_visible(True)
            plt.title("Right")
        elif i==13:
            ax.get_xaxis().set_visible(True)


    plt.suptitle(file_name + " Unit: " +str(unit).zfill(2) + " Channel: " + str(channel).zfill(2)+
    " No. of spikes: "+ str(len(spikes)))

    #npy_file_name = file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+plot_string+'_'+str(window)+'sec_window'

    plt.savefig(file_dir + file_name + '/time_course_plot_' + file_name+'_'+plot_string+ '_'+str(window)+'sec_window_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'.png', fontsize = 20)
    plt.close()

    #Add max_index and min_index information to the end of the time_course_data*.txt file
    for i in range(len(Max_index)):
        temp_array.append(Max_index[i])
        temp_array.append(Min_index[i])

    #Save time_course_data
    with open(file_dir + file_name + '/time_course_data_' + file_name+'_'+plot_string+ '_'+str(window)+'sec_window_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+"_"+area_names[0], "w") as f:
        writer = csv.writer(f)
        writer.writerows(temp_array)


def Plot_matrix_maps(average_areas, file_dir, file_name, area_names, img_rate, unit, spikes, channel, ptp):
    
    #fig = plt.figure()
    plt.close()
    gs = gridspec.GridSpec(2,6)
    
    fig, ax = plt.subplots(nrows=1, ncols=2)

    #Compute time courses
    
    img1 = []
    for k in range(len(area_names)):
       img1.append(np.float16(average_areas[k*2])*1E2)
    v_abs1 = max(np.max(np.array(img1)), -np.min(np.array(img1)))

    img2 = []
    for k in range(len(area_names)):
       img2.append(np.float32(average_areas[k*2+1])*1E2)
    v_abs2 = max(np.max(np.array(img2)), - np.min(np.array(img2)))

    global_max = max(v_abs1,v_abs2)

    #****** Plot left hemisphere
    #ax = plt.subplot(gs[0,0:1])
    ax1 = plt.subplot(121)

    im  = ax1.imshow(img1, aspect='auto', cmap=plt.get_cmap('jet'), vmin=-global_max, vmax=global_max, interpolation='none')
    
    yy = np.arange(0,len(area_names),1)
    #labels = [ 'hindlimb', 'forelimb', 'barrel', 'motor', 'visual', 'retrosplenial', 'acc']
    plt.yticks(yy, area_names) 
    plt.tick_params(axis='both', which='both', labelsize=8)

    original_xx = np.arange(0,len(average_areas[0])+2,img_rate)
    xx = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
    plt.ylim(-0.5,-0.5+len(area_names))
    plt.xlim(0,len(average_areas[0]))
    plt.xticks(original_xx, xx)
    
    plt.plot([3*img_rate,3*img_rate],[-0.5,-0.5+len(area_names)], color='black', linewidth=2)
    for i in range(len(area_names)):
        plt.plot([0,len(average_areas[0])+2],[-0.5+i,-0.5+i], 'r--', color='black', linewidth=1)
    for i in range(7):
        plt.plot([0+i*img_rate,0+i*img_rate],[-0.5,-0.5+len(area_names)], 'r--', color='black', linewidth=1)
    plt.title("Left (DF/F max: "+str(round(np.max(img1),2))+"  min: "+str(round(np.min(img1),2))+")", fontsize=10) 

    np.save(file_dir+file_name+'/'+file_name+'_matrix_unit_'+str(unit).zfill(2)+'_left', img1)
    
    
    #****** Plot right hemisphere
    #ax = plt.subplot(gs[0,1:2])
    ax2 = plt.subplot(122)
    im = ax2.imshow(img2, aspect='auto', cmap=plt.get_cmap('jet'), vmin=-global_max, vmax=global_max, interpolation='none')
    plt.tick_params(axis='both', which='both', labelsize=8)
    ax2.get_yaxis().set_visible(False)

    plt.ylim(-0.5,-0.5+len(area_names))
    plt.xlim(0,len(average_areas[0]))
    plt.xticks(original_xx, xx)
    plt.plot([3*img_rate,3*img_rate],[-0.5,-0.5+len(area_names)], color='black', linewidth=2)
    for i in range(len(area_names)):
        plt.plot([0,len(average_areas[0])+2],[-0.5+i,-0.5+i], 'r--', color='black', linewidth=1)
    for i in range(7):
        plt.plot([0+i*img_rate,0+i*img_rate],[-0.5,-0.5+len(area_names)], 'r--', color='black', linewidth=1)

    plt.title("Right (DF/F max: "+str(round(np.max(img2),2))+"  min: "+str(round(np.min(img2),2))+")", fontsize=10) 


    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    plt.suptitle(file_name + "  unit: " + str(unit)+ "  ptp: " + str(ptp)+ "uV  spikes: "+str(len(spikes))+ "  ch: "+str(channel))

    np.save(file_dir+file_name+'/'+file_name+'_matrix_unit_'+str(unit).zfill(2)+'_right', img2)


    if (ptp>35) and (len(spikes)>128) and (global_max>0.5):
        plt.savefig(file_dir + file_name + '/'+file_name+'_matrix_unit_'+str(unit).zfill(2)+'.png', dpi=100)#, fontsize = 20)
        plt.close()
    else:
        plt.close()
    


def Animate_images(unit, channel, window, img_rate, Images_areas, file_dir, file_name, n_pixels, spikes, plot_string, n_procs, generic_mask_indexes, 
    Max_plot, Min_plot, Max_pixel_value, Min_pixel_value, Max_index, Min_index, area_names, sides, depth):
    
    print "Generating Animation. No. of Frames: ", len(Images_areas[0])
    
    colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','lightsalmon','pink','darkolivegreen','blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','lightsalmon','pink','darkolivegreen']

    #***************** Process plots *****************
    plot_frames = []
    x = np.arange(-window, window, 1./img_rate)
    x = x[:len(Max_plot[0])]
    Max_plot = np.array(Max_plot)

    print "Number of area recordings: ", len(Max_plot)
    ranges = np.arange(0,len(x),1)

    #Divide data into nchunks for some reason works faster; TO DO: PARALLELIZE
    n_chunks = 8
    range_chunks = np.array_split(ranges,n_chunks)
    
    plt.close()
    
    #*********** PLOT CURVE GRAPHS ************************
    if True: #PLOT ALL CURVES FOR ALL AREAS
        print "... processing plot frames..."
        for chunk in range(n_chunks):
            fig = plt.figure()
            fig.add_subplot(111)
            fig.set_size_inches(5, 20)
            #fig.canvas.draw()
            
            #FIXED LINES FOR EACH PLOT
            for i in range(len(Max_plot)+2): #need to make figs at top of graph also                
                plt.plot(x, [int(i/2)*0.10]*len(x), color='black', linewidth = 2)
                plt.plot(x, [int(i/2)*0.10+0.03]*len(x), 'r--', color=colors[int(i/2)], alpha=0.5)
                plt.plot(x, [int(i/2)*0.10-0.03]*len(x), 'r--', color=colors[int(i/2)], alpha=0.5)
                plt.axhspan(int(i/2)*0.10+0.03, int(i/2)*0.10-0.03, color='black', alpha=0.2, lw=0)   
            plt.xlabel("Time (s)")
            plt.ylim(-0.06, 0.90)
            plt.plot([-window,window],[-100,100], color='black')

            plt.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='off') # 

            #DRAW CURVES OVER TIME
            for k in range_chunks[chunk]:
                for i in range(len(Max_plot)):

                    plt.plot(x[:k], Max_plot[i][:k]+(int(i/2)+1)*0.10, color=colors[i], linewidth=2) #Plot individual overlayed left/rigth curves
                    
                fig.savefig(file_dir + file_name+'/figs_'+str(k).zfill(3)+'.jpg', dpi=40)
                #fig.canvas.draw()
                #fig.canvas.flush_events()
                #data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
                #data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
                #scipy.misc.toimage(data).save(file_dir + file_name+'/figs_'+str(k).zfill(3)+'.jpg')
            plt.close(fig)

        #Make vid1 containing time course curves
        devnull = open(os.devnull, 'w')
        subprocess.call("ffmpeg -f image2 -r 15 -i " + file_dir+file_name+ "/figs_%03d.jpg  -vcodec libx264 -y -vf scale=480:960 "+file_dir+file_name+'/vid1.mp4', shell=True,stdout=devnull, stderr=devnull)
        os.system('rm '+file_dir+file_name+'/figs*')


    #***************** Process images ****************
    Images_areas = np.array(Images_areas)
    Images_areas = np.swapaxes(Images_areas, 0, 1)
    
    #Generate borders around recorded areas and mask them out in order for color to work properly
    img=Image.new("RGBA", (n_pixels,n_pixels),(0,0,0))
    if True:
        if True:
            draw = ImageDraw.Draw(img)
            counter=0
            for area in area_names:
                for side in sides:
                    area_file = file_dir+ depth+'_' + area+'_'+side+'_contour'
                    if (os.path.exists(area_file+'.npy')==True):
                        area_contour = np.load(area_file+'.npy')     #Load coordinates; 2D vectors stacked; Mask out everything else

                    for i in range(len(area_contour)-1):  #Skip last border which represents all cortex
                        draw.line((area_contour[i][1],area_contour[i][0],area_contour[i+1][1],area_contour[i+1][0]), fill = colors[counter], width = 1)
                counter+=1
        borders_=asarray(img)

        if True:
            #Generate borders mask to zero out [Ca] image pixels to get color right
            borders_mask = np.zeros((n_pixels,n_pixels),dtype=bool)
            for i in range(n_pixels):
                for j in range(n_pixels):
                    if np.sum(borders_[i][j])>255:
                        borders_mask[i][j]=True

        #if False:
            #draw = ImageDraw.Draw(img)
            #for area in area_names:
                #for side in sides:
                    #area_file = file_dir+ depth+'_' + area+'_'+side
                    #if (os.path.exists(area_file+'.npy')==True):
                        #area_mask = np.load(area_file+'.npy')     #Load coordinates; 2D vectors stacked; Mask out everything else

                    ##for i in range(len(area_mask)):  #Skip last border which represents all cortex
                    ##    for j in range(len(area_mask[[0]])):
                    ##        if area_mask[i][j]==True:
                    ##            print i,j
                    ##            draw.point((i,j), fill = 'white')

                
    #Generate time index text and convert it to array for adding to frames
    font = ImageFont.truetype("/usr/share/fonts/truetype/ttf-dejavu/DejaVuSansMono.ttf",9)
    time_text=[]
    for i in range(len(Images_areas)):
        img=Image.new("RGBA", (n_pixels,n_pixels),(0,0,0))
        draw = ImageDraw.Draw(img)
        #if n_pixels>128:
        #    draw.text((0, n_pixels-10),"Spikes: " + str(len(spikes))+ " Time: "+str(format((float(i)/len(Images_areas)-.5)*6.0,'.2f')),(255,255,255),font=font)
        #else:
        draw.text((10, 10),"Spikes: " + str(len(spikes))+ " Time: "+str(format((float(i)/len(Images_areas)-.5)*6.0,'.2f')),(255,255,255),font=font)
        draw = ImageDraw.Draw(img)
        time_text.append(asarray(img))

    label_text=[]
    for i in range(len(Images_areas)):
        img=Image.new("RGBA", (n_pixels,n_pixels),(0,0,0))
        draw = ImageDraw.Draw(img)
        #if n_pixels>128:
        #    draw.text((0, n_pixels-20),file_name + " Unit: "+str(unit),(255,255,255),font=font)
        #else:
        draw.text((10, 0),file_name + " Unit: "+str(unit),(255,255,255),font=font)
        draw = ImageDraw.Draw(img)
        label_text.append(asarray(img))
    
    #Generate image frames as arrays and save to .pngs
    images_out = []
    my_cmap = matplotlib.cm.get_cmap('jet')
    v_min = min(Min_pixel_value)
    v_max = max(Max_pixel_value)

    print "... processing images frames..."
    for i in range(len(Images_areas)):
        temp_img = np.float32(Images_areas[i][0])
        temp_img = ndimage.gaussian_filter(temp_img, sigma=.5)
        temp_img = (temp_img - v_min)/(v_max-v_min)
        masked_data = np.ma.array(temp_img, mask=generic_mask_indexes)
        
        if True: masked_data = np.ma.array(masked_data, mask=borders_mask)
        
        images_out = my_cmap(masked_data, bytes=True) #, vmin=min(Min_pixel_v), vmax=max(Max_pixel_v), origin='lower') #cmap=my_cmap, clim=[0.9, 1]) #, cmap=cmap, interpolation='nearest')
        
        if True: images_out = borders_+ images_out + time_text[i] + label_text[i]

        scipy.misc.imsave(file_dir + file_name+'/figs_'+str(i).zfill(3)+'.png', images_out)
        #scipy.misc.toimage(images_out, cmin=v_min, cmax=v_max).save(file_dir + file_name+'/figs_'+str(i).zfill(3)+'.png')
   
    print "... making vids ..."

    #Make video 2 - [Ca] imaging vid
    devnull = open(os.devnull, 'w')
    subprocess.call("ffmpeg -f image2 -r 15 -i " + file_dir+file_name+ "/figs_%03d.png  -vcodec libx264 -y -vf scale=960:960 "+file_dir+file_name+'/vid2.mp4', shell=True, stdout=devnull, stderr=devnull)
    os.system('rm '+file_dir+file_name+'/figs*')

    #Combine videos into 960 x 480 vid
    devnull = open(os.devnull, 'w')
    subprocess.call("ffmpeg -i "+ file_dir+file_name+"/vid2.mp4 -i " + file_dir+file_name+"/vid1.mp4 -filter_complex '[0:v]pad=iw*2:ih[int];[int][1:v]overlay=W/2:0[vid]' -map [vid] -c:v libx264 -crf 23 -preset veryfast -y "
    + file_dir+file_name+'/'+file_name+'_unit'+str(unit)+'_'+plot_string+'_'+str(window)+'sec_window1.mp4', shell=True,stdout=devnull, stderr=devnull)

    #Crop video
    devnull = open(os.devnull, 'w')
    subprocess.call("ffmpeg -i " + file_dir+file_name+'/'+file_name+'_unit'+str(unit)+'_'+plot_string+'_'+str(window)+"sec_window1.mp4 -filter:v 'crop=1440:960:0:0' -c:a copy -y " 
    + file_dir+file_name+'/'+file_name+'_unit'+str(unit)+'_'+plot_string+'_'+str(window)+'sec_window.mp4', shell=True,stdout=devnull, stderr=devnull)

    #Delete old videos
    subprocess.call("rm " + file_dir + file_name+ "/vid1.mp4 "+ file_dir + file_name+ "/vid2.mp4 "+ file_dir+file_name+'/'+file_name+'_unit'+str(unit)+'_'+plot_string+'_'+str(window)+'sec_window1.mp4', shell=True, stdout=devnull, stderr=devnull)

    print "Finished unit: ", unit
    print ""
    print ""
    print ""
    #quit()


def Multitaper_specgram_allfreqs(time_series):

    s = time_series
    
    #print "Length of recording: ", len(s)
    
    #******************************************************
           
    #Multi taper function parameters
    nfft = 5000
    shift = 10
    Np = 1000
    k = 6
    tm = 6.0
    
    #Computing multi taper trf_specgram
    spec = tfr_spec(s, nfft, shift, Np, k, tm)

    #Martin changes
    zis = np.where(spec == 0.0) # row and column indices where P has zero power
    if len(zis[0]) > 0: # at least one hit
        spec[zis] = np.finfo(np.float64).max # temporarily replace zeros with max float
        minnzval = spec.min() # get minimum nonzero value
        spec[zis] = minnzval # replace with min nonzero values
    spec = 10. * np.log10(spec) # convert power to dB wrt 1 mV^2?

    p0=-40
    p1=None
    if p0 != None:
        spec[spec < p0] = p0
    if p1 != None:
        spec[spec > p1] = p1

    return spec, nfft
    
def Compute_lpf_static_maps(file_dir, file_name, images_raw, img_start, img_end, len_frame, img_rate, n_pixels, n_frames, img_times):

    #Load Generic Mask
    coords_generic=[]
    if (os.path.exists(file_dir + 'genericmask.txt')==True):
        print "Loading existing generic mask"
        generic_mask_file = file_dir + 'genericmask.txt'
        coords_generic = np.loadtxt(generic_mask_file)

    #Load Artifact Mask
    coords_artifact=[]
    if (os.path.exists(file_dir + 'artifactmask.txt')==True):
        print "Loading existing artifact mask"
        artifact_mask_file = file_dir + 'artifactmask.txt'
        coords_artifact = np.loadtxt(artifact_mask_file)

    #Load Bregma Mask
    coords_bregma=[]
    if (os.path.exists(file_dir + 'bregmamask.txt')==True):
        print "Loading existing bregma mask"
        bregma_mask_file = file_dir + 'bregmamask.txt'
        coords_bregma = np.loadtxt(bregma_mask_file)

    #PROCESS IMAGE TIMES FIRST
    #Test only on images in first X secs
    time_length = 120    #No. of secs of analysis
    temp0 = np.where(np.logical_and(img_times>=img_times[0], img_times<=img_times[0]+time_length))[0]
    img_times = img_times[temp0]

    print "No. images: ", len(img_times)
    print "First and last img time: ", img_times[0], img_times[-1]
    #print len(img_times)

    #selected_images = np.array(images_raw[temp0[:-1]],dtype=np.float32) #Remove last image as it is not being use
    selected_images = np.array(images_raw[temp0],dtype=np.float32) #Remove last image as it is not being use

    #remove baseline for [Ca]; but not VSD?
    if True:
        print "Removing baseline over all data"
        baseline = np.mean(selected_images, axis=0)
        selected_images = (selected_images - baseline)/baseline
    
    #PROCESS ephys data;
    sim_dir = file_dir
    sorted_file = file_name
    
    low_pass = False
    raw = True
    mua_load = False
    
    if low_pass:
        f = open(sim_dir + sorted_file+'/' + sorted_file + '_lp.tsf', "rb")
        tsf = Tsf_file(f, sim_dir)  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others
        tsf.sim_dir = sim_dir
        tsf.tsf_name = sim_dir + sorted_file+'/' + sorted_file + '.tsf'
        tsf.fname = sim_dir + sorted_file+'/' + sorted_file + '.tsf'
        f.close()
    
    if raw:
        f = open(sim_dir + sorted_file+'/' + sorted_file + '_raw.tsf', "rb")
        tsf = Tsf_file(f, sim_dir)  #Auto load tsf file attributes: n_electrodes, ec_traces, SampleFreqeuncy and others
        tsf.sim_dir = sim_dir
        tsf.tsf_name = sim_dir + sorted_file+'/' + sorted_file + '.tsf'
        tsf.fname = sim_dir + sorted_file+'/' + sorted_file + '.tsf'
        f.close()
    
    if mua_load:
        #Load Sorted data
        work_dir = sim_dir + file_name + "/"
        file_name = sorted_file + '_hp'
        ptcs_flag = 0
        Sort = Loadptcs(file_name, work_dir, ptcs_flag, save_timestamps=False)
        Sort.name=file_name
        Sort.filename=file_name
        Sort.directory=work_dir
       
        
    print "Sample freq: ", tsf.SampleFrequency
    print "Loaded .tsf for xcorrelation"
    
    #Load only ephys data within imaging window.
    #Find time indexes in ephys data that fall within imaging window - temp1

    if True:  #If using ephys data
        print "downsampling ephys to 1khz...", len(tsf.ec_traces[0])
        ephys_times = np.arange(0,float(img_times[-1])*1E3,1.E3/float(tsf.SampleFrequency))*1.E-3
        print "finding matching imaging frames for ephys times..."
        temp1 = np.where(np.logical_and(ephys_times>=img_times[0], ephys_times<=img_times[-1]))[0]
        print "Ephys time indexes: ", np.float32(temp1)/tsf.SampleFrequency

        #computing split_array used to chunk the ephys data: used to obtain average LFP during single img frame; timesteps of original data
        split_array = []
        for i in range(0,len(img_times)-1,1):
            split_array.append(int(len_frame*i*tsf.SampleFrequency))

    bands = [[0.1, 4.0], [4.0, 8.0], [8.0, 12.0], [12.0, 25.0], [25.0, 100.0], [500.0, 5000.0]]
    #bands = [[0.1,5000.]]
    bands_names = ['Delta (0.1-4.0)', 'Theta (4.0-8.0)', 'Alpha (8.0-12.0)', 'Beta (12.0-25.0)', 'Gamma (25.0-150.0)', 'High (500-5000)']
    #bands_names = ["all"]
    bands_out = ['delta', 'theta', 'alpha', 'beta', 'gamma', 'high']
    #bands_out = ['all']
    
    counter = 0
    top_row = True
    #electrodes = np.arange(0,16,1)
    electrodes = np.arange(0,16,1)
    #electrodes = [0,5,10,15]
    
    n_electrodes = len(electrodes)
    
    lfp_correlation_method = True
    lfp_power_method = False    #Requires split array +1 and selected_images - 1... not clear why.
    
    for q in electrodes:
        if lfp_correlation_method:
            for b in range(len(bands)):
                print "Channel: ", q, " band: ", bands_names[b]
                
                #Use LFP Correlation method: average lfp signal during each frame and mutiply by frame
                ephys_data = tsf.ec_traces[q][temp1].copy() * 1.0       #Can look at only postivie or neative data
                fs = 20000 
                lowcut, highcut = bands[b]
                ephys_data = butter_bandpass_filter(ephys_data, lowcut, highcut, fs, order = 2)

                ephys_data = np.clip(ephys_data, 0, 1E10) #Look only at positive power ;
                #ephys_data = np.clip(ephys_data, -1E10, 0) #Look only at negative power;
                
                #print splitting ephys_data into chunks triggered on imaging times: to get average value of LFP during single img frame
                ephys_split = np.split(ephys_data, split_array) #
                
                ephys_split_mean = []
                for i in range(len(ephys_split)):
                    ephys_split_mean.append(np.mean(ephys_split[i]))

                    
                ##Visualize img frame-to-ephys matching: LFP data overlayed with averaged chunks triggered on image times
                #if True:
                    #width = []
                    #for i in range(len(img_times)):
                        #width.append(10)
                    #plt.bar(img_times, width, .01, color='black', alpha=0.65)
                    #plt.plot(ephys_times[temp1], ephys_data, color='blue')
                    #print len(img_times), ephys_split_mean.shape
                    #plt.bar(img_times, ephys_split_mean, .03, color='pink', alpha=0.45)
                    #plt.show()
                
                #print "Computing mean of ephys_split data"
                ephys_split_mean = np.clip(ephys_split_mean,0,1000000) #remove negative values;
                ephys_split_mean = np.float32(ephys_split_mean)
                ephys_split_mean_max = np.nanmax(ephys_split_mean)     #Normalize to largest LFP value (negative or positive)
                ephys_split_mean = ephys_split_mean/ephys_split_mean_max    
                #ephys_split_mean = np.clip(ephys_split_mean,0,1) #remove negative values;
                ephys_split_mean = np.nan_to_num(ephys_split_mean)

                #Compute lfp average triggered images
                image_lfp = np.einsum('m,mdr->mdr',ephys_split_mean, selected_images)
                image_lfp = np.mean(image_lfp, axis=0)

                np.save(file_dir + file_name+'/'+file_name+'_lfpmap_band_'+bands_out[b]+'_channel_'+str(q).zfill(2), image_lfp)

        if lfp_power_method:
            ephys_data = tsf.ec_traces[q][temp1].copy()     #This selects only part of recording
            fs = 20000 
            

            #ephys_temp = butter_bandpass_filter(ephys_data, lowcut, highcut, fs, order = 2)
            print "Computing time-frequency reassignment specgram channel: ", q
            
            tfr_file = file_dir + file_name+'/'+file_name+'_tfr_channel_'+str(q).zfill(2)
            if os.path.exists(tfr_file+'.npz'): 
                data = np.load(tfr_file+'.npz')
                mt_specgram = data['mt_specgram']
                nfft = data['nfft']
            else: 
                mt_specgram, nfft = Multitaper_specgram_allfreqs(ephys_data)
                np.savez(tfr_file, mt_specgram=mt_specgram, nfft=nfft)
            
            for b in range(len(bands)):
                f0, f1 = bands[b]
                lo = int(float(nfft)/1000. * f0)
                hi = int(float(nfft)/1000. * f1)
                mt_specgram_temp = mt_specgram[lo:hi][::-1] #Take frequency band slice only; also invert the data

                #print f0, f1
                #plt.imshow(mt_specgram_temp, extent = [0,f1, len(mt_specgram_temp),0], origin='upper', aspect='auto', interpolation='sinc')
                #plt.show()
            
                spec_split = np.array(split_array)*len(mt_specgram_temp[0])/split_array[-1] #Normalize splitting array to length of spectrogram
                #print spec_split
                spec_split_mean = []
                for ss in range(len(spec_split)-1):
                    spec_split_mean.append(np.mean(mt_specgram_temp[:,spec_split[ss]:spec_split[ss+1]]))

                #print "Computing mean of ephys_split data"
                spec_split_mean = np.array(spec_split_mean)
                spec_split_mean_max = np.max(np.abs(spec_split_mean))     #Normalize to largest LFP value (negative or positive); THIS CAN"T BE NEGATIVE!
                spec_split_mean_min = np.min(np.abs(spec_split_mean))     #Normalize to largest LFP value (negative or positive); THIS CAN"T BE NEGATIVE!
                spec_split_mean = (spec_split_mean-spec_split_mean_min)/(spec_split_mean_max-spec_split_mean_min)
                #plt.plot(spec_split_mean)
                #plt.show()
                #quit()
                #spec_split_mean = np.clip(spec_split_mean,0,1) #remove negative values; NOT REQUIRED

                #Compute lfp average triggered images
                print spec_split_mean.shape, selected_images.shape
                image_lfp = np.einsum('m,mdr->mdr', spec_split_mean, selected_images)
                image_lfp = np.mean(image_lfp, axis=0)

                np.save(file_dir + file_name+'/'+file_name+'_powermap_band_'+bands_out[b]+'_channel_'+str(q).zfill(2), image_lfp)

            #print "Plotting image_lfp"
            ax = plt.subplot(n_electrodes, 7, counter+1)
            ax.get_xaxis().set_visible(False)
            ax.set_yticklabels([])
            if top_row:
                plt.title(bands_names[counter])

            #Ablate generic map
            min_pixel = np.min(image_lfp)
            for j in range(len(coords_generic)):
                image_lfp[min(n_pixels,int(coords_generic[j][0]))][min(n_pixels,int(coords_generic[j][1]))]=min_pixel
            for j in range(len(coords_artifact)):
                image_lfp[min(n_pixels,int(coords_artifact[j][0]))][min(n_pixels,int(coords_artifact[j][1]))]=min_pixel
            for j in range(len(coords_bregma)):
                for k in range(n_pixels):
                    for l in range(7):
                        image_lfp[k][min(n_pixels-1,int(coords_bregma[1])-3+l)]=min_pixel
                        
            plt.imshow(image_lfp, origin='lower')
            if counter%7==0:
                plt.ylabel("Ch: "+str(q))
            plt.xlim(0,n_pixels-1)
            plt.ylim(n_pixels-1, 0)
            counter+=1


        counter+=1
        top_row=False
    


def animate_data(data):
    '''Make movies from data matrix
    '''
    v_max=np.nanmax(data)
    v_min=np.nanmin(data)
    print v_max, v_min
    
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)

    fig = plt.figure() # make figure

    # make axesimage object
    im = plt.imshow(data[0], cmap=plt.get_cmap('jet'), vmin=v_min, vmax=v_max, interpolation='none')#, vmin=0, vmax=v_max)
    # function to update figure
    def updatefig(j):
        # set the data in the axesimage object
        im.set_array(data[j])
        plt.title("Frame: "+str(j)+"\n"+str(round(float(j)/150,2))+"sec")
        # return the artists set
        return im,
    # kick off the animation
    ani = animation.FuncAnimation(fig, updatefig, frames=range(len(data)), interval=100, blit=False, repeat=True)

    if False:
        ani.save(fname+'.mp4', writer=writer)
    
    plt.show()

def mask_data_single_frame(data, main_dir, midline_mask):
    
    n_pixels = len(data)
            
    #Load General mask (removes background)
    generic_mask_file = []
    generic_mask_file = main_dir + 'genericmask.txt'        #Mask is in 256 x 256 resolution
    if (os.path.exists(generic_mask_file)==False):
        generic_coords = Define_generic_mask(data, main_dir)
    else:
        generic_coords = np.int32(np.loadtxt(generic_mask_file))
        
    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)):
        generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True

    if n_pixels != 256:                                     #Subsample mask if required
        generic_mask_indexes = scipy.misc.imresize(generic_mask_indexes,n_pixels/256.)

    #Load midline mask
    for i in range(midline_mask):
        generic_mask_indexes[:,n_pixels/2+int(midline_mask/2)-i]=True
        
    #temp_array = np.ma.array(np.zeros((len(data),n_pixels,n_pixels),dtype=np.float32), mask=True)
    #Mask all frames; NB: PROBABLY FASTER METHOD
    #for i in range(0, len(data),1):
    temp_array = np.ma.masked_array(data, mask=generic_mask_indexes)
    
    return temp_array
    
def mask_data(data, main_dir, midline_mask):
    
    n_pixels = len(data[0])
            
    #Load General mask (removes background)
    generic_mask_file = []
    generic_mask_file = main_dir + 'genericmask.txt'        #Mask is in 256 x 256 resolution
    if (os.path.exists(generic_mask_file)==False):
        generic_coords = Define_generic_mask(data, main_dir)
    else:
        generic_coords = np.int32(np.loadtxt(generic_mask_file))
        
    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)):
        generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True

    if n_pixels != 256:                                     #Subsample mask if required
        generic_mask_indexes = scipy.misc.imresize(generic_mask_indexes,n_pixels/256.)

    #Load midline mask
    for i in range(midline_mask):
        generic_mask_indexes[:,n_pixels/2+int(midline_mask/2)-i]=True
        
    temp_array = np.ma.array(np.zeros((len(data),n_pixels,n_pixels),dtype=np.float32), mask=True)
    #Mask all frames; NB: PROBABLY FASTER METHOD
    for i in range(0, len(data),1):
        temp_array[i] = np.ma.masked_array(data[i], mask=generic_mask_indexes)
    
    return temp_array
    


def PCA(X, n_components):
    from sklearn import decomposition

    pca = decomposition.PCA(n_components)
    pca.fit(X)
    X=pca.transform(X)

    coords = []
    for i in range(len(X)):
         coords.append([X[i][0], X[i][1], X[i][2]])
    
    return X, np.array(coords).T














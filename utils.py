import subprocess

import numpy as np
import math
import os.path
import multiprocessing as mp
from matplotlib import animation
from matplotlib.path import Path
import glob 
import skimage
from skimage import data
from skimage.transform import rotate
            
from pylab import *
import csv

import matplotlib.patches as mpatches

from scipy import ndimage

import PIL
from PIL import ImageFont
from PIL import Image
from PIL import ImageDraw


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


def Load_images_start_end(file_dir, file_name, images_raw):
    ''' This method aligns imaging record to electrophysiology by identifying the on and off 
        trigger times from the original Multi-Channel-Systems .mcd 17th channel. The data has already been parsed and saved as a txt file: ephys_times.txt
        
        Other electrophysiology files will require different code for identifying the trigger time.
    '''
    
    ephys_times_file = file_dir + file_name+ '/ephys_times.txt'
    temp_data= np.loadtxt(ephys_times_file)

    start_array = temp_data[0]
    end_array = temp_data[-1]

    #Find star/end for multi-part recordings
    img_start = start_array
    img_end = end_array
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
    
    return len_frame, img_rate, n_pixels, img_times 
    
    
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
    
def Compute_sta_motif(unit, channel, all_spikes, window, img_rate, img_times, n_pixels, images_aligned, file_dir, file_name, n_procs, overwrite, stm_types, random_flag, spiking_modes):
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
            images_temp = np.array(images_aligned.copy(), dtype=np.float32) #This is an inneficient and memory consuming approach; should modify eventually

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
                    print file_out[:-4]+"_imagingspikes_grouped_spikes.txt"
                    
                    #Here loading the spikes for each mode; Each row contains the spikes for each mode in txt format to be human readable
                    if os.path.exists(file_out[:-4]+"_imagingspikes_grouped_spikes.txt")==False:
                        print "... spiking mode file '*_grouped_spikes.txt' missing..."
                        print "... can only generate 'all' spiking modes..."
                        print "See sta_maps.py code for details on how to generate mode spiking file..."
                        quit()
                        
                    mode_spikes = []
                    mode_names = []
                    with open(file_out[:-4]+"_imagingspikes_grouped_spikes.txt", 'rt') as inputfile:
                        reader = csv.reader(inputfile)
                        for row in reader:
                            if len(row)<20: 
                                mode_names.append(row[0])               #Read names of modes; They should all be less than 20 characters
                            else:
                                mode_spikes.append(np.float32(row))     #Break the loop, first save the data as spiking data.
                                break
                                
                        for row in reader:
                            mode_spikes.append(np.float32(row))
                    
                    #Load spiking modes from disk
                    for n in range(len(spiking_modes)): 
                        for m in range(len(mode_spikes)): 
                            if mode_names[m]!=spiking_modes[n]: continue       #Search spiking rasters text file for matching mode indicated in spiking_modes file
    
                            npy_file_name = file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+ \
                                            stm_type+'_'+str(window)+'sec_window_'+str(len(spikes)).zfill(5)+"_spikes_"+spiking_modes[n]
                            
                            if os.path.exists(npy_file_name+'.npy')==False: 
                                            
                                temp3 = []
                                for i in mode_spikes[m]:        #temp3 contains all the frame indexes from img_times for each spike in raster; e.g. 180 frames for each spike automatically aligned
                                    temp3.append(np.where(np.logical_and(img_times>=i-window, img_times<=i+window))[0][0:2*int(window*img_rate)]) #Fixed this value or could be off +/-1 frame
                        
                                #Initialize list to capture data and store frame sequences;  Make MP Pool for spiking mode
                                print "... computing spiking mode: ", spiking_modes[n], "  #spikes: ", len(mode_spikes[m])
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
   
def View_sta_motif(unit, main_dir, file_dir, file_name, stm_types, img_rate, spiking_modes, n_spikes):
    
    print "... viewing STA motifs..."
    
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

        #Plot color bar
        cbar = fig.colorbar(im, ticks = [v_min, 0, v_max], ax=ax, fraction=0.02, pad=0.05, aspect=3)
        cbar.ax.set_yticklabels([str(round(v_min*100,1))+"%", '0'+"%", str(round(v_max*100,1))+"%"])  # vertically oriented colorbar
        cbar.ax.tick_params(labelsize=15) 
        
    plt.suptitle("Experiment: " + file_name + "  unit: "+str(unit) + "   #spikes: "+str(n_spikes), fontsize=15)
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



def Compute_STM(img_rate, window, n_procs, main_dir, file_dir, file_name, n_pixels, unit, channel, n_spikes, ptp, stm_types, spiking_modes):

    ''' Compute STM using pre- and post- spike intervals to obtain static map; See also the "max" version of this function'''
    
    print "Computing STMs"
    
    plot_strings=[]
    stm_filenames=[]

    for mode_ in stm_types:
        if mode_=="modes":
            plot_strings.extend(spiking_modes)    #Need these new modes for the revision analysis
            
            for spiking_mode in spiking_modes:
                temp_name = glob.glob(file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_modes_*'+spiking_mode+'.npy')
    
                if len(temp_name)!=1:
                    print "...missing motif file (or too many files). Must run 'compute_stm_motif' flag to generate sta motif first..."
                    return
                else: 
                    stm_filenames.append(temp_name[0])
        else:
            plot_strings.extend(mode_)
            temp_name = glob.glob(file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_all_*.npy')
            if len(temp_name)!=1:
                print "...missing motif file (or too many files). Must run 'compute_stm_motif' flag to generate sta motif first..."
                return
            else: 
                stm_filenames.append(temp_name[0])
        
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

        #Compute Maxmaps    (can also compute amin maps)
        temp_array = []
        for i in range(-int(img_rate), int(img_rate),1):    #Load -1sec ... + 1sec data and look for max value for each pixel
            temp = np.ma.array(images_areas[int(window*img_rate)+i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
            temp_array.append(temp)

        temp_array = np.float32(temp_array)
        images_out2 = np.amax(temp_array,axis=0) #Search for max-pixel
        images_out2 = np.ma.array(images_out2, mask=generic_mask_indexes, fill_value = 0)
        images_out2 = np.ma.filled(images_out2, 0.0)
        images_out2 = np.nan_to_num(images_out2)

        #Display max pixel maps
        midline_mask = 3        #Number of pixels to mask at the midline

        fig, ax = plt.subplots()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        image_loaded = images_out2 #np.load(stm_filename[:-4]+"_maxmaps.npy")
        image_loaded = mask_data_single_frame(image_loaded, main_dir, midline_mask)

        #Process image loaded
        v_max = np.nanmax(np.abs(image_loaded)); v_min = -v_max

        im = plt.imshow(image_loaded, vmin=v_min, vmax=v_max)

        plt.title("Experiment: " + file_name + "   unit:    " + str(unit) + "   #spikes: "+str(n_spikes), fontsize=20)

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
    ''' Paralellized spaghetti code to manually mask out 
    '''
    
    global coords_temp
    
    
    for j in range(len(coords_temp)):
        if temp_array1[coords_temp[j][0]][coords_temp[j][1]] == 0:  #Only overwrite non-masked values
            temp_array1[coords_temp[j][0]][coords_temp[j][1]] = -1000

    return temp_array1

    
def Mp_min((temp_array2)):
    global coords_temp

    for j in range(len(coords_temp)):
        if temp_array1[coords_temp[j][0]][coords_temp[j][1]] == 0:  #Only overwrite non-masked values
            temp_array2[coords_temp[j][0]][coords_temp[j][1]] = 1000

    return temp_array2


def Compute_STMTD(unit, channel, spikes, file_dir, file_name, img_rate, window, n_procs, area_names, depth, sides, stm_types, spiking_modes):

    ''' This routine loads specific ROIs of interest, identifies the highest and lowest active pixel in each ROI during the motif time (e.g. 6sec)
        and then saves the activity at that pixel for the duration of the motif.
        - each ROI STMTD is computed for both sides of the brain and for both the max_value pixel and min_value pixel.
        - data is saved into a readable .txt file with # spikes at the top, STMTD names, STMTD traces, and the location in pixels of each of the STMTDs.
    
        STMTDs are computed in regions of interest (ROIs) that are defined for each cortical and subcortical recording and saved in the /roi/ directory.
        - each ROI is defined as 256 x 256 single frame with the ROI coordinates set to 0 and the non-ROI coordinates set to 1.
        - each ROI has a contour defined which just outlines the edges of the ROI and is saved as an .npy file with several coordinates in it. These contours
          are used in the animations to indicate where each ROI was approximately computed. The code should run without these 
    
    '''

    global coords_temp

    #Convert firing modes to a list of strings
    plot_strings=[]
    stm_filenames=[]
    for mode_ in stm_types:
        if mode_=="modes":
            plot_strings.append(spiking_modes)    #Need these new modes for the revision analysis
            for spiking_mode in spiking_modes:
                stm_filenames.append(glob.glob(file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_modes_*'+spiking_mode+'.npy')[0])
        else:
            plot_strings.append(mode_)
            stm_filenames.append(glob.glob(file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_all_*.npy')[0])
                    
    
    for plot_string, stm_filename in zip(plot_strings, stm_filenames):
    
        images_processed = np.load(stm_filename)
        n_pixels = images_processed.shape[1]

        out_file = file_dir + file_name + '/time_course_data_' + file_name+'_'+plot_string+ '_'+str(window)+'sec_window_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'.txt'
        
        if os.path.exists(out_file): continue

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
                
                #Not sure all this masked array architecture is required.
                temp_max_array = []
                temp_min_array = []
                temp_array1 = np.ma.array(np.zeros((int(img_rate*2),n_pixels,n_pixels)), mask=True)
                temp_array2 = np.ma.array(np.zeros((int(img_rate*2),n_pixels,n_pixels)), mask=True)
                for i in range(int(img_rate*2)):     #This searches +/- 1 sec from time = 0 sec; OR Maximum to begining of data;
                    temp_array1[i] = np.ma.array(temp_array[int(window*img_rate)-int(img_rate)+i]).copy()
                    temp_array2[i] = np.ma.array(temp_array[int(window*img_rate)-int(img_rate)+i]).copy()


                #temp_array1._sharedmask=False
                #temp_array2._sharedmask=False

                temp_max_array=temp_array1
                temp_min_array=temp_array2
                
                #Search for max value using 1D unravel of temp_arrays; assign location of max/min index; detect max/min values overall;
                temp1_array = np.ma.array(np.zeros((int(img_rate*2),n_pixels,n_pixels)), mask=True)
                for k in range(len(temp_max_array)):
                    temp1_array[k] = np.ma.array(temp_max_array[k], mask=area_mask) #np.ma.array(temp_max_array)

                temp_max_array = temp_array1
                max_index = np.unravel_index(np.nanargmax(temp_max_array), temp_max_array.shape)
                max_pixel_value = temp_array[int(window*img_rate)-int(img_rate)+max_index[0]][max_index[1]][max_index[2]]

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
                
                counter+=1  #Keep track of each STMTD in a different line to be saved below


        #Make save array and save number of spikes
        temp_array = []
        temp_array.append([len(spikes)])  #Save # spikes within imaging period
        
        #Save each time course label:
        for area in area_names:
            for side in sides: 
                for maxmin in ['max_value', 'min_value']:
                    temp_array.append([area+'_'+side+"_"+maxmin])           
                    
        #Add all time courses to the file in same order as above.
        for i in range(len(Max_plot)):
            max_plot=np.array(Max_plot[i])
            temp_array.append(max_plot)

            min_plot=np.array(Min_plot[i])
            temp_array.append(min_plot)
            
            
        #Add max_index and min_index information to the end of the time_course_data*.txt file
        for i in range(len(Max_index)):
            temp_array.append(Max_index[i])
            temp_array.append(Min_index[i])


        #Save time_course_data
        with open(out_file, "w") as f:
            writer = csv.writer(f)
            writer.writerows(temp_array)
        
        np.save(out_file[:-4]+"_images_areas", Images_areas)

    #return Max_plot, Min_plot, Max_pixel_value, Min_pixel_value, Max_index, Min_index, area_names, Images_areas

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
    


def View_STMTD(unit, channel, spikes, window, len_frame, file_dir, file_name, area_names, sides, stm_types, spiking_modes):
               
    ''' Loads and displays STMTDs already computed    
    '''
    
    print "...viewing STMTDs..."
    
    colours = ['blue', 'red','green', 'magenta', 'cyan', 'orange', 'brown', 'pink', 'yellow', 'grey']
    labelsize = 30

    #Convert firing modes to a list of strings
    plot_strings=[]
    for mode_ in stm_types:
        if mode_=="modes":
            plot_strings.extend(spiking_modes)      #spiking modes are 'first', 'last', 'burst', 'tonic'
        else:
            plot_strings.extend([mode_])            #other spiking modes are 'all'; Others can be added
    

    for plot_string in plot_strings:
        in_file = file_dir + file_name + '/time_course_data_' + file_name+'_'+plot_string+ '_'+str(window)+'sec_window_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'.txt'
        #Save time_course_data
        with open(in_file, "r") as f:
            rows = list(csv.reader(f))
            n_spikes = rows[0]
            del rows[0]
            
            #Load time course labels
            roi_names=[]; ctr=0
            for area in area_names:
                for side in sides: 
                    for maxmin in ['max_value', 'min_value']:
                        roi_names.append(rows[ctr])
                        ctr+=1
            
            #Load time course labels
            roi_stmtds=[]
            for area in area_names:
                for side in sides: 
                    for maxmin in ['max_value', 'min_value']:
                        roi_stmtds.append(rows[ctr])
                        ctr+=1
            
        
        #Plot left side STMTDs
        ax = plt.subplot(1,2,1)
        for i in range(len(area_names)):
            print "...plotting left side max value STMTDs..."
            
            xx = np.linspace(-window, window, len(roi_stmtds[0]))
            plt.plot(xx, np.float32(roi_stmtds[i*4])*100., color=colours[i], linewidth=3, alpha=0.8)
            

        ax.plot([-3,3],[0,0], color='black')
        ax.plot([0,0],[-5.0,5.0], color='black')
        ax.set_ylim(-5, 5)
        ax.tick_params(axis='both', which='major', labelsize=labelsize)
        ax.yaxis.set_ticks(np.arange(-5,5,1))
        
        plt.title("Left", fontsize = labelsize)
        ax.get_yaxis().set_visible(True)
        #plt.ylabel(area_names[int(i/2)], fontsize=labelsize)
        plt.ylabel("DF/F %", fontsize = labelsize)
        plt.xlabel("Time (sec)", fontsize = labelsize)

        #Plot Legend
        legend_array = []
        labels = []
        for ctr, area_name in enumerate(area_names): 
            legend_array.append(mpatches.Patch(facecolor = colours[ctr], edgecolor="black"))
            labels.append(area_name)

        ax.legend(legend_array, labels, fontsize=12, loc=0)

        #Plot right side STMTDs
        ax = plt.subplot(1,2,2)
        for i in range(len(area_names)):
            print "...plotting right side max value STMTDs..."
            
            xx = np.linspace(-window, window, len(roi_stmtds[0]))
            plt.plot(xx, np.float32(roi_stmtds[i*4+2])*100., color=colours[i], linewidth=3, alpha=0.8)
            

        ax.plot([-3,3],[0,0], color='black')
        ax.plot([0,0],[-5.0,5.0], color='black')
        ax.set_ylim(-5, 5)
        ax.tick_params(axis='both', which='major', labelsize=labelsize)
        ax.yaxis.set_ticks(np.arange(-5,5,1))
        
        plt.title("Right", fontsize = labelsize)
        ax.get_yaxis().set_visible(True)
        plt.ylabel("DF/F %", fontsize = labelsize)
        plt.xlabel("Time (sec)", fontsize = labelsize)


        plt.suptitle(file_name + " Unit: " + str(unit).zfill(2) + " Channel: " + str(channel).zfill(2) + " No. of spikes: " + str(n_spikes), fontsize=labelsize)

        #npy_file_name = file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+plot_string+'_'+str(window)+'sec_window'

        plt.show()
        
        #plt.savefig(file_dir + file_name + '/time_course_plot_' + file_name+'_'+plot_string+ '_'+str(window)+'sec_window_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'.png', fontsize = 20)
        
        #plt.close()



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
    


def Animate_images(unit, channel, window, img_rate, main_dir, file_dir, file_name, n_pixels, spikes, plot_string, n_procs, area_names, sides, depth, stm_types, spiking_modes):
    
    #(unit, channel, window, img_rate, Images_areas, file_dir, file_name, n_pixels, spikes, plot_string, n_procs, 
    #Max_plot

    colors = ['blue', 'red','green', 'magenta', 'cyan', 'orange', 'brown', 'pink', 'yellow', 'grey']

    #print "Generating Animation. No. of Frames: ", len(Images_areas[0])
    
    #Load General mask (removes background)
    generic_mask_file = []
    generic_mask_file = main_dir + 'genericmask.txt'
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.int16(np.loadtxt(generic_mask_file))
    
    print "... n_pixels: ", n_pixels
    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)):
        generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True
        

    #Convert firing modes to a list of strings
    plot_strings=[]
    for mode_ in stm_types:
        if mode_=="modes":
            plot_strings.extend(spiking_modes)      #spiking modes are 'first', 'last', 'burst', 'tonic'
        else:
            plot_strings.extend([mode_])            #other spiking modes are 'all'; Others can be added
        
        
    #Load Max_plot and other values.
    for plot_string in plot_strings:
        in_file = file_dir + file_name + '/time_course_data_' + file_name+'_'+plot_string+ '_'+str(window)+'sec_window_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'.txt'
        #Save time_course_data
        with open(in_file, "r") as f:
            rows = list(csv.reader(f))
            n_spikes = rows[0]
            del rows[0]

            #Load time course labels
            roi_names=[]; ctr=0
            for area in area_names:
                for side in sides: 
                    for maxmin in ['max_value', 'min_value']:
                        roi_names.append(rows[ctr])
                        ctr+=1

            #Load time course labels
            Max_plot=[]
            for area in area_names:
                for side in sides: 
                    for maxmin in ['max_value', 'min_value']:
                        Max_plot.append(rows[ctr])
                        ctr+=1

    Max_plot_left = np.float32(Max_plot)[::4] *100.  #Convert data to percent
    Max_plot_right = np.float32(Max_plot)[2::4] *100.  #Convert data to percent

    #Process plots
    plot_frames = []
    x = np.arange(-window, window, 1./img_rate)
    x = x[:len(Max_plot_left[0])]
    #Max_plot = np.array(Max_plot)

    print "Number of ROIs: ", len(Max_plot_left)
    ranges = np.arange(0,len(x),1)


    n_chunks = 8
    range_chunks = np.array_split(ranges,n_chunks)
    
    plt.close()
    
    #PLOT CURVE GRAPHS
    if False: #PLOT ALL CURVES FOR ALL AREAS
        print "... processing plot frames..."
        for time_step in range(len(Max_plot_left[0])):
            fig = plt.figure()
            fig.add_subplot(211)
            fig.set_size_inches(5, 20)
            
            plt.xlabel("Time (s)")

            plt.plot([-window,window],[0,0], color='black')
            plt.plot([0,0],[-5,5], color='black')
            plt.title("Left Hemisphere", fontsize=20)

            plt.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='off') # 

            #DRAW CURVES OVER TIME
            for i in range(len(Max_plot_left)):
                plt.plot(x[:time_step], Max_plot_left[i][:time_step], color=colors[i], linewidth=2) #Plot individual overlayed left/rigth curves
            
            plt.ylim(-5,5)

            fig.add_subplot(212)
            fig.set_size_inches(5, 20)

            plt.title("Right Hemisphere", fontsize=20)
            plt.xlabel("Time (s)")
            plt.plot([-window,window],[0,0], color='black')
            plt.plot([0,0],[-5,5], color='black')

            plt.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='off') # 

            #DRAW CURVES OVER TIME
            for i in range(len(Max_plot_right)):
                plt.plot(x[:time_step], Max_plot_right[i][:time_step], color=colors[i], linewidth=2) #Plot individual overlayed left/rigth curves

            plt.ylim(-5,5)
                
            #Save each figure to disk; make movie below
            fig.savefig(file_dir + file_name+'/figs_'+str(time_step).zfill(3)+'.jpg', dpi=40)

            plt.close(fig)

        #Make vid1 containing time course curves
        devnull = open(os.devnull, 'w')
        subprocess.call("ffmpeg -f image2 -r 15 -i " + file_dir+file_name+ "/figs_%03d.jpg  -vcodec libx264 -y -vf scale=480:960 "+file_dir+file_name+'/vid1.mp4', shell=True,stdout=devnull, stderr=devnull)
        os.system('rm '+file_dir+file_name+'/figs*')


    #***************** Process images ****************
    
    Images_areas = np.load(in_file[:-4]+"_images_areas.npy")

    Images_areas = Images_areas*100.    #Convert to % scale
    Images_areas = np.swapaxes(Images_areas, 0, 1)
    
    
    #Generate borders around recorded areas and mask them out in order for color to work properly
    img=Image.new("RGBA", (n_pixels,n_pixels),(0,0,0))
    if True:
        draw = ImageDraw.Draw(img)
        counter=0
        for area in area_names:
            for side in sides:
                area_file = file_dir+ 'roi/'+depth+'_' + area+'_'+side+"_contour"
                print area_file
                if (os.path.exists(area_file+'.npy')==True):
                    area_contour = np.load(area_file+'.npy')     #Load coordinates; 2D vectors stacked; Mask out everything else

                for i in range(len(area_contour)-1):
                    if area == 'retrosplenial':  print area_contour[i][1],area_contour[i][0],area_contour[i+1][1],area_contour[i+1][0]
                    draw.line((area_contour[i][1],area_contour[i][0],area_contour[i+1][1],area_contour[i+1][0]), fill = colors[counter], width = 1)
                    #draw.line((area_contour[i][1],area_contour[i][0],area_contour[i+1][1],area_contour[i+1][0]), fill = 'white', width = 1)
            
            counter+=1
            
        borders_=asarray(img)

        if True:
            #Generate borders mask to zero out [Ca] image pixels to get color right
            borders_mask = np.zeros((n_pixels,n_pixels),dtype=bool)
            for i in range(n_pixels):
                for j in range(n_pixels):
                    if np.sum(borders_[i][j])>255:
                        borders_mask[i][j]=True

                
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
    
    #REDO THESE*******************************************************************
    v_min = np.amin(Max_plot_left)
    v_max = np.amax(Max_plot_left)
    print v_min, v_max

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














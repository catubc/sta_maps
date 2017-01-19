import numpy as np
import time, math
import os.path
import multiprocessing as mp
import re
import sys
import struct, array, csv
import scipy.optimize 
import matplotlib.mlab as mlab

from scipy import stats
from scipy import signal
from pylab import *
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from matplotlib import animation
from PIL import Image

from sta_utils import *

#***********************************************
np.seterr(divide='ignore', invalid='ignore')        #Ignore division by zero due to divide frames that were rotated and contain unmasked pixels

#******************************** SET LOADING DIRECTORIES DEFAULTS *************
main_dir = ''

file_dirs = []
file_names = []

#############**************** LOAD ALL FILES DATA IN EXPERIMENT DIRECTORY ****************
file_dirs.append('2015-7-22/') 
file_names.append([
'2015-7-22-1',     #cortical recording          #Number of cells: 9
#'2015-7-22-4',     #subcortical recording      #Number of cells: 30
])

#Select 
cell_list = [1]

#******Processing info******
window = 3      #Window of STM motif/montages
n_procs=14      #Number of processors to use for parallel sections of code

#Meta data for generating videos and STMTDS post STM processing
area_names = ['hindlimb', 'forelimb', 'barrel','retrosplenial','visual', 'motor', 'pta', 'acc'] 
sides = ['left','right']

#Spiking modes to use for computing STMs
spiking_modes = ['burst', 'last', 'tonic', 'first']
stm_types = ["all"]     #Types of firing modes ["all", "modes"] #Other types not in use;

    #NB: for 'modes' data need to generate an input text file with extention "_grouped_spikes.txt" containing spike times (in seconds) from each mode on separate lines; 
    #For example: filename:  "unit_01_channel_06_ptp_037_imagingspikes_grouped_spikes.txt"
    #line 1: "28.072, 29.431, 30.626, 30.640 ..."
    #line 2: "14.944, 15.889, 28.089, 29.490 ..."
    #line 3: "8.657, 9.998, 12.134, 12.835 ..."
    #line 4: "13.547, 14.98, 26.120, 28.308 ..."

    #Additional modes can be defined by changing the spiking modes name in the code; Or by simply re-purposing existing labels to a different class of spikes

#Flast for computing sta maps;
sta_maps = True
overwrite = False                        #Overwrite existing STM data if already computed
view_stms = False            #View each unit's STM after processing
search_max_min = True

#Set this flag along with the correct unit number to compute STMs for all motifs in a recording; This is used for other processing
random_flag = False  

#Static maps - SUA and LFP
sua_static_maps = False             #make maxmap and minmap data files 
lfp_static_maps = False

#****************************************************************************************
#Loop over experiments 
for dir_counter, file_dir in enumerate(file_dirs):

    for file_name in file_names[dir_counter]:

        #Load units from .csv file name; 
        files = os.listdir(file_dir+file_name)
        temp_names = []
        for file_ in files:
            if (file_[:5]=="unit_") and (".csv" in file_): temp_names.append(file_)
        
        #Save individual unit names, channels, ptps
        units = []
        channels = []
        ptps = []
        for file_ in temp_names:
            units.append(int(file_[5:7]))
            channels.append(int(file_[16:18]))
            ptps.append(int(file_[23:26]))

        #Load depth of recording and anesthetic state:
        infile = file_dir+file_name+"/state.txt"
        with open(infile, "r") as f:
            for line in f:
                state = line.rstrip()
        infile = file_dir+file_name+"/depth.txt"
        with open(infile, "r") as f:
            for line in f:
                depth = line.rstrip()
            
        #Load raw images and rotate + align them
        file_name_aligned = file_dir + file_name+'/'+file_name+'_images_aligned.npy'
        images_aligned = np.load(file_name_aligned)

        #Load start and end times of imaging system and compute required parameters
        img_start, img_end, len_frame, img_rate, n_pixels, n_frames, img_times = Load_images_start_end(file_dir, file_name, images_aligned)
        print "Imaging rate: ", img_rate
            

        #********************************COMPUTE SPIKE-TRIGGERED-MOTIFS********************************************
        #Compute and plot motifs: multi-frame representations of activity locked to spiking
        if sta_maps:
            print "No. units: ", len(units)
            print "OVERWRITE: ", overwrite

            #**************** LOOP OVER ALL UNITS IN EXPERIMENT **************************
            for i in range(len(units)):
                print "*Processing ", file_name, " unit: ", i+1, " of ", len(units), " unit name: ", units[i]
                unit=units[i]
                channel=channels[i]
                ptp=ptps[i]
                
                #Process specific units only from cell_list defined above
                if True:
                    if unit not in cell_list: continue            
    
                #Load sorted spikes for unit
                spikes = np.loadtxt(file_dir+file_name+'/unit_'+str(unit).zfill(2)+ '_channel_' + str(channel).zfill(2) + '_ptp_'+str(ptp).zfill(3)+'.csv')
                
                #Determine if there are spikes in recording window to continue; Required as some units are sparsely firing
                temp0 = np.where(np.logical_and(spikes>=img_times[0]+window, spikes<=img_times[-1]-window))[0]
                spikes_in_window = spikes[temp0]
                if len(spikes_in_window) == 0: 
                    print "Zero spikes in window - skip unit"; continue
                
                n_spikes = len(spikes_in_window)  #NB: Save n_spikes value just for current epoch for continuous recordings

                #Compute and view motifs
                images_processed, spikes, plot_string = Compute_spike_triggered_average(unit, channel, spikes, window, img_rate, img_times, n_pixels, 
                                                            images_aligned, file_dir, file_name, n_procs, overwrite, stm_types, random_flag, spiking_modes)

                #Preview STM for each cell
                if view_stms:
                    view_static_stm(unit, main_dir, file_dir, file_name, stm_types, img_rate, spiking_modes)


                #Search Max/Min, Save Time Courses
                if search_max_min:
                    
                    #Search for Min and Max pixels
                    Max_plot, Min_plot, Max_pixel_value, Min_pixel_value, Max_index, Min_index, area_names, images_areas = \
                        Search_max_min(unit, channel, spikes, file_dir, file_name, img_rate, window, n_procs, area_names, depth, sides, stm_types, spiking_modes)

                    #Save time course .npy file and figures
                    print "Saving time course for unit : ", unit, " ..."
                    Save_time_course(unit, channel, spikes, Max_plot, Min_plot, Max_index, Min_index, window, len_frame, file_dir, file_name, area_names, sides, stm_types, spiking_modes)


                if animate_images: 
                    Animate_images(unit, channel, window, img_rate, images_areas, file_dir, file_name, n_pixels, spikes, plot_string, n_procs, generic_mask_indexes, 
                        Max_plot, Min_plot, Max_pixel_value, Min_pixel_value, Max_index, Min_index, area_names, sides, depth)
                
                
                #*********** Generate ROI matrix and time courses
                if False:
                    
                    images_areas = Load_areas_and_mask(depth, unit, channel, n_pixels, main_dir, file_dir, file_name, images_processed, area_names, sides)
                    
                    average_areas = Average_roi(images_areas, img_rate, window, n_procs, area_names, sides)

                    Plot_matrix_maps(average_areas, file_dir, file_name, area_names, img_rate, unit, spikes, channel, ptp)


                print "...done."
                print " "
       
        #************************************COMPUTE STMS********************************************
        #Compute and plot Spike-Triggered-Maps (STMS): single frame representations of activity locked to spiking
        #Need to compute motifs first (see above) 
        elif sua_static_maps:
            
            for i in range(len(units)):
                print "*Processing ", file_name, " unit: ", i+1, " of ", len(units), " unit name: ", units[i]
                unit=units[i]
                channel=channels[i]
                ptp=ptps[i]
                
                #Process specific units only, otherwise process all cells from recording
                if True:
                    if unit not in cell_list: continue            #Process unit 2 only
    
                #*********** Load sorted spikes for unit
                spikes = np.loadtxt(file_dir+file_name+'/unit_'+str(unit).zfill(2)+ '_channel_' + str(channel).zfill(2) + '_ptp_'+str(ptp).zfill(3)+'.csv')
                
                #Determine if any spikes in recording window to continue:
                temp0 = np.where(np.logical_and(spikes>=img_times[0]+window, spikes<=img_times[-1]-window))[0]
                spikes_in_window = spikes[temp0]
                if len(spikes_in_window) == 0: 
                    print "Zero spikes in window - skip unit"
                    continue
                n_spikes = len(spikes_in_window)  #NB: Save n_spikes value just for current epoch for continuous recordings

                #*********** Compute static maps +/- 1sec from spike rasters
                plotting = True     #View STM maps;
                Compute_static_maps_max(img_rate, window, n_procs, main_dir, file_dir, file_name, n_pixels, unit, channel, n_spikes, ptp, stm_types, plotting, spiking_modes)
                
       
        #************************************COMPUTE LFP TRIGGERED MAPS********************************************
        #Compute and plot LFP Amplitude generated maps
        if lfp_static_maps:
            n_pixels = 256

            file_name_aligned = file_dir + file_name+'/'+file_name+'_images_aligned.npy'
            if (os.path.exists(file_name_aligned)==False):
                images_raw = Load_images(file_dir, file_name)
                images_rotated = Rotate_images(images_raw, file_dir, file_name, overwrite_shifted)
                images_aligned = Define_bregma_lambda(images_rotated, main_dir, file_dir, file_name)
            else:
                images_aligned = np.load(file_name_aligned)
            
            Compute_lpf_static_maps(file_dir, file_name, images_aligned, img_start, img_end, len_frame, img_rate, n_pixels, n_frames, img_times)


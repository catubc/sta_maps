import numpy as np
import time, math
import os.path
import multiprocessing as mp
import re
import sys
import struct, array, csv
import scipy.optimize 

from matplotlib import animation
import matplotlib.mlab as mlab

from scipy import stats
from scipy import signal
from pylab import *

from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline

from PIL import Image

from sta_utils import *

#***********************************************
np.seterr(divide='ignore', invalid='ignore')        #Ignore division by zero due to divide frames that were rotated and contain unmasked pixels

#******************************** SET LOADING DIRECTORIES DEFAULTS *************
main_dir = ''

file_dirs = []
file_names = []

#**************** LOAD ALL FILES DATA IN EXPERIMENT DIRECTORY ****************
file_dirs.append('2015-7-22/') 
file_names.append([
'2015-7-22-1',     #cortical recording          #Number of cells: 9
#'2015-7-22-4',     #subcortical recording      #Number of cells: 30
])

#Select 
cell_list = [1]

#******Processing info******
window = 3      #Window of STM motif/montages
n_procs=14      #Number of processors to use for parallel sections of code; Ensure not using more than available cores or significant slow down may occur (or crashes)

#Meta data for generating videos and STMTDS post STM processing
area_names = ['hindlimb', 'forelimb', 'barrel','retrosplenial','visual', 'motor', 'pta', 'acc'] 
sides = ['left','right']

#Spiking modes to use for computing STMs;  #Current available modes: ['burst', 'last', 'tonic', 'first']
spiking_modes = ['last', 'burst']

''' NB: for 'modes' option, need to generate an input text file with extention "_grouped_spikes.txt" containing each spiking mode in simple text format, 
        followed by spike times (in seconds) from each mode on separate lines.
    
    For example: filename:  "unit_01_channel_06_ptp_037_imagingspikes_grouped_spikes.txt".
    - line 1: "burst" 
    - line 2: "last"
    - line 3: "tonic"
    - line 4: "first"
    - line 5: "28.072, 29.431, 30.626, 30.640 ..."  <- spikes belonging to spikes in the middle of a burst
    - line 6: "14.944, 15.889, 28.089, 29.490 ..."  <- spikes belonging to last spike in a burst
    - line 7: "8.657, 9.998, 12.134, 12.835 ..."    <- spikes belonging to tonic modes
    - line 8: "13.547, 14.98, 26.120, 28.308 ..."   <- spikes belonging to first spike in a burst

    Additional modes can be defined by changing the spiking modes name in the code; or by simply re-purposing existing labels to a different class of spikes.
'''

stm_types = ["all"]     # Must insert spiking type to process data below; current options are: "all" or "modes"; 
                          # 'modes' option is defined above and requires a separate text file 


#Flags for computing sta motifs;
compute_sta_motif = True
overwrite = False                       #Overwrite existing STM data if already computed
view_sta_motif = False                  #View each unit's STM after processing


#Flags for computing STMs
compute_stm = True


#Flags for STMTD computations
compute_stmtd = False
view_stmtd = False


#Flag for making videos of STM and STMTDs
animate_images = False


#Set this flag along with the correct unit number to compute STMs for all motifs in a recording; This is used for other processing (not included here).
random_flag = False  


#****************************************************************************************
#Loop over experiments 
for dir_counter, file_dir in enumerate(file_dirs):

    for file_name in file_names[dir_counter]:
        ''' These loops search for and load spike rasters saved in .csv files. The raster filenames contain metadata such as # of spikes, channel, and PTP max.
        '''
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
            

        ''' After loading spike rasters and aligned imaging data, can proceed to compute and visualize sta-motifs, STMTDs, or STMs.
        '''
        
        #********************************COMPUTE SPIKE-TRIGGERED-MOTIFS********************************************
        #Compute and plot motifs: multi-frame representations of activity locked to spiking
        if compute_sta_motif or view_sta_motif or compute_stmtd or view_stmtd or animate_images:
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
    
                #Determine if there are spikes in recording window to continue; Required as some units are sparsely firing
                spikes = np.loadtxt(file_dir+file_name+'/unit_'+str(unit).zfill(2)+ '_channel_' + str(channel).zfill(2) + '_ptp_'+str(ptp).zfill(3)+'.csv')
                temp0 = np.where(np.logical_and(spikes>=img_times[0]+window, spikes<=img_times[-1]-window))[0]
                spikes_in_window = spikes[temp0]
                if len(spikes_in_window) == 0:  print "Zero spikes in window - skip unit"; continue
                
                n_spikes = len(spikes_in_window)  #NB: Save n_spikes value just for current epoch for continuous recordings

                #Compute and view motifs
                Compute_sta_motif(unit, channel, spikes, window, img_rate, img_times, n_pixels, images_aligned, file_dir, file_name, n_procs, 
                                  overwrite, stm_types, random_flag, spiking_modes)

                #Preview STM for each cell
                if view_sta_motif:
                    View_sta_motif(unit, main_dir, file_dir, file_name, stm_types, img_rate, spiking_modes)


                #Search Max/Min, Save Time Courses
                if compute_stmtd:
                    #Search for Min and Max pixels
                    Compute_STMTD(unit, channel, spikes, file_dir, file_name, img_rate, window, n_procs, area_names, depth, sides, stm_types, spiking_modes)

                #Plot STMTDs if computed
                if view_stmtd:
                    View_STMTD(unit, channel, spikes, window, len_frame, file_dir, file_name, area_names, sides, stm_types, spiking_modes)

                #Make animated STMTD videos
                if animate_images: 
                    Animate_images(unit, channel, window, img_rate, main_dir, file_dir, file_name, n_pixels, spikes, plot_string, n_procs, area_names, sides, depth, stm_types, spiking_modes)
                
      
        #************************************COMPUTE STMS********************************************
        #Compute and plot Spike-Triggered-Maps (STMS): single frame representations of activity locked to spiking (see Xiao et al Methods)
        #Need to compute motifs first (see above) 
        if compute_stm:
            
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
                Compute_STM(img_rate, window, n_procs, main_dir, file_dir, file_name, n_pixels, unit, channel, n_spikes, ptp, stm_types, spiking_modes)
                
       

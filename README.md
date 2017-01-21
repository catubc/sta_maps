# sta_maps
##**Spike Triggered Average Mapping Code (Python 2.7)**

Python code for accompanying manuscript (in review): "Mapping cortical mesoscopic networks of single spiking cortical or sub-cortical neurons", Xiao D, Vanni MP, Mitelut C, Chan AW, LeDue JM, Xie Y, Chen AC, Swindale, NV, and Murphy TH. (eLife; in review).

The script computes spike-triggered-average motifs (or montages), spike-triggered-maps (STMs) and spike-triggered-map temporal dynamics (STMTDs). 



###**Dependencies and Source Data**

Processing requirements are at least 32GB of ram and preferably a multi-core processor. 

There are several python module dependencies including: numpy, math, mutiplrociseing, scipy, PIL, subprocess and others. They are all standard packages that can be installed using pip or other common methods.

The code was written to load Multi-Channel-Systems ephys files that align to imaging data. We provide examples of data recorded cortically and subcortically here: [source data on dropbox](https://www.dropbox.com/sh/chet957crw41267/AADgke5NMnM__f4L4PDaK4QHa?dl=0).

The imaging data files (converted and aligned) should be downloaded and saved into the appropriate directories.

###**Instructions for Running Code**

The code is a simplified version of our scripts. It should be run using python 2.7 using the command:
python sta_maps.py

There are several option flags in the sta_maps.py file which can be set to compute and display motifs, STMs and STMTDs. Additional, video animations can be generated using the "animate_images" flag.

Assistance can be provided by email: cat@alumni.ubc.ca.


###**Matlab Code**

Alternative Matlab code with other examples is also available here [Matlab code and source data] (https://www.dropbox.com/home/nmpaper/Figure%201-source%20data%201).

# sta_maps
##**Spike Triggered Average Mapping Code (Python 2.7)**

Python code for accompanying manuscript (in review): "Mapping cortical mesoscopic networks of single spiking cortical or sub-cortical neurons", Xiao D, Vanni MP, Mitelut C, Chan AW, LeDue JM, Xie Y, Chen AC, Swindale, NV, and Murphy TH. (eLife; in review).

The script computes spike-triggered-average motifs (or montages), spike-triggered-maps (STMs) and spike-triggered-map temporal dynamics (STMTDs). 



###**Dependencies and Source Data**

Processing requirements are at least 32GB of ram and preferably a multi-core processor. There are several python module dependencies including: `numpy, math, multiprocessing, scipy, PIL, subprocess` and others. These are all common packages that can be installed using pip or other standard methods.

The code was written to load Multi-Channel-Systems ephys files that align to imaging data. The imaging data files (converted and aligned) should be downloaded and saved into the appropriate directories. We provide examples of cortical and subcortical recordings: [source data](https://www.dropbox.com/sh/chet957crw41267/AADgke5NMnM__f4L4PDaK4QHa?dl=0).


###**Instructions for Running Code**

The code is a simplified version of our scripts. It should be run using python 2.7 using the command:

`python main.py`

There are several option flags in the main.py file which can be set to compute and display: motifs, STMs and STMTDs. 

Assistance with module installation and running can be provided on github via the "Issues" tab.


###**Matlab Code**

Alternative Matlab code with additional examples (and source data for Figure 1 from above reference publication) is also available here: [Matlab code and source data] (https://www.dropbox.com/sh/l7nk4f6ula0mmte/AAC2YkMT7Z6AE8VAUxDtduiya?dl=0).

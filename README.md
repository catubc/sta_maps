# sta_maps
Spike Triggered Average Mapping Code (python)

Python code for accompanying manuscript (in review): "Mapping cortical mesoscopic networks of single spiking cortical or sub-cortical neurons", Xiao D, Vanni MP, Mitelut C, Chan AW, LeDue JM, Xie Y, Chen AC, Swindale, NV, and Murphy TH. (eLife; in review).

The script computes spike-triggered-average motifs (or montages), spike-triggered-maps (STMs) and spike-triggered-map temporal dynamics (STMTDs). 

Processing requirements are at least 32GB of ram and preferably a multi-core processor. 

There are several python module dependencies including: numpy, math, mutiplrociseing, scipy, PIL, subprocess and others. They are all standard packages that can be installed using pip or other common methods.

The code was written to load Multi-Channel-Systems ephys files that align to imaging data. We provide examples of data recorded cortically and subcortically here: https://www.dropbox.com/sh/chet957crw41267/AADgke5NMnM__f4L4PDaK4QHa?dl=0


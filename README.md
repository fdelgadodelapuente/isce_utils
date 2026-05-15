ISCE2 MATLAB/Python/GMT utilities for loading interferograms in the JPL RMG file format, mask the data, remove ramps with or without an empirical correlation between phase and topography, and for exporting to GMT grids. These utilities require to requiere to download [grdwrite2](https://www.mathworks.com/matlabcentral/fileexchange/26290-grdwrite2) and [grdread2](https://www.mathworks.com/matlabcentral/fileexchange/25683-grdread2). They include ALOS-2 interferograms as sample data sets.  

alos2_wd1_deramp_example.zip is the swath 4 of an ALOS-2 ScanSAR interferogram in range azimuth coordinates over Aniakchak Crater in Alaska. For removing the ramp uncomment lines 93-95 in deramp.py. These lines will mask unwrapping errors. The code will generate the following image and export the corrected interferogram.

<img width="1920" height="1440" alt="isce_deramp" src="https://github.com/user-attachments/assets/450df6ee-1b20-40d3-a15f-e61e21cc3bec" />


For exporting an ISCE interferogram to GMT, unzip the alos2_sm3_220711_230612.zip file, first run isce2gmt.m in MATLAB and then run isce2gmt.sh in the terminal (it is C Shell despite the .sh extension).  The code will generate the following image.

<img width="1068" height="869" alt="aniakchak_alos2_sm3_2022_2023" src="https://github.com/user-attachments/assets/beb5df29-f929-4857-ba8a-d25756dc5c4a" />


If you use these codes, please cite the following paper.

Delgado, F. (2026). Abnormally large magma flux does not lead to eruption in subduction zone calderas: The 2022–2023 episode of uplift of Aniakchak Crater (Aleutians). Geophysical Research Letters, 53, e2025GL117786. https://doi.org/10.1029/2025GL117786



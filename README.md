ISCE2 MATLAB/Python/GMT utilities for loading interferograms in the JPL RMG file format, mask the data, remove ramps with or without an empirical correlation between phase and topography, and for exporting to GMT grids. They include ALOS-2 interferograms as sample data sets.

These utilities require to requiere to download [grdwrite2](https://www.mathworks.com/matlabcentral/fileexchange/26290-grdwrite2)) and [grdwrite2](https://www.mathworks.com/matlabcentral/fileexchange/25683-grdread2)) 






alos2_wd1_deramp_example.zip is the swath 4 of an ALOS-2 ScanSAR interfeorgram in range azimuth coordinates over Aniakchak Crater in Alaska. For removing the ramp uncomment lines 93-95 in deramp.py. These lines will mask unwrapping errors. The code will generate the following image and export the corrected interferogram

<img width="1920" height="1440" alt="isce_deramp" src="https://github.com/user-attachments/assets/450df6ee-1b20-40d3-a15f-e61e21cc3bec" />




ISCE2 MATLAB/Python/GMT utilities for loading interferograms in the JPL RMG file format, mask the data, remove ramps with or without an empirical correlation between phase and topography, and for exporting to GMT grids. These utilities require to requiere to download [grdwrite2](https://www.mathworks.com/matlabcentral/fileexchange/26290-grdwrite2) and [grdread2](https://www.mathworks.com/matlabcentral/fileexchange/25683-grdread2). They include ALOS-2 interferograms as sample data sets.  

alos2_wd1_deramp_example.zip is the swath 4 of an ALOS-2 ScanSAR interferogram in range azimuth coordinates over Aniakchak Crater in Alaska. For removing the ramp uncomment lines 93-95 in deramp.py. These lines will mask unwrapping errors. The code will generate the following image and export the corrected interferogram.

<img width="1920" height="1440" alt="isce_deramp" src="https://github.com/user-attachments/assets/450df6ee-1b20-40d3-a15f-e61e21cc3bec" />


For exporting an ISCE interferogram to GMT, unzip the alos2_sm3_220711_230612.zip file (ALOS-2 stripmap SM3 over Aniakchak Crater) , run isce2gmt.m in MATLAB and then run isce2gmt.sh in the terminal (it is C Shell despite the .sh extension).  The code will generate the following image.

<img width="534" height="435" alt="aniakchak_alos2_sm3_2022_2023" src="https://github.com/user-attachments/assets/beb5df29-f929-4857-ba8a-d25756dc5c4a" />


If you use these codes, please cite the following paper.

Delgado, F. (2026). Abnormally large magma flux does not lead to eruption in subduction zone calderas: The 2022–2023 episode of uplift of Aniakchak Crater (Aleutians). Geophysical Research Letters, 53, e2025GL117786. https://doi.org/10.1029/2025GL117786




For compiling the software in macOS 
```
### Feb 22 2019 MacMini2014/MacBookAir2015, High Sierra and Mojave, gcc7
### Nov XY 2021 MacMini2014 Monterey, python37 and gcc11 
### Feb 06 2025 MacBookAir2015 Monterey, python312 and gcc13
### Feb 07 2025 MacBookAir2015 Monterey, python39 and gcc11 after python312 failure in running stripmapApp. Default compiler in the mp system is gcc13, though
### Aug 07 2025 MacMini M4 Pro Sequoia, python313 and gcc13
### May 28 2026 MacBook Neo Tahoe, python313 and gcc13


#### first reinstall macports for the new OSX version, then
port -qv installed > myports.txt ### backup existing ports
port echo requested | cut -d ' ' -f 1 > requested.txt ### backup another thing
sudo port -f uninstall installed ### remove old ports
sudo rm -rf /opt/local/var/macports/build/*  ### remove old ports
curl --location --remote-name 
https://github.com/macports/macports-contrib/raw/reference/restore_ports/restore_ports.tcl
chmod +x restore_ports.tcl

#if starting from a new MacOS, install macports and then self update
sudo port selfupdate

### install Xcode, then install the Xcode tools with
xcode-select --install

### ISCE is written in FORTRAN, C, C++. Use the GNU compiler gcc 
sudo port install gcc13
sudo port select gcc mp-gcc13    
port select --list gcc ## check the versions of gcc installed
## restart terminal after this

### Install python3 and many libraries. Never use the default Python included in macOS
sudo port install python313
sudo port select python3 python313
sudo ln -s /opt/local/Library/Frameworks/Python.framework/Versions/3.13/include/python3.13m /opt/local/include/python3.13m  #not neede for python3.13 in MacMini Neo

sudo port install xorg-libXt +flat_namespace   #for mdx
sudo port install freetype tiff openmotif      #for mdx
sudo port install fftw-3 +gcc13   
sudo port install fftw-3-single +gcc13    
sudo port install hdf5 +gcc13    #FAILED, then skipped to install h5py Feb 2025
sudo port install hdfeos5 h5utils #FAILED, then skipped to install h5py Feb 2025
sudo port install py313-numpy +gcc13
sudo port install py313-scipy +gcc13 #FAILED, installed with +gcc flag
sudo port install py313-matplotlib +cairo
sudo port install py313-pandas    ###not sure if I really needed
sudo port install py313-cython
ln  -s /opt/local/bin/cython-3.13 /opt/local/bin/cython3
cython3 -V #### check cython3 version  > 0.28
sudo port install py313-h5py  #python binding for h5py
sudo port install py313-matplotlib-basemap #flag not existent anymore with 3.12 version
####sudo port install opencv +python36 #for ionospheric correction only
#### for Monterey update to opencv3, still called cv2 in python, though
#sudo port install opencv3 +python313
#sudo port install py313-opencv3   
#### for Monterey and python3.12 update to opencv4, still loaded as import cv2, though.
sudo port install opencv4 +python313
sudo port install py313-opencv4     #python binding for opencv4, it doesnt exist for opencv3 in the MacBookNeo

#### ---
sudo port install py313-ipython 
sudo port select --set ipython3 py313-ipython
sudo port install py313-jupyter 
sudo port install hdf4 hdfeos
######sudo port install postgresql95 postgresql95-server    #this wasn't installed
#sudo port install gdal +curl +expat +geos +hdf4 +hdf5 +netcdf  +openjpeg +postgresql95 +sqlite3
sudo port install gdal +curl +expat +geos +gdal-hdf5 +netcdf  +openjpeg +postgresql95 +sqlite3
#2026/05/28 MacNeo. Error: The '+netcdf' variant has been removed and replaced by the 'gdal-netcdf' subport. I then used
sudo port install gdal +curl +expat +geos +gdal-hdf5   +openjpeg +postgresql95 +sqlite3
export GDAL_DATA=/opt/local/share/gdal
sudo port install scons
sudo port install py313-gdal  #python binding for gdal

sudo port install gmt5 +fftw3 #only for gmt plots
sudo port install py313-rasterio  #for loading and exporting ifgs in python
sudo port install ImageMagick #for exporting interferograms to .kml

```


Create /Applicatons/isce/SConfigISCE file for scons
```
PRJ_SCONS_BUILD   = /Applications/isce/isce2-2.6.4/build/isce
PRJ_SCONS_INSTALL = /Applications/isce/isce2-2.6.4/install/isce


LIBPATH = /opt/local/lib
#the last path in CPP is new for autoRIFT in 2.4.x version
CPPPATH = 
/opt/local/Library/Frameworks/Python.framework/Versions/3.13/include/python3.13 
/opt/local/include /opt/local/include/opencv4  /opt/local/lib/opencv4
/opt/local/Library/Frameworks/Python.framework/Versions/3.13/lib/python3.13/site-packages/numpy/core/include
FORTRANPATH = /opt/local/include
FORTRAN = /opt/local/bin/gfortran
CC = /opt/local/bin/gcc
CXX = /opt/local/bin/g++

#libraries needed for mdx display utility
MOTIFLIBPATH = /opt/local/lib       # path to libXm.dylib
X11LIBPATH = /opt/local/lib         # path to libXt.dylib
MOTIFINCPATH = /opt/local/include   # path to location of the Xm
                                    # directory with various include files (.h)
X11INCPATH = /opt/local/include     # path to location of the X11 directory
                                    # with various include files

# turn off CUDA code on this computer
ENABLE_CUDA = FALSE
```

Now install it in Applicatons/isce/isce2-2.6.4
```
rm -rf config.log .sconfig.dblite .sconf_temp .sconsign.dblite; SCONS_CONFIG_DIR=/Applications/insar_software/isce scons install  
```

Source it with shell script called insar.sh. Inp is the software version (2.2.0, 2.5.1, 2.6.4, etc)
```
#!/bin/sh

inp="$1"
echo "Loading ISCE $1"
 
export PYTHONPATH=/Applications/insar_software/isce/isce-$1/install:$PYTHONPATH
export PATH=/Applications/insar_software/isce/isce-$1/install/isce/bin:$PATH
export PATH=/Applications/insar_software/isce/isce-$1/install/isce/applications:$PATH
export ISCE_HOME=/Applications/insar_software/isce/isce-$1/install/isce
 
export PATH=/Applications/insar_software/isce/isce-$1/contrib/stack/stripmapStack:$PATH
#export PATH=/Applications/isce/isce$1/contrib/stack/topsStack:$PATH
  
export GDAL_DATA=/opt/local/share/gdal
```




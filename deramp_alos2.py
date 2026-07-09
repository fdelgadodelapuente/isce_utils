#!/usr/bin/env python

# ramp removal after phase unwrapping 
# Python3 translation of very old load_isce.m
# Francisco Delgado, May 2025

from osgeo import gdal
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import isce
import warnings
import sys
import os
import glob


warnings.filterwarnings("ignore")
plt.rcParams['font.family'] = 'Helvetica'

path_unw=glob.glob("insar/filt_*msk.unw")[0]
##print(path_unw[1:-4])

####cmd = 'grep snaphu *.xml | wc -l > kk'  # check snaphu or icu
###os.system(cmd)
###file = open('kk','r')
###os.system('rm kk')
###f = int(file.read())
###file.close()

###if f == 1:
#    print('Interferogram unwrapped with SNAPHU_MCF, masking with filt_topophase.unw.conncomp')
#    path_con='interferogram/filt_topophase.unw.conncomp'
###else:
   ### print('Interferogram unwrapped with ICU, masking with filt_topophase.conncomp')
   ### path_con='interferogram/filt_topophase.conncomp'
    
cax_min=-2*2*np.pi
cax_max=-cax_min

cax_min=60
cax_max=90

####box = [50,570,386,767] #x1,x2,y1,y2

## backup original file
path_unwb = path_unw + '.orig'
tmp = os.path.exists(path_unwb)

if tmp == False:
    print('Backup original file')
    cmd = 'cp ' +path_unw+ ' ' +path_unw+'.orig'
    os.system(cmd)
    cmd = 'cp ' +path_unw+'.xml'+ ' ' +path_unw+'.orig.xml'
    os.system(cmd)
    cmd = 'cp ' +path_unw+'.vrt'+ ' ' +path_unw+'.orig.vrt'  
    os.system(cmd)

cmd = 'cp ' +path_unwb+'.vrt'+ ' ' +path_unw+'.vrt'  
os.system(cmd)
cmd = 'cp ' +path_unwb+'.xml'+ ' ' +path_unw+'.xml'  
os.system(cmd)

# ---- reemplazo rasterio -> gdal ----
ds = gdal.Open(path_unwb, gdal.GA_ReadOnly)
mag = ds.GetRasterBand(1).ReadAsArray()
unw = ds.GetRasterBand(2).ReadAsArray()
ds = None

unw = np.flipud(unw)
unw[unw == 0] = np.nan

ny, nx = np.shape(unw)
unw[0:1500,0:nx]=0 #mask unw error bottom
unw[4050:ny,0:nx]=0 #mask unw error top
unw[ny-495:ny-375,800:nx]=0 #mask unw error island


###ds = gdal.Open(path_con, gdal.GA_ReadOnly)
###con = ds.GetRasterBand(1).ReadAsArray()
###ds = None

###con = np.flipud(con)
###unw[con == 0] = np.nan
#unw[mag == 0] = np.nan
#unw[mag<.009] = np.nan

phs_wr = np.arctan2(np.sin(unw), np.cos(unw))

[ny,nx]=np.shape(unw)

x = np.arange(nx)
y = np.arange(ny)

Xv, Yv = np.meshgrid(x, y)

unw[unw==0] = np.nan

# estimate and remove ramp by linear least squares
d = unw.reshape(-1,1)
g = ~np.isnan(d)

xr = Xv.reshape(-1,1)
yr = Yv.reshape(-1,1)
xryr = np.multiply(xr,yr)
dcr = xr*0+1

G = np.transpose([xr,yr,dcr])
Gg = np.transpose([xr[g],yr[g],dcr[g]])
G = np.transpose([xr,yr,xryr,dcr])
Gg = np.transpose([xr[g],yr[g],xryr[g],dcr[g]])
dg = d[g]

m = np.linalg.lstsq(Gg, dg)[0]
ramp = np.matmul(G, m)
ramp = ramp.reshape(ny,nx)

unw_flat = unw - ramp
ramp[np.isnan(unw) == 1] = np.nan

dxplot = 1
wrapped = (unw + np.pi) % (2 * np.pi) - np.pi

# Figures
fig = plt.figure(1)
fig.subplots_adjust(wspace=0.4)

ax = fig.add_subplot(2,3,1)
cax = ax.pcolormesh(wrapped[::dxplot,::dxplot], vmin = -np.pi , vmax = np.pi, cmap = 'jet')
ax.set_title("wrapped interferogram [rads]", fontsize=6)
cbar = fig.colorbar(cax, orientation='vertical',fraction=0.03, pad=0.04)
cbar.ax.tick_params(labelsize=6) 
ax.set_aspect('equal', adjustable='box')
ax.tick_params(axis='both', labelsize=6)

ax = fig.add_subplot(2,3,2)
cax = ax.pcolormesh(unw[::dxplot,::dxplot], vmin = cax_min , vmax = cax_max, cmap = 'jet')
ax.set_title("Unwrapped interferogram [rads]", fontsize=6)
cbar = fig.colorbar(cax, orientation='vertical',fraction=0.03, pad=0.04)
cbar.ax.tick_params(labelsize=6) 
ax.set_aspect('equal', adjustable='box')
ax.tick_params(axis='both', labelsize=6)

ax = fig.add_subplot(2,3,3)
cax = ax.pcolormesh(ramp[::dxplot,::dxplot], vmin = cax_min , vmax = cax_max, cmap = 'jet')
ax.set_title("Least-squares ramp [rads]", fontsize=6)
cbar = fig.colorbar(cax, orientation='vertical',fraction=0.03, pad=0.04)
cbar.ax.tick_params(labelsize=6) 
ax.set_aspect('equal', adjustable='box')
ax.tick_params(axis='both', labelsize=6)

ax = fig.add_subplot(2,3,5)
cax = ax.pcolormesh(unw_flat[::dxplot,::dxplot], vmin = -3*np.pi , vmax = 3*np.pi, cmap = 'jet')
ax.set_title("Unwrapped flattened interferogram [rads]", fontsize=6)
cbar = fig.colorbar(cax, orientation='vertical',fraction=0.03, pad=0.04)
cbar.ax.tick_params(labelsize=6) 
ax.set_aspect('equal', adjustable='box')
ax.tick_params(axis='both', labelsize=6)

plt.savefig("isce_deramp.png", dpi=300)

unw_flat = np.flipud(unw_flat)
unw_flat[np.isnan(unw_flat)==1] = 0

# ---- reemplazo rasterio -> gdal (write) ----
ds = gdal.Open(path_unw, gdal.GA_Update)
ds.GetRasterBand(1).WriteArray(mag)
ds.GetRasterBand(2).WriteArray(unw_flat)
ds = None

cmd = 'cp ' +path_unw+'.orig.xml'+ ' ' +path_unw+'.xml'
os.system(cmd)
cmd = 'cp ' +path_unw+'.orig.vrt'+ ' ' +path_unw+'.vrt'  
os.system(cmd)

os.system('xdg-open isce_deramp.png')

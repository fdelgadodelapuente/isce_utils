#!/usr/bin/env python

#ramp removal after phase unwrapping 
#Python3 translation of very old load_isce.m
#Francisco Delgado, May 2025
#Updated May 2026 for stripmapApp, topsApp and alos2App file formats


from osgeo import gdal
import rasterio
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
#matplotlib.pyplot.set_loglevel (level = 'warning')
#np.set_printoptions(threshold=sys.maxsize)

path_unw=glob.glob("*/*.unw")
path_con=glob.glob("*/filt*conncomp")[0] #only for SNAPHU_MCF, automatically masked by ICU
if len(path_unw)==2:   #for ALOS-2
    path_unw=glob.glob("*/*.unw")[1]
else:
    path_unw=glob.glob("*/*.unw")[0]
 
box = [50,570,386,767] #x1,x2,y1,y2

##backup original file
path_unwb = path_unw + '.orig'
tmp = os.path.exists(path_unwb)
#######print(tmp)
if  tmp == False:
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

with rasterio.open(path_unwb) as src:
    profile = src.profile.copy()
    mag = src.read(1)
    unw = src.read(2)
unw = np.flipud(unw)
unw[unw == 0] = np.nan

if path_con != []:
    with rasterio.open(path_con) as src:
        con = src.read(1)
    con = np.flipud(con)
    unw[con == 0] = np.nan

#unw[mag == 0] = np.nan
#unw[mag<.009] = np.nan


ds = gdal.Open(path_unw, gdal.GA_ReadOnly)
mag = ds.GetRasterBand(1).ReadAsArray()

phs_wr = np.arctan2(np.sin(unw), np.cos(unw))

#box = unw[boxx1:boxx2,boxy1:boxy2]
#unw = unw - np.nanmean(box)


[ny,nx]=np.shape(unw)

x = np.arange(nx)
y = np.arange(ny)
nDim = nx*ny;

Xv, Yv = np.meshgrid(x, y)

unw[unw==0] = np.nan
####unw[:,1:300] = np.nan; %ALOS-4 fixed PRF
###unw[0:400,0:nx]=np.nan #mask unw
###unw[0:2500,0:nx]=np.nan #mask unw Kilauea WD1

###aniakchak 2022-2023 WD1 1 swath test
##unw[0:1500,0:nx]=np.nan #mask unw error bottom
##unw[4050:ny,0:nx]=np.nan #mask unw error top
##unw[ny-495:ny-375,800:nx]=np.nan #mask unw error island

#estimate and remove ramp by linear least squares
d = unw
d = d.reshape(-1,1)
g = ~np.isnan(d)
xr = Xv.reshape(-1,1)
yr = Yv.reshape(-1,1)
xryr = np.multiply(xr,yr)
dcr = xr*0+1
G = np.transpose( [xr,yr,dcr] )
Gg = np.transpose( [xr[g],yr[g],dcr[g]] )
#G = np.transpose( [xr,yr,xryr,dcr] )
#Gg = np.transpose( [xr[g],yr[g],xryr[g],dcr[g]] )
dg = d[g]
m = np.linalg.lstsq(Gg, dg)[0]; #first element is the structure with the inverted model parameters
ramp = np.matmul(G, m) #estimate ramp
ramp = ramp.reshape(ny,nx)
unw_flat = unw-ramp
ramp[np.isnan(unw) == 1] = np.nan

dxplot = 1

wrapped = (unw + np.pi) % (2 * np.pi) - np.pi


#Figures

cax_min=0.95 * np.nanmin(unw)
cax_max=0.95 * np.nanmax(unw)

fig = plt.figure(1) ##plt.figure(figsize=(14,12))
fig.subplots_adjust(wspace=0.4)
ax = fig.add_subplot(2,3,1)
cax = ax.pcolormesh(wrapped[::dxplot,::dxplot], vmin = -np.pi , vmax = np.pi, cmap = 'jet')
ax.set_title("wrapped interferogram [rads]", fontsize=6)
cbar = fig.colorbar(cax, orientation='vertical',fraction=0.03, pad=0.04)
cbar.ax.tick_params(labelsize=6) 
#plt.colorbar(shrink=0.5)
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
#plt.xlim(xu_min,xu_max)
#plt.ylim(yu_min,yu_max)
ax.tick_params(axis='both', labelsize=6)

ax = fig.add_subplot(2,3,2)
cax = ax.pcolormesh(unw[::dxplot,::dxplot], vmin = cax_min , vmax = cax_max, cmap = 'jet')
ax.set_title("Unwrapped interferogram [rads]", fontsize=6)
cbar = fig.colorbar(cax, orientation='vertical',fraction=0.03, pad=0.04)
cbar.ax.tick_params(labelsize=6) 
#plt.colorbar(shrink=0.5)
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
#plt.xlim(xu_min,xu_max)
#plt.ylim(yu_min,yu_max)
ax.tick_params(axis='both', labelsize=6)

ax = fig.add_subplot(2,3,3)
cax = ax.pcolormesh(ramp[::dxplot,::dxplot], vmin = cax_min , vmax = cax_max, cmap = 'jet')
ax.set_title("Least-squares ramp [rads]", fontsize=6)
cbar = fig.colorbar(cax, orientation='vertical',fraction=0.03, pad=0.04)
cbar.ax.tick_params(labelsize=6) 
#plt.colorbar(shrink=0.5)
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
#plt.xlim(xu_min,xu_max)
#plt.ylim(yu_min,yu_max)
ax.tick_params(axis='both', labelsize=6)

ax = fig.add_subplot(2,3,5)
cax = ax.pcolormesh(unw_flat[::dxplot,::dxplot], vmin = 0.8*np.nanmin(unw_flat) , vmax = 0.8*np.nanmax(unw_flat), cmap = 'jet')
ax.set_title("Unwrapped flattened interferogram [rads]", fontsize=6)
cbar = fig.colorbar(cax, orientation='vertical',fraction=0.03, pad=0.04)
cbar.ax.tick_params(labelsize=6) 
#plt.colorbar(shrink=0.5)
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
#plt.xlim(xu_min,xu_max)
#plt.ylim(yu_min,yu_max)
ax.tick_params(axis='both', labelsize=6)

plt.savefig("isce_deramp.png", dpi=300) 

#plt.show()

unw_flat = np.flipud(unw_flat)
unw_flat[np.isnan(unw_flat)==1] = 0
profile['scheme'] = 'BIL'
with rasterio.open(path_unw, 'w', **profile) as dst:
    dst.write_band(1, mag)
    dst.write_band(2, unw_flat)

cmd = 'cp ' +path_unw+'.orig.xml'+ ' ' +path_unw+'.xml'
os.system(cmd)
cmd = 'cp ' +path_unw+'.orig.vrt'+ ' ' +path_unw+'.vrt'  
os.system(cmd)

#fig = plt.figure(2)
#cax = ax.pcolormesh(unw_flat[::dxplot], vmin = cax_min , vmax = cax_max, cmap = 'jet')
#plt.show()

###cmd ='fixImageXml.py -f -i interferogram/filt_topophase.unw'
###os.system(cmd)
os.system('open isce_deramp.png')
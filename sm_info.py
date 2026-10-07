#!/usr/bin/python3

#%posting calculation ISCE from StripmapProc.py, also from jupyter notebook
#https://github.com/isce-framework/isce2-docs/blob/master/Notebooks/UNAVCO_2020/Stripmap/stripmapApp.ipynb
#https://github.com/isceplus/2026-isceplus/blob/main/S03_Stripmap_data_processing_with_stripmapApp/stripmapApp.ipynb

#F.D. 2026/10/06-07 updated to plot swath dimensions

import isce
import isceobj
import isceobj.StripmapProc.StripmapProc as St
from isceobj.Planet.Planet import Planet
from isceobj.Planet.AstronomicalHandbook import c as SPEED_OF_LIGHT
import numpy as np
import logging
logging.getLogger("matplotlib").setLevel(logging.WARNING)
import matplotlib.pyplot as plt
import os
import sys
from pprint import pprint
#import warnings
#warnings.filterwarnings("ignore", module="matplotlib")
#warnings.filterwarnings("ignore")

stObj = St()
stObj.configure()
frame = stObj.loadProduct("reference_slc.xml")
mission = frame.catalog["instrument"]["platform"]["mission"]
print("Mission = {0}".format(mission))
print("Speed of light = {0} m/s".format(SPEED_OF_LIGHT))
print("Radar Wavelength = {0} m".format(frame.radarWavelegth))
print("Slant Range Pixel Size = {0} m".format(frame.instrument.rangePixelSize))
print("Radar frequency = {0} GHz".format(frame.instrument.radarFrequency/1e9))
print("Incidence Angle = {0}º".format(frame.instrument.incidenceAngle))
print("Chirp Slope = {0} ".format(frame.instrument.chirpSlope))
print("Range Sampling Rate = {0} MHz".format(frame.instrument.rangeSamplingRate/1e6))
print("Pulse Repetition Frequency = {0} Hz".format(frame.instrument.PRF))
print("Difference Slant Range Pixel Size = {0} m".format(frame.instrument.rangePixelSize - SPEED_OF_LIGHT/(2*frame.instrument.rangeSamplingRate)))
#print ('')
#print((frame.instrument))
#print((frame.instrument.__dict__))
#pprint(frame.instrument.__dict__, sort_dicts=True)
#print(list(frame.catalog.keys()))
#print ('')

#For azimuth pixel size we need to multiply azimuth time interval by the platform velocity along the track
t_mid = frame.sensingMid # the acquisition time at the middle of the scene
t_stop=frame.sensingStop
t_start=frame.sensingStart
st_mid=frame.orbit.interpolateOrbit(t_mid) #get the orbit for t_mid
st_stop=frame.orbit.interpolateOrbit(t_stop)
st_start=frame.orbit.interpolateOrbit(t_start)
Vs = st_mid.getScalarVelocity() # platform velocity
vels = [st_start.getScalarVelocity(), st_mid.getScalarVelocity(), st_stop.getScalarVelocity()]
print('Vels',vels)
prf = frame.instrument.PRF # pulse repitition frequency
#print('PRF',prf,' Hz')
ATI = 1.0/prf #Azimuth time interval 
az_pixel_size = ATI*Vs #Azimuth Pixel size
print("Azimuth Pixel Size = {0} m".format(az_pixel_size))
##print(vels[0]*ATI, vels[1]*ATI, vels[2]*ATI)


#Range Pixel size
r0 = frame.startingRange #near range
rmax = frame.getFarRange() #far range
rng =(r0+rmax)/2 #mid range

print("Near range slant range = {0} m".format(r0))
print("Mid range slant range = {0} m".format(rng))
print("Far range slant range = {0} m".format(rmax))

elp = Planet(pname='Earth').ellipsoid
tmid = frame.sensingMid

sv = frame.orbit.interpolateOrbit( tmid, method='hermite') #.getPosition()
llh = elp.xyz_to_llh(sv.getPosition())
print(sv)
print(llh)

hdg = frame.orbit.getENUHeading(tmid)
elp.setSCH(llh[0], llh[1], hdg)
sch, vsch = elp.xyzdot_to_schdot(sv.getPosition(), sv.getVelocity()) #position and velocity in SCH

print('State vectors')
print('SCH',sch)
print('VSCH',vsch)


Re = elp.pegRadCur
print('Re=', Re , 'km')
H = sch[2] #elevation in SCH
print('H=', H , 'km')
cos_beta_e = (Re**2 + (Re + H)**2 -rng**2)/(2*Re*(Re+H)) #Cosine theorem
sin_bet_e = np.sqrt(1 - cos_beta_e**2)
sin_theta_i = sin_bet_e*(Re + H)/rng #sine theorem
print("incidence angle at the middle of the swath: ", np.arcsin(sin_theta_i)*180.0/np.pi)
groundRangeRes = frame.instrument.rangePixelSize/sin_theta_i
print("Ground range pixel size: {0} m ".format(groundRangeRes))
pix_ratio=groundRangeRes/az_pixel_size
print("Pixel ratio: {0} ".format(pix_ratio))

# the azimuth pixel size calculated in Example 1 was at the orbit altitude
# here we do an approximate calculation for the pixel spacing on the ground
az_pixel_ground = az_pixel_size*Re/(Re+H)
print("Azimuth Pixel Size on ground = {0} m".format(az_pixel_ground))
pix_ratio=groundRangeRes/az_pixel_ground
print("Pixel ratio ground: {0} ".format(pix_ratio))



#calculate the viewing geometry
img = isceobj.createSlcImage()
img.load("reference_slc/reference.slc.xml")

n  = r0
p  = frame.instrument.rangePixelSize
nx = img.getWidth()
f = rmax

h=0; #local elevation
Re2 = Re + h;
Rs  = Re + H;
r_n = n
r_f = rmax
r_f2 = r_n + p*(nx-1)
print("Difference Far Field Slant Range = {0} km".format((r_f2 - r_f)/1e3))

if mission =='ALOS4': #ALOS-4 parsing code has a bug on the Far Field Slant Range
    r_f=r_f2  

psi_nr = np.arccos(  (Rs**2 + Re2**2 -r_n**2)/(2*Rs*Re2)  ); #near range
psi_fr = np.arccos(  (Rs**2 + Re2**2 -r_f**2)/(2*Rs*Re2)  ); #far range
S = Re2*(psi_fr-psi_nr)/1e3; #swath width
print("Swath Width on ground = {0} km".format(S))


##plot
npoints =300;
thetas = np.linspace(psi_nr,psi_fr,npoints); #angle arc
x = Re2*np.sin(thetas)/1000; #Swath Tierra, basic trigonometry
y = Re2*(np.cos(thetas)-1)/1000; #Swath Tierra, basic trigonometry
x_mid = Re2*np.sin(thetas[npoints//2])/1000; #%Swath Tierra mid
y_mid = Re2*(np.cos(thetas[npoints//2])-1)/1000; #Swath mid

thetas_E = np.linspace(0,.12,360);
E_surf = [Re/1e3*np.sin(thetas_E), Re/1e3*np.cos(thetas_E)-Re/1e3]; #Displaced Earth's surface, the top is zero, average elevation

plt.plot(E_surf[0],E_surf[1])
plt.plot(0,H/1e3,color='black', marker='s', markersize=5);
plt.plot([0,0],[0,np.max(H)/1e3+50],'--',color='black') #Nadir line
plt.plot([0,x[-1]], [H/1e3,y[-1]],color='black')
plt.plot([0,x_mid], [H/1e3,y_mid], '--',  color='black')
plt.plot([0,x[0]], [H/1e3,y[0]], color='black')
plt.plot(x,y, color='red', linewidth=4.5)
plt.ylabel('Elevation [km]')
plt.xlabel('Distance along Earth surface [km]')
plt.text(x_mid-20,-80,f'{S:.1f} km' , fontsize=14)
plt.text(25,200,f'H = {H/1e3:.0f} km' , fontsize=10)
plt.axis('equal')
plt.text(500,650,f'near  = {r_n/1e3:.0f} km' , fontsize=10)
plt.text(500,600,f'far = {r_f/1e3:.0f} km' , fontsize=10)
#plt.show()
plt.savefig('swath.pdf', bbox_inches='tight')

if sys.platform == "linux" or sys.platform == "linux2":
    print("Running on Linux")
    os.system('xdg-open swath.pdf')
elif sys.platform == "darwin":
    print("Running on macOS")
    os.system('open swath.pdf')
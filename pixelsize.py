#!/usr/bin/python3

#%posting calculation ISCE from StripmapProc.py, also from jupyter notebook
#https://github.com/isce-framework/isce2-docs/blob/master/Notebooks/UNAVCO_2020/Stripmap/stripmapApp.ipynb

import isce
import isceobj
import isceobj.StripmapProc.StripmapProc as St
from isceobj.Planet.Planet import Planet
import numpy as np

stObj = St()
stObj.configure()
frame = stObj.loadProduct("reference_slc.xml")
#frame = stObj.loadProduct("20221106.slc_slc.xml")
print("Wavelength = {0} m".format(frame.radarWavelegth))
print("Slant Range Pixel Size = {0} m".format(frame.instrument.rangePixelSize))
print ('')
print((frame.instrument))
print ('')

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
print('PRF',prf,' Hz')
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

hdg = frame.orbit.getENUHeading(tmid)
elp.setSCH(llh[0], llh[1], hdg)
sch, vsch = elp.xyzdot_to_schdot(sv.getPosition(), sv.getVelocity()) #position and velocity in SCH

print('State vectors')
print('SCH',sch)
print('VSCH',vsch)


Re = elp.pegRadCur
H = sch[2] #elevation in SCH
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
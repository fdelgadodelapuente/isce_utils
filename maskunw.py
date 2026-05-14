#!/usr/bin/env python3

#Requires rasterio library
#OSX sudo port install py36-rasterio
#Ubuntu sudo apt-get install rasterio

import os
import numpy as np
import rasterio

corth = 0.1
f1 = os.popen("ls f*snaphu.unw").read()
f1=f1[:-1] #remove last backslash
f2 = os.popen("ls f*.cor").read()
f2=f2[:-1] #remove last backslash
f3 = os.popen("ls f*snaphu.unw.conncomp").read()
f3 = f3[:-1]
f4 = f1[:-8] + '.masked.unw'
f4 = 'filt_topophase.unw'

##backup original file

f1b = f1 + '.orig'
tmp = os.path.exists(f1b)
#######print(tmp)
if  tmp == False:
    print('Backup original file')
    cmd = 'cp ' +f1+ ' ' +f1+'.orig'
    os.system(cmd)
    cmd = 'cp ' +f1+'.xml'+ ' ' +f1+'.orig.xml'
    os.system(cmd)
    cmd = 'cp ' +f1+'.vrt'+ ' ' +f1+'.orig.vrt'  
    os.system(cmd)

with rasterio.open(f1b) as src:
    profile = src.profile.copy()
    amp = src.read(1)
    phs = src.read(2)
    
with rasterio.open(f2) as src:  #same as phsig.cor, no real coherence file
    cor = src.read(1)
#with rasterio.open('topophase.cor') as src:
#    ampc = src.read(1)
#    cort = src.read(2)
    
with rasterio.open(f3) as src:
    concomp = src.read(1)
    
###cor=cort  #switch topophase coherence instead of phsig to mask
phs[(amp<1e7)] = 0 #remove band to the right 
phsm = phs*0;   

cc=[1,3,2] ##
nc=np.size(cc)
##phsm[(concomp==1) & (cor>corth)] = phs[(concomp==1) & (cor>corth)]
for x in cc:
  print('Adding connected component ' + str(x))
  phsm[(concomp==x) & (cor>corth)] = phs[(concomp==x) & (cor>corth)]
  phsm[(concomp==4) & (cor>corth)] = phs[(concomp==4) & (cor>corth)] + 2*np.pi #manually check how many cycles of phase


profile['scheme'] = 'BIL'
with rasterio.open(f4, 'w', **profile) as dst:
      dst.write_band(1, amp)
      dst.write_band(2, phsm)

cmd = 'mv ' +f4+  ' ' +f1
os.system(cmd)
cmd = 'mdx.py '+f1+'.orig' + ' '+f3+' '+f1+' &'
print(cmd)   
os.system(cmd)


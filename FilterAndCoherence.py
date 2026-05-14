#!/usr/bin/env python3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Copyright 2012, by the California Institute of Technology. ALL RIGHTS RESERVED.
# United States Government Sponsorship acknowledged. Any commercial use must be
# negotiated with the Office of Technology Transfer at the California Institute of
# Technology.  This software is subject to U.S. export control laws and regulations
# and has been classified as EAR99.  By accepting this software, the user agrees to
# comply with all applicable U.S. export laws and regulations.  User has the
# responsibility to obtain export licenses, or other export authority as may be
# required before exporting such information to foreign countries or providing
# access to foreign persons.
#
# Author: Brett George
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Modified to apply the layover mask before filtering
#Francisco Delgado, IPGP, 2019/02/12

import logging
import isce
import isceobj
import argparse
import os
logger = logging.getLogger('isce.tops.runFilter')

def runFilter(infile, outfile, filterStrength):
    from mroipac.filter.Filter import Filter
    logger.info("Applying power-spectral filter")

    # Initialize the flattened interferogram
    topoflatIntFilename = infile
    intImage = isceobj.createIntImage()
    intImage.load( infile + '.xml')
    intImage.setAccessMode('read')
    intImage.createImage()

    # Create the filtered interferogram
    filtImage = isceobj.createIntImage()
    filtImage.setFilename(outfile)
    filtImage.setWidth(intImage.getWidth())
    filtImage.setAccessMode('write')
    filtImage.createImage()

    objFilter = Filter()
    objFilter.wireInputPort(name='interferogram',object=intImage)
    objFilter.wireOutputPort(name='filtered interferogram',object=filtImage)
    objFilter.goldsteinWerner(alpha=filterStrength)

    intImage.finalizeImage()
    filtImage.finalizeImage()
   

####---------------
###Francisco Delgado, IPGP, 2019/02/12, 2026/03/15
def mask_layover(igram_in, mask, igram_out):
    #######imageMath.py -e='a*(b<2)' --a=topoflatIntFilename --b=Igrams/shadowlayovermask.15alks_14rlks.rdr -t CFLOAT -o=igram_out
    #cmd = "imageMath.py -e='a*(b<3)' --a={0} --b={1} -o={2} -t CFLOAT".format(igram_in,mask,igram_out)
    cmd = "imageMath.py -e='a*(b==0)' --a={0} --b={1} -o={2} -t CFLOAT".format(igram_in,mask,igram_out)
    print(cmd)
    os.system(cmd)
def mask_water(igram_in, mask, igram_out):
    cmd = "imageMath.py -e='a*(b!=0)' --a={0} --b={1} -o={2} -t CFLOAT".format(igram_in,mask,igram_out)
    print(cmd)
    os.system(cmd)    
####---------------

def estCoherence(outfile, corfile):
    from mroipac.icu.Icu import Icu

    #Create phase sigma correlation file here
    filtImage = isceobj.createIntImage()
    filtImage.load( outfile + '.xml')
    filtImage.setAccessMode('read')
    filtImage.createImage()

    phsigImage = isceobj.createImage()
    phsigImage.dataType='FLOAT'
    phsigImage.bands = 1
    phsigImage.setWidth(filtImage.getWidth())
    phsigImage.setFilename(corfile)
    phsigImage.setAccessMode('write')
    phsigImage.createImage()

    
    icuObj = Icu(name='sentinel_filter_icu')
    icuObj.configure()
    icuObj.unwrappingFlag = False
    icuObj.useAmplitudeFlag = False

    icuObj.icu(intImage = filtImage,  phsigImage=phsigImage)
    phsigImage.renderHdr()

    filtImage.finalizeImage()
    phsigImage.finalizeImage()


def createParser():
    '''
    Create command line parser.
    '''

    parser = argparse.ArgumentParser(description='Filter interferogram and generated coherence layer.')
    parser.add_argument('-i','--input', type=str, required=True, help='Input interferogram',
            dest='infile')
    parser.add_argument('-f','--filt', type=str, default=None, help='Ouput filtered interferogram',
            dest='filtfile')
    parser.add_argument('-c', '--coh', type=str, default='phsig.cor', help='Coherence file',
            dest='cohfile')
    parser.add_argument('-s', '--strength', type=float, default=0.5, help='Filter strength',
            dest='filterstrength')

    return parser

def cmdLineParse(iargs=None):
    parser = createParser()
    return parser.parse_args(args=iargs)


def main(iargs=None):
    inps = cmdLineParse(iargs)

    if inps.filtfile is None:
        inps.filtfile = 'filt_' + inps.infile

####---------------
###Francisco Delgado, IPGP, 2019/02/12, 2026/03/15
    outfile= inps.infile[0:-4] + '.mask.int';
    infilebak0= inps.infile[0:-4] + '.orig.int';
    infilebak1= inps.infile[0:-4] + '.orig.int.vrt';
    infilebak2= inps.infile[0:-4] + '.orig.int.xml';
    print('')
    print('')
    import glob
    ######dirs=glob.glob('configs/config_igram*');
    ######dirs=dirs[0]
    ######dirs=open(dirs, 'r').readlines()
    ######alks = dirs[8]
    ######rlks = dirs[9]
    ######alks = alks[7:-1]
    ######rlks = rlks[7:-1]
    ######print(rlks+' range looks, '+alks+' azimuth looks')
    ######cmd = 'looks.py -i merged/geom_master/shadowMask.rdr -a %s -r %s' %(alks,rlks)
    ######print(cmd)
    ######os.system(cmd)
    ######maskLayover = 'merged/geom_master/shadowMask.%salks_%srlks.rdr' %(alks,rlks)
    maskLayover = 'geom_reference/shadowMask.rdr'
    print(maskLayover)
    mask_layover(inps.infile, maskLayover, outfile)
    print('')
    print('')
    cmd = 'cp %s %s' %(inps.infile,infilebak0)
    os.system(cmd)
    print(cmd)
    cmd = 'cp %s %s' %(inps.infile + '.vrt',infilebak1)
    os.system(cmd)
    print(cmd)
    cmd = 'cp %s %s' %(inps.infile + '.xml',infilebak2)
    os.system(cmd)
    print(cmd)
    cmd = 'mv %s %s' %(outfile,inps.infile)
    os.system(cmd)
    print(cmd)
    print('')
    print('')
    waterMask = 'geom_reference/waterMask.rdr'
    outfile2= inps.infile[0:-4] + '.mskw.int';
    mask_water(inps.infile, waterMask, outfile2)
    cmd = 'mv %s %s' %(outfile2,inps.infile)
    os.system(cmd)
####---------------

    runFilter(inps.infile, inps.filtfile, inps.filterstrength)

    estCoherence(inps.filtfile, inps.cohfile)

if __name__ == '__main__':
    
    main()
#!/bin/csh

set dir1 = ${1}
set dir2 = ${2}
set dirn = quad_${dir1}_${dir2}
set dem = ~/insarproc/dem/copernicus/cop_dem_glo30m_wgs84_icefields.dem
set rlks = 8
set alks = 2
set bbox = "-49.88 -47.9 -74.71 -72.17"

echo $bbox

mkdir ${dirn}
#imageMath.py -e='a*conj(b)' --a=merged/interferograms/${dir1}/fine.int --b=merged/interferograms/${dir2}/fine.int -t CFLOAT -o merged/interferograms/${dirn}/quad.int
#imageMath.py -e='a*(b==0)' --a=merged/interferograms/${dirn}/quad.int --b=merged/geom_reference/shadowMask.rdr -o=merged/interferograms/${dirn}/quad_mask.int -t CFLOAT
#FilterAndCoherence.py -i merged/interferograms/${dirn}/quad_mask.int -s 0.3 -f merged/interferograms/${dirn}/filt_quad.int -c merged/interferograms/${dirn}/quad_phsig.cor
##unwrap.py -i merged/interferograms/${dirn}/filt_quad.int -u merged/interferograms/${dirn}/filt_quad -c merged/interferograms/${dirn}/quad_phsig.cor -m snaphu_mcf
#geocodeIsce.py -f merged/interferograms/${dirn}/filt_quad.unw -d ${dem} -r ${rlks} -a ${rlks} -m reference -s reference -b "${bbox}"

#echo "geocodeIsce.py -f merged/interferograms/${dirn}/filt_quad.int -d ${dem} -r ${rlks} -a ${alks} -m reference -s reference -b "${bbox}""
geocodeIsce.py -f merged/interferograms/${dirn}/filt_quad.int -d ${dem} -r ${rlks} -a ${alks} -m reference -s reference -b "${bbox}"
geocodeIsce.py -f merged/interferograms/${dir1}/filt_fine.int -d ${dem} -r ${rlks} -a ${alks} -m reference -s reference -b "${bbox}"
geocodeIsce.py -f merged/interferograms/${dir2}/filt_fine.int -d ${dem} -r ${rlks} -a ${alks} -m reference -s reference -b "${bbox}"
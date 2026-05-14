#!/bin/csh

# Extract phase fron an snaphu unwrapped interferogram with the ISCE stack proc and export it to a png
# Francisco Delgado, IPGP, 2019/06/04

mkdir pngs
rm list_ifgs
#ls 20*/*unw -d1 > list_ifgs
#ls 20*/fine.int -d1 > list_ifgs
ls 20*/2*.int -d1 > list_ifgs
foreach ifg (`cat list_ifgs`)
	echo "INTERFEROGRAM $ifg"
	set size  = `grep size -A 1 ${ifg}.xml | head -2 | tail -1 | sed "s/[^0-9]//g"`
	#mdx -s $size $ifg -c8pha -rmg -CW -RMG-Hgt -wrap 6.28 -P; convert out.ppm -resize 50% ${ifg}.png
	mdx -s $size $ifg -c8pha -wrap 6.28 -cmap CMY -P;         convert out.ppm -resize 50% ${ifg}.png
	set fname = `echo $ifg | awk '{print "ifg_"substr($0,1,17)}'`;
	mv ${ifg}*png  ${fname}.png
	#mdx -s $size $ifg -CW -unw -r4 -rhdr 4640 -cmap cmy -wrap 20 ${ifg}.conncomp -s $size -ch1 -i1 -P; convert out.ppm ${ifg}_rewrapped.png
	mv ${fname}*png pngs/.
	##cd ..
end 

#!/bin/tcsh -f
set proj = "-JM8"
set psfile = "aniakchak_alos2_sm3_2022_2023"
clear all

gmt gmtset PS_PAGE_ORIENTATION LANDSCAPE
gmt gmtset FORMAT_GEO_MAP DF
gmt gmtset COLOR_MODEL RGB
gmt gmtset FONT_ANNOT_PRIMARY=8p
gmt gmtset COLOR_NAN 128/128/128
gmt gmtset COLOR_NAN 255/255/255
gmt gmtset COLOR_BACKGROUND 128/128/128
gmt gmtset COLOR_BACKGROUND 180/180/180
gmt gmtset COLOR_FOREGROUND 165/42/42 
gmt gmtset MAP_FRAME_TYPE plain 



set fname = "filt_220711-230612_8rlks_16alks_msk.unw.geo"
set tmp  = `gmt grdinfo -C ${fname}.grd`
set xmin = `echo $tmp | awk '{print $2}'`
set xmax = `echo $tmp | awk '{print $3}'`
set ymin = `echo $tmp | awk '{print $4}'`
set ymax = `echo $tmp | awk '{print $5}'`
set dx = `echo $tmp | awk '{print $8}'`
set dy = `echo $tmp | awk '{print $9}'`
set bounds = "-R$xmin/$xmax/$ymin/$ymax"
#set bounds = "-R-72.4/-72/-40.65/-40.35"
###set bounds = '-R-72.52500139/-71.85833499/-40.81666354/-40.23333044'
set ticks = "25"

gmt makecpt -Cseis -I -T-5/75/.1 -D211/211/211 -N > test.cpt
gmt grdgradient dem.grd -A45 -Gintens.grd -Ne1/0.6 -fg -V

set misc = "-V -K"
gmt gmtset MAP_TICK_LENGTH_PRIMARY -8p
#gmt grdimage stripmap/filt_topophase.grd -B0.1SWne -nn $proj $bounds $misc -X2 -Y5 -Iintens.grd -Ctest.cpt > $psfile.ps
gmt grdimage ${fname}.grd -B0.3SWne -nn $proj $bounds $misc -X3 -Y5 -Iintens.grd -Ctest.cpt > $psfile.ps
gmt gmtset MAP_TICK_LENGTH_PRIMARY 8p

set misc = "-V -K -O"
#echo "-72.1749 -40.4833" | gmt psxy $proj $bounds -Ss0.15 -G100/100/100 -W1,black $misc >> $psfile.ps
#echo "-72.2625 -40.4877" | gmt psxy $proj $bounds -Ss0.15 -G100/100/100 -W1,black $misc >> $psfile.ps
#echo "-72.2175 -40.4619" | gmt psxy $proj $bounds -Ss0.15 -G100/100/100 -W1,black $misc >> $psfile.ps
#echo "-72.2056 -40.5198" | gmt psxy $proj $bounds -Ss0.15 -G100/100/100 -W1,black $misc >> $psfile.ps

gmt psxy -J -R -L -G255 -W0p,black $misc << EOF >> $psfile.ps
-72.40 -40.65 
-72.33 -40.65
-72.33 -40.56
-72.40 -40.56
-72.40 -40.65 
EOF
gmt psscale -Ctest.cpt $misc -D0.2/1/1.5/0.25 -B${ticks}::/:cm: >> $psfile.ps
gmt psxy -J -R -L -G255 -W0p,black $misc << EOF >> $psfile.ps
-72.065 -40.64 
-72.010 -40.64
-72.010 -40.60
-72.065 -40.60
-72.065 -40.64 
EOF
echo "-157.8 56.4 101 0.75" | gmt psxy $proj $bounds -Sv0.15i+ea -G0/0/0 -W1,black $misc >> $psfile.ps
echo "-157.8 56.4 011 0.75" | gmt psxy $proj $bounds -Sv0.15i+ea -G170/170/170 -W1,black $misc >> $psfile.ps
gmt pstext -J -R $misc -F+f+a+j -G255/255/255 << END >> $psfile.ps
-159 57.05 12,1,black 0 BL ALOS-2 SM3 2022/07/11 - 2023/06/12
-72.04 -40.62 7,1,black 0 BL 37
END
echo "-72.146405 -40.523152 "| gmt psxy $proj $bounds -St0.2  -W1,black $misc >> $psfile.ps #2011 eruptive vent


set misc = "-V -O"
gmt psbasemap -L-158/56.42/56.42/10/10+l  $bounds $proj $misc >> $psfile.ps


gmt psconvert -A -Tf $psfile.ps
gmt psconvert -A -Tj $psfile.ps
rm intens.grd *cpt $psfile.ps gmt.history gmt.conf topo_crop.grd RHAT.pw alos2_170107_160625_nonan* 
open $psfile.pdf

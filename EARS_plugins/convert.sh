#!/bin/bash

# Assuming that we are starting with a file in lat, lon, value 
# grid to 3 arc minutes using surface (set .gmtdefaults to 6 decimal points)
# For Tahiry
awk '{print $1,$2,$3}' $1 > tmp
gmt surface tmp -: -G$2.grd -I0.00436332312998582 -R-0.523599/1.22173/0.907571/2.44346
gmt grd2xyz $2.grd > $2.tmp
# grd2xyz will put in lon lat value format 
awk '{print $2,$1,$3}' $2.tmp | sort -k1,1nr -k2,2 > $2

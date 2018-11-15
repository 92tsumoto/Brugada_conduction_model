#!/bin/sh

cat ndata_final.out | awk '{if(NR%3==2) print}' | cat -n | awk '{print $1,$7}' > tmp
awk -f cv.awk tmp | awk '{print 100, 0.1/(0.001*($2-$1))}' > 100-0.dat


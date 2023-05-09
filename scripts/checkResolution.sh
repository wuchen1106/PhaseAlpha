#!/bin/bash

#for i in anaout/*.hit.root
for i in anaout/run01802-01881.hit.root anaout/run01600-01699.hit.root
do
    runname=`echo $i | sed 's/.*\/\(.*\)\.hit.root/\1/'`
    for sizeName in width20 width30 width40 width50 width60 width70 width80 width90 width100 fwhm area
    do
        ./BinaryFiles/bin/checkResolution $i $runname $sizeName
    done
done

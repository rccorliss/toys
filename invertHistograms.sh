#!/bin/bash

#usage:  ./invertHistograms.sh [searchstring] [outputdir]

#check arguments:
if [ $# -lt 2 ]
   then echo "usage: ./invertHistograms.sh [searchstring] [outputdir]"
fi

#echo $# minus 1 is `expr $# - 1`
nFiles=`expr $# - 1`
echo found $nFiles files.
outputdir=${!#}
echo outputdir is $outputdir


for i in $(seq 1 $nFiles)
do
    #echo processing file $i
    inputname=${!i}
    basefilename=`basename $inputname`
    outputname=`echo "${basefilename/.distortion_map.hist.root/.correction_map.hist.root}"`
    outputname=${outputdir}/${outputname}
    #echo $inputname becomes $outputname
    root -l invertHistograms.C\(\"$inputname\",\"$outputname\"\)
done

exit

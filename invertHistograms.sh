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
    outputname=`echo "${basefilename/.distortion_map.hist.root/.invert.distortion_map.hist.root}"`
    outputname=${outputdir}/${outputname}
    #echo $inputname becomes $outputname
    echo root -l invertHistograms.C\(\"$inputname\",\"hugoJan2022/static_only_inverted_10-new.root\",\"hugoJan2022/static_only_closure_test.distortion_map.hist.root\",1\)
done

exit

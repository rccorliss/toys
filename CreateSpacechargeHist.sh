#!/bin/bash

#at home:
#dir="/osx_sphenix/"
#basename="G4Hits_sHijing_9-11fm_"

#at rcas:
dir="/sphenix/sim/sim01/sphnxpro/sHijing/2019-07-28/fm_9-11/"
basename="G4Hits_sHijing_9-11fm_"
maxend=3000
for ((r=0;r<=maxend;r=r+100)); do
    rp=$(($r+100))
    fname=$basename`printf "%05d\n" $r`_`printf "%05d\n" $rp`.root
    root -b -q CreateSpacechargeHist.C\(\"$dir\",\"$fname\",$r\)
done

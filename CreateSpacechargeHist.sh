#!/bin/bash

#at home: dir="/osx_sphenix/"
#basename="G4Hits_sHijing_9-11fm_"

#at rcas:
dir="/sphenix/sim/sim01/production/sHijing/new/fm_0-12/"
basename="G4Hits_sHijing_0-12fm_"
maxend=200
for ((r=0;r<=maxend;r=r+100)); do
    rp=$(($r+100))
    fname=$basename`printf "%05d\n" $r`_`printf "%05d\n" $rp`.root
    root -b -q CreateSpacechargeHist.C\(\"$dir\",\"$fname\",$r\)
done

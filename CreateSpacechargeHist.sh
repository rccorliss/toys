#!/bin/bash

#at home:
#dir="/osx_sphenix/"
#basename="G4Hits_sHijing_9-11fm_"

#at rcas:
#dir="/sphenix/sim/sim01/sphnxpro/sHijing/2019-07-28/fm_9-11/"
#basename="G4Hits_sHijing_9-11fm_"
#note that these two sims have different naming conventions for the TPC hits!  Must change the .C code to match.
#new sim:
dir="/sphenix/sim/sim01/sphnxpro/Geant4-10.05.p01/fm_0-12/FTFP_BERT_HP/"
basename="G4Hits_sHijing_0-12fm_"

maxend=1450
for ((r=0;r<=maxend;r=r+50)); do
    rp=$(($r+50))
    fname=$basename`printf "%05d\n" $r`_`printf "%05d\n" $rp`.root
    root -b -q CreateSpacechargeHist.C\(\"$dir\",\"$fname\",$r\)
done

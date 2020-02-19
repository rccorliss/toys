#!/bin/bash

#at home:
#dir="/osx_sphenix/"
#basename="G4Hits_sHijing_9-11fm_"

##########
#
# This macro runs 'CreateSpacechargeHist.C', which reads Hijing events from the
# specified directory, putting one event at every bunch crossing, and repeating
# that placement with a tilesize defined by maxend.  Bunch crossings occur at
# freq in kHz.  Ions drift at a speed hard-coded into the root macro.
#
# A more detailed version would pull a poisson distribution # of events each xing.
#
##########




#at rcas:
#dir="/sphenix/sim/sim01/sphnxpro/sHijing/2019-07-28/fm_9-11/"
#basename="G4Hits_sHijing_9-11fm_"
#note that the old sim has different naming conventions for the TPC hits (delimited by '.' instead of '#')!  Must change the .C code to match.
#new sim:
dir="/sphenix/sim/sim01/sphnxpro/Geant4-10.05.p01/fm_0-12/FTFP_BERT_HP/"
basename="G4Hits_sHijing_0-12fm_"

maxend=1000
eveperfile=50
#freq=200
for freq in 200 100 20 10
  do
      for ((r=0;r<=maxend-eveperfile;r=r+eveperfile)); do
	  rp=$(($r+$eveperfile))
	  fname=$basename`printf "%05d\n" $r`_`printf "%05d\n" $rp`.root
	  root -b -q CreateSpacechargeHist.C\(\"$dir\",\"$fname\",$r,$maxend,$freq,0\)
      done
  done

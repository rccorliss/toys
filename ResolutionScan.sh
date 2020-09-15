#!/bin/bash

chargemapname=evgeny_sept/Summary_bX1508071_0_10_events.root
np0=10
nr0=10
nz0=22

for ((offset=0;offset<=14;offset=offset+1)); do
    np=$(($np0+$offset))
    nr=$(($nr0+$offset))
    nz=$(($nz0+$offset))
    echo running:  root -b -q quick_distortion.C\\\($nr,$np,$nz,\\\"$chargemapname\\\",\\\"res_scan/\\\"\\\)
    root -b -q quick_distortion.C\($nr,$np,$nz,\"$chargemapname\",\"res_scan/\"\)
  done

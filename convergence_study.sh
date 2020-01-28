#!/bin/bash

maxred=20
for ((r=1;r<=maxred;r++)); do
    root -b -q digital_current_macro_alice.C\($r\)
    mv last_macro_output.ttree.root pre-hybrid_fixed_reduction_${r}.ttree.root
done

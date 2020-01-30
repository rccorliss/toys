#!/bin/bash

maxred=40
for ((r=1;r<=maxred;r++)); do
    root -b -q digital_current_macro_alice.C\($r\)
    mv last_macro_output.ttree.root analytic_fixed_reduction_1e8scale_${r}.ttree.root
done

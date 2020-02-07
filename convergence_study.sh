#!/bin/bash

maxred=99
for ((r=1;r<=maxred;r++)); do
    root -b -q digital_current_macro_alice.C\(-$r\)
    mv last_macro_output.ttree.root output/analytic_fixed_fraction_f${r}.ttree.root
done

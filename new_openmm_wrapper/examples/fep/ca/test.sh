#!/bin/bash

for f in sp_vdwl.yaml sp_coul.yaml sp.yaml;do
    runOpenMM_lite $f |grep "Total Potential Energy"
done

runOpenMM_lite fep_vdwl.yaml | grep "Energy for λ"
runOpenMM_lite fep_coul.yaml | grep "Energy for λ"

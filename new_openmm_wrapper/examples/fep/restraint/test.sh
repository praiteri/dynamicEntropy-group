#!/bin/bash

for f in sp0.0.yaml sp0.5.yaml sp1.0.yaml;do
    runOpenMM_lite $f |grep "Total Potential Energy"
done

runOpenMM_lite fep.yaml | grep "Energy for λ"

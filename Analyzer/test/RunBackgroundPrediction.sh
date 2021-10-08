#!/bin/bash

#cd <Path to>/CMSSW_10_6_20/src/SUSYBSMAnalysis # uncomment this line for batch job

# Compile executable
echo ">>> Compile BackgroundPrediction executable ..."
COMPILER=$(root-config --cxx)
FLAGS=$(root-config --cflags --libs)
time $COMPILER -g -O3 -Wall -Wextra -Wpedantic -o RunBackgroundPrediction Analyzer/bin/BackgroundPrediction.cpp $FLAGS

if [ "$#" -ne 0 ]; then
	echo ""
	./RunBackgroundPrediction --inputFiles $1
fi

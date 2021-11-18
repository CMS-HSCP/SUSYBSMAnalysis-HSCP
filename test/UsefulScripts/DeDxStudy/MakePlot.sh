#!/bin/bash

if [ -z "$2" ]; then
   arg2="EMPTY"
else
   arg2="$2"
fi

root -l -b << EOF
  TString makeshared(gSystem->GetMakeSharedLib());
  TString dummy = makeshared.ReplaceAll("-Wshadow ", "-D__USE_XOPEN2K8 ");
  gSystem->SetMakeSharedLib(makeshared);
  gSystem->Load("libFWCoreFWLite");
  FWLiteEnabler::enable();
  gSystem->Load("libDataFormatsFWLite.so");
  gSystem->Load("libDataFormatsHepMCCandidate.so");
  gSystem->Load("libDataFormatsCommon.so");
  gSystem->Load("libDataFormatsTrackerRecHit2D.so");
  gSystem->Load("libAnalysisDataFormatsSUSYBSMObjects.so");
  gInterpreter->SetClassAutoparsing(false);
  .x MakePlot.C+();
EOF

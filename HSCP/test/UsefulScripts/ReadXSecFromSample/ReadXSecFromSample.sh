#!/bin/bash
root -l -b << EOF
  TString makeshared(gSystem->GetMakeSharedLib());
  makeshared.ReplaceAll("-W ", "-Wno-deprecated-declarations -Wno-deprecated -Wno-unused-local-typedefs -Wno-attributes ");
  makeshared.ReplaceAll("-Woverloaded-virtual ", " ");
  makeshared.ReplaceAll("-Wshadow ", " -std=c++0x -D__USE_XOPEN2K8 ");
  cout << "Compilling with the following arguments: " << makeshared << endl;
  gSystem->SetMakeSharedLib(makeshared);
  gSystem->Load("libFWCoreFWLite");
  FWLiteEnabler::enable();
  gSystem->Load("libDataFormatsFWLite.so");
  gSystem->Load("libDataFormatsHepMCCandidate.so");
  gSystem->Load("libDataFormatsCommon.so");
  gSystem->Load("libDataFormatsTrackerRecHit2D.so");
  gSystem->Load("libAnalysisDataFormatsSUSYBSMObjects.so");
  gInterpreter->SetClassAutoparsing(false);
  .x ReadXSecFromSample.C+
EOF

#!/bin/bash
root -l -b << EOF
  TString makeshared(gSystem->GetMakeSharedLib());
  TString dummy = makeshared.ReplaceAll("-W ", "-Wno-deprecated-declarations -Wno-deprecated -D__USE_XOPEN2K8 ");
  TString dummy = makeshared.ReplaceAll("-Wshadow ", " -std=c++0x ");
  gSystem->SetMakeSharedLib(makeshared);
  gSystem->Load("libFWCoreFWLite");
  FWLiteEnabler::enable();
  gSystem->Load("libDataFormatsFWLite.so");
  gSystem->Load("libAnalysisDataFormatsSUSYBSMObjects.so");
  gSystem->Load("libDataFormatsVertexReco.so");
  gSystem->Load("libDataFormatsCommon.so");
  gSystem->Load("libDataFormatsHepMCCandidate.so");
  gSystem->Load("libPhysicsToolsUtilities.so");
  gInterpreter->SetClassAutoparsing(false);
  .x GetEfficiencyMaps.C+("COMPILE");
  //.x MakePlot.C+
EOF

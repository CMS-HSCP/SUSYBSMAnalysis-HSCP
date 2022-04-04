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
  gSystem->SetIncludePath( "-I$ROOFITSYS/include" );
  gInterpreter->SetClassAutoparsing(false);
  .x ModelIndependent_Acceptance.C+("pictures/PPStau494.txt","root://eoscms//eos/cms//store/cmst3/user/querten/12_08_30_HSCP_EDMFiles/GMStau_8TeV_M494.root");
  //.x ModelIndependent_Acceptance.C+("Analyze");
  //.x MakePlot.C+
EOF

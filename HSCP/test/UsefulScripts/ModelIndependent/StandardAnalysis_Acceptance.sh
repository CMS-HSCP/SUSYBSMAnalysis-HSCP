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
  .x StandardAnalysis_Acceptance.C++("COMPILE");
  .x StandardAnalysis_Acceptance.C+ ("PPStau_8TeV_M100", 0, 10);
  .x StandardAnalysis_Acceptance.C+ ("PPStau_8TeV_M126", 0, 20);
  .x StandardAnalysis_Acceptance.C+ ("PPStau_8TeV_M156", 0, 50);
  .x StandardAnalysis_Acceptance.C+ ("PPStau_8TeV_M200", 0, 90);
  .x StandardAnalysis_Acceptance.C+ ("PPStau_8TeV_M247", 0,130);
  .x StandardAnalysis_Acceptance.C+ ("PPStau_8TeV_M308", 0,180);
  .x StandardAnalysis_Acceptance.C+ ("PPStau_8TeV_M370", 0,230);
  .x StandardAnalysis_Acceptance.C+ ("PPStau_8TeV_M432", 0,280);
  .x StandardAnalysis_Acceptance.C+ ("PPStau_8TeV_M494", 0,320);
  .x StandardAnalysis_Acceptance.C+ ("PPStau_8TeV_M557", 0,370);
  //.x MakePlot.C+
EOF

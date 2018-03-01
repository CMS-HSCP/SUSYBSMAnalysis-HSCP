#!/bin/bash
root -l -b << EOF
  TString makeshared(gSystem->GetMakeSharedLib());
  TString dummy = makeshared.ReplaceAll("-W ", "-Wno-shadow -D__USE_XOPEN2K8 -Wno-unused-local-typedefs ");
  gSystem->SetMakeSharedLib(makeshared);
  gSystem->Load("libFWCoreFWLite");
  FWLiteEnabler::enable();
  gSystem->Load("libDataFormatsFWLite.so");
  gSystem->Load("libDataFormatsCommon.so");
  gInterpreter->SetClassAutoparsing(false);
  .x analyzeEDM.C++("SkimEff_Gluino300_pT20.0_nH3_dEdx2.8")
EOF

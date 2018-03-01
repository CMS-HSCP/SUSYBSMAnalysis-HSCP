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


  .x TriggerStudy.C+("PPStauM1599_13TeV16", "/storage/data/cms/store/user/jozobec/HSCP2016/PPStau_13TeV_M1599.root")
EOF

# FIRST BATCH
#  .x TriggerStudy.C+("Run2016_273158WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_273158.root")
#  .x TriggerStudy.C+("Run2016_273302WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_273302.root")
#  .x TriggerStudy.C+("Run2016_273402WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_273402.root")
#  .x TriggerStudy.C+("Run2016_273403WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_273403.root")
#  .x TriggerStudy.C+("Run2016_273404WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_273404.root")
#  .x TriggerStudy.C+("Run2016_273405WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_273405.root")
#  .x TriggerStudy.C+("Run2016_273406WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_273406.root")
#  .x TriggerStudy.C+("Run2016_273408WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_273408.root")
#  .x TriggerStudy.C+("Run2016_273409WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_273409.root")
#  .x TriggerStudy.C+("Run2016_273410WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_273410.root")
#  .x TriggerStudy.C+("Run2016_273411WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_273411.root")
#  .x TriggerStudy.C+("Run2016_273425WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_273425.root")
#  .x TriggerStudy.C+("Run2016_273446WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_273446.root")
#  .x TriggerStudy.C+("Run2016_273447WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_273447.root")
#  .x TriggerStudy.C+("Run2016_273448WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_273448.root")
#  .x TriggerStudy.C+("Run2016_274159WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274159.root")
#  .x TriggerStudy.C+("Run2016_274160WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274160.root")
#  .x TriggerStudy.C+("Run2016_274161WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274161.root")
#  .x TriggerStudy.C+("Run2016_274172WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274172.root")
#  .x TriggerStudy.C+("Run2016_274198WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274198.root")
#  .x TriggerStudy.C+("Run2016_274199WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274199.root")
#  .x TriggerStudy.C+("Run2016_274200WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274200.root")
#  .x TriggerStudy.C+("Run2016_274240WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274240.root")
#
# SECOND BATCH
#  .x TriggerStudy.C+("Run2016_274240WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274244.root")
#  .x TriggerStudy.C+("Run2016_274240WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274250.root")
#  .x TriggerStudy.C+("Run2016_274240WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274251.root")
#  .x TriggerStudy.C+("Run2016_274240WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274283.root")
#  .x TriggerStudy.C+("Run2016_274240WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274284.root")
#  .x TriggerStudy.C+("Run2016_274240WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274286.root")
#  .x TriggerStudy.C+("Run2016_274240WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274314.root")
#  .x TriggerStudy.C+("Run2016_274240WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274315.root")
#  .x TriggerStudy.C+("Run2016_274240WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274316.root")
#  .x TriggerStudy.C+("Run2016_274240WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274317.root")
#  .x TriggerStudy.C+("Run2016_274240WiTS", "/storage/data/cms/store/user/jozobec/HSCP2016/Run2016_274319.root")
# 
# 
# 
# 

#!/bin/bash
executable=`echo $0 | sed 's/.sh/.C/'` #assume the .sh and .C file have the same name
arguments=''
if [ $# -ge 1 ]; then arguments=$arguments"(\"`dirname $0`\"" ;fi
for var in "$@"
do
    arguments=$arguments",\"$var\"";
done
if [ $# -ge 1 ]; then arguments=$arguments");" ;fi

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
  gSystem->Load("libAnalysisDataFormatsSUSYBSMObjects.so");
  gSystem->Load("libDataFormatsVertexReco.so");
  gSystem->Load("libDataFormatsCommon.so");
  gSystem->Load("libDataFormatsHepMCCandidate.so");
  gSystem->Load("libPhysicsToolsUtilities.so");
  gInterpreter->SetClassAutoparsing(false);
  .x $executable+$arguments
EOF


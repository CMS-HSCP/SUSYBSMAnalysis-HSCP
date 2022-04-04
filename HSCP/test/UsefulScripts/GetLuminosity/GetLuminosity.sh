#!/bin/bash

outJSON=FinalBatch.json

export X509_USER_PROXY=~/x509_user_proxy/x509_proxy; voms-proxy-init --noregen;
export REMOTESTORAGEPATH=/store/user/jozobec/HSCP2016/

#root -l -b << EOF
#  TString makeshared(gSystem->GetMakeSharedLib());
#  TString dummy = makeshared.ReplaceAll("-W ", "-D__USE_XOPEN2K8 ");
#  gSystem->SetMakeSharedLib(makeshared);
#  gSystem->Load("libFWCoreFWLite");
#  FWLiteEnabler::enable();
#  gSystem->Load("libDataFormatsFWLite.so");
#  gSystem->Load("libDataFormatsCommon.so");
#  gInterpreter->SetClassAutoparsing(false);
#  .x GetLuminosity.C+
#EOF
export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.0.3/bin:$PATH
pip install --install-option="--prefix=$HOME/.local" brilws &> /dev/null #will be installed only the first time
brilcalc lumi -b "STABLE BEAMS" --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i $outJSON -u /pb -o LUMI_TABLE
cat LUMI_TABLE

pileupCalc.py -i $outJSON --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69000 --maxPileupBin 100 --numPileupBins 100  DataPileupHistogram.root
pileupCalc.py -i $outJSON --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 72475 --maxPileupBin 100 --numPileupBins 100  DataPileupUPHistogram.root
pileupCalc.py -i $outJSON --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 65550 --maxPileupBin 100 --numPileupBins 100  DataPileupDNHistogram.root

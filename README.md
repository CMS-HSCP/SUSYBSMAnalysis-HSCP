# Heavy Stable Charged Particle

## Setup working area

```bash
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_6_20
cd CMSSW_10_6_20/src/
cmsenv
```

For the following step you should have a ssh key associated to your GitLab account

```bash
git clone -b dev git@github.com:enibigir/SUSYBSMAnalysis-HSCP.git SUSYBSMAnalysis 
```
<!--
# Before compile, hide BigNTuplizer 
pushd SUSYBSMAnalysis/HSCP/plugins
mv BigNtuplizer.cc BigNtuplizer.cc.bkp
popd
-->

To compile the code, run
```bash
cd SUSYBSMAnalysis
scram b -j8
```

## Run the code

Get the scripts first:
```bash
cp HSCP/test/HSCParticleProducer_cfg.py .
cp HSCP/test/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt .
cp HSCP/test/crabConfig_Data.py .
```

Get proxy:
```bash
voms-proxy-init --voms cms -valid 192:00
```

A quick test:
```bash
#file exists?
dasgoclient -query="site file=/store/data/Run2017B/MET/AOD/09Aug2019_UL2017_rsb-v1/00000/AA1FC1E6-1E88-204D-B867-4637AEAC4BEA.root"
#if file exists on DISK (not on Tape), then do the following
cmsRun HSCParticleProducer_cfg.py LUMITOPROCESS=Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt inputFiles=root://cms-xrd-global.cern.ch//store/data/Run2017B/MET/AOD/09Aug2019_UL2017_rsb-v1/00000/AA1FC1E6-1E88-204D-B867-4637AEAC4BEA.root
```

Then run crab:
```bash
source /cvmfs/cms.cern.ch/common/crab-setup.sh
python crabConfig_Data.py
```


## EDAnalyzer on top of EDM files (created during the previous step)

```bash
cp Analyzer/test/HSCParticleAnalyzer_cfg.py .
cmsRun HSCParticleAnalyzer_cfg.py inputFiles=file:HSCP.root maxEvents=100
```

## Production of EDM files (on-fly) and run of the EDAnalyzer

```bash
cp  Analyzer/test/HSCParticleProducerAnalyzer_cfg.py .
#just like previously...
cmsRun HSCParticleProducerAnalyzer_cfg.py LUMITOPROCESS=Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt inputFiles=root://cms-xrd-global.cern.ch//store/data/Run2017B/MET/AOD/09Aug2019_UL2017_rsb-v1/00000/AA1FC1E6-1E88-204D-B867-4637AEAC4BEA.root
```
Then run crab:
```bash
cp Analyzer/test/crabConfigProdAnalyzer.py
python crabConfigProdAnalyzer_cfg.py
```

## Check of the EDAnalyzer (comparison with the old workflow)
Use the script:
```bash
Analyzer/test/compareRootFiles.py
```
This script takes two root files (to set in the file) and compares their histograms with a Kolmogorov test. Any difference is saved in:
```bash
Analyzer/test/differences.txt
```

## Background prediction

```bash
cp Analyzer/test/RunBackgroundPrediction.sh .
```
List your root files in a single text file, e.g `input.txt`

### Run locally

First compile: 
```bash
source RunBackgroundPrediction.sh
```

Then run
```bash
./RunBackgroundPrediction -h # for help
```

### Run on HTCondor

Uncomment and change the commented line in `RunBackgroundPrediction.sh`

Get submit file:
```bash
cp Analyzer/test/batch.sub .
```

Run:
```bash
condor_submit batch.sub
```

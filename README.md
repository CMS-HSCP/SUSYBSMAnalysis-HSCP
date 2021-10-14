# Heavy Stable Charged Particle

## Table of Contents

1.  [Setup working area](#setup-working-area)
1.  [Run the code](#run-the-code)
    * [Step 0](#step-0)
    * [Step 1](#step-1)
    * [Step 2](#step-2)
1.  [Pileup reweighting](#pileup-reweighting)

## Setup working area

```bash
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_6_20
cd CMSSW_10_6_20/src/
cmsenv
```

For the following step you should have a ssh key associated to your GitHub account

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

|Steps:  |                                                             |
|:---    |:------                                                      |
|Step 0  |  Produce edm files from AOD                                 |
|Step 1  |  Produce histograms and trees from step 0                   |
|Step 2  |  Estimate Background using Step 1 histograms                |
|Step 3  |  Make plots                                                 | 
|Step 4  |  Compute Limits                                             | 

### Step 0

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


### Step 1

|Analysis Type:  | |
|:---    |:------  |
|Type 0  |  Tk only |
|Type 1  |  Tk+Muon |
|Type 2  |  Tk+TOF  |
|Type 3  |  TOF only | 
|Type 4  |  Q>1 | 
|Type 5  |  Q>1 | 

#### EDAnalyzer on top of EDM files (created during the previous step)

```bash
cp Analyzer/test/HSCParticleAnalyzer_cfg.py .
cmsRun HSCParticleAnalyzer_cfg.py inputFiles=file:HSCP.root maxEvents=100
```

#### Production of EDM files (on-fly) and run of the EDAnalyzer

```bash
cp  Analyzer/test/HSCParticleProducerAnalyzer_cfg.py .
#just like previously...
cmsRun HSCParticleProducerAnalyzer_cfg.py LUMITOPROCESS=Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt inputFiles=root://cms-xrd-global.cern.ch//store/data/Run2017B/MET/AOD/09Aug2019_UL2017_rsb-v1/00000/AA1FC1E6-1E88-204D-B867-4637AEAC4BEA.root
```
Then run crab:
```bash
cp Analyzer/test/crabConfigProdAnalyzer_Data.py .
python crabConfigProdAnalyzer_Data.py
```

#### Check of the EDAnalyzer (comparison with the old workflow)
Use the script:
```bash
Analyzer/test/compareRootFiles.py
```
This script takes two root files (to set in the file) and compares their histograms with a Kolmogorov test. Any difference is saved in:
```bash
Analyzer/test/differences.txt
```

### Step 2

#### Background prediction

<!--
```bash
cp Analyzer/test/RunBackgroundPrediction.sh .
```
-->
List your root files in a single text file, e.g `input.txt`

#### Run locally

```bash
BackgroundPrediction -h # for help
BackgroundPrediction -f input.txt
```

#### Run on HTCondor

<!--
Uncomment and change the commented line in `RunBackgroundPrediction.sh`

Get submit file:
```bash
cp Analyzer/test/batch.sub .
```

Run:
```bash
condor_submit batch.sub
```
-->

## Pileup reweighting

# Heavy Stable Charged Particle

## Table of Contents

1.  [Setup working area](#setup-working-area)
1.  [Run the code](#run-the-code)
    * [Step 0](#step-0)
    * [Step 1](#step-1)
    * [Step 2](#step-2)
1.  [Compute Lumi](#compute-lumi)
<!--1.  [Pileup reweighting](#pileup-reweighting)-->

## Setup working area

```bash
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_6_27
cd CMSSW_10_6_27/src/
cmsenv
```

For the following step you should have a ssh key associated to your GitHub account.
For more information, see [connecting-to-github-with-ssh-key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent).

```bash

git clone https://github.com/CMS-HSCP/SUSYBSMAnalysis-HSCP.git SUSYBSMAnalysis 

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
|Step 0  |  Produce EDM-tuples from AOD                                 |
|Step 1  |  Produce bare ROOT-tuples (histograms and trees) from step 0                   |
|Step 2  |  Estimate Background using histograms from Step 1  |
|Step 3  |  Make plots                                                 | 
|Step 4  |  Compute Limits                                             | 

### Step 0

**Get the main script first**
```bash
cp HSCP/test/HSCParticleProducer_cfg.py .
```
Have a look at `HSCParticleProducer_cfg.py` to see all available options.

**Get proxy**
```bash
voms-proxy-init --voms cms -valid 192:00
```

**Run locally**

Before running locally:
```bash
#file exists? Note that site of type "TAPE" has no user access
dasgoclient -query="site file=/store/data/Run2017B/MET/AOD/09Aug2019_UL2017_rsb-v1/00000/AA1FC1E6-1E88-204D-B867-4637AEAC4BEA.root"
```

How to run
```bash
cmsRun HSCParticleProducer_cfg.py LUMITOPROCESS=HSCP/test/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt inputFiles=root://cms-xrd-global.cern.ch//store/data/Run2017B/MET/AOD/09Aug2019_UL2017_rsb-v1/00000/AA1FC1E6-1E88-204D-B867-4637AEAC4BEA.root
```
You can also use `inputFiles_load=input.txt`, where `input.txt` contains a list of files.

**Run on Grid using crab**

Setup CRAB environment
```bash
source /cvmfs/cms.cern.ch/common/crab-setup.sh  
```

Get CRAB configuration file to run on Data: 
```bash
cp HSCP/test/submitToCrab.py .
python submitToCrab.py -h #for help
```

**Important** Replace `config.Site.storageSite = 'T2_FR_IPHC'` with the site where you have permission to write.
You can check permission by running the following command: `crab checkwrite --site=<site-name>`. 

You can now run on Dataset (or input files)
```bash
python submitToCrab.py --dataset <Dataset> --name <request-name> --sample <isData> --lumiToProcess <JSON file>
```
The following directory will be created: `crab_projects/crab_<request-name>`.
- To get status: `crab status -d crab_projects/crab_<request-name>`
- To resubmit (killed and failed jobs): `crab resubmit -d crab_projects/crab_<request-name>`
- To retrieve the output: `crab getoutput -d crab_projects/crab_<request-name> [--jobids id1,id2]`
- To get report (processed lumi json): `crab report -d crab_projects/crab_<request-name>`
- How to get corresponding integrated lumi, see section [Compute Lumi](#compute-lumi).

<!-- -->
**Important** Running crab with `--dryrun` option provides 
1. jobs splitting
1. estimates of runtime and memory consumption. To show only splitting results, add the `--skip-estimates` option.
You can then run `crab proceed -d crab_projects/crab_<request-name>`. 
<!-- -->


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

Main script: `cp Analyzer/test/HSCParticleAnalyzer_cfg.py .`

```bash
cmsRun HSCParticleAnalyzer_cfg.py inputFiles=file:HSCP.root maxEvents=100
```
Or, if there are many files: 
```bash
ls HSCP*.root|sed 's/^/file:/'>list.txt
cmsRun HSCParticleAnalyzer_cfg.py inputFiles_load=list.txt
```

#### Production of EDM files (on-fly) and run of the EDAnalyzer

Main scripts: 
```bash
cp Analyzer/test/HSCParticleProducerAnalyzer_cfg.py .
cp Analyzer/test/crabConfigProdAnalyzer_Data.py .
```

```bash
#just like previously...
cmsRun HSCParticleProducerAnalyzer_cfg.py LUMITOPROCESS=HSCP/test/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt inputFiles=root://cms-xrd-global.cern.ch//store/data/Run2017B/MET/AOD/09Aug2019_UL2017_rsb-v1/00000/AA1FC1E6-1E88-204D-B867-4637AEAC4BEA.root
```
To use crab:
```bash
python crabConfigProdAnalyzer_Data.py
```

#### Check of the EDAnalyzer (comparison with the old workflow)
Copy the script:
```bash
cp Analyzer/test/compareRootFiles.py .
```
For more details, execute the following `python compareRootFiles.py -h`.
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

## Compute Lumi

From a new terminal

1. Setup the CMS environment 
1. Install brilws (in your home directory) as follows
   ```bash
   export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH
   pip install --user brilws
   ```
   Some tests
   ```bash
   brilcalc --version # to check the installation
   brilcalc lumi --help # for help
   ```
1. Compute lumi
   ```bash
   brilcalc lumi -c web -i Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt --output-style csv -u /pb > result.txt
   ```

<!--## Pileup reweighting

[https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData](https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData)-->

# Support, Suggestions ?
> For any support/suggestions, mail to emery.nibigiraSPAMNOT@cern.ch

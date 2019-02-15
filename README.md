# Heavy Stable Charged Particle

## Setup working area

```bash
cmsrel CMSSW_8_0_30
cd CMSSW_9_4_3/src/
cmsenv

git cms-init
# for the following step you should have a GitHub's ssh key
git clone git@github.com:CMS-HSCP/SUSYBSMAnalysis-HSCP.git SUSYBSMAnalysis/HSCP 

scram b -j8
```

-------------------------------

## Run the code

Steps:

1. Make a first skim from the AOD to some lighter AOD with only the collections needed:

	Directory : `/test/MakeEDMtuples`
	Code: 
	- HSCParticleProducerSingleMu2017_cfg.py and HSCParticleProducerMET2017_cfg.py:  configuration files to run on the SingleMuon and MET datasets, respectively
	- crab_cfg.py crab configuration file

2. test/UsefulScripts/DeDxStudy. Go to the [README](./test/UsefulScripts/DeDxStudy/README.md)

3. test/UsefulScripts/MuonTimingStudy

4. ...





# Heavy Stable Charged Particle

## Setup working area

Instruction for Combine installation taken from [Combine page](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/#slc6cc7-release-cmssw_8_1_x).

Be sure to use SL6.

```bash
cmsrel CMSSW_8_0_30
cd CMSSW_8_0_30/src/
cmsenv

git cms-init
# for the following step you should have a GitHub's ssh key
git clone git@github.com:CMS-HSCP/SUSYBSMAnalysis-HSCP.git SUSYBSMAnalysis/HSCP 

# Pick the 2016 branch
pushd SUSYBSMAnalysis/HSCP/
git checkout Run2_2016
popd

scram b -j8

#Installing Combine

git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit

cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v7.0.13

scramv1 b clean; scramv1 b

```
-------------------------------

## Run the code

### Step 1
```
cd SUSYBSMAnalysis/HSCP/test/AnalysisCode/
python Launch.py 1
```





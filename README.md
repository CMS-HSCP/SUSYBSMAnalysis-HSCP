# Heavy Stable Charged Particle

## Setup working area

```bash
cmsrel CMSSW_8_0_30
cd CMSSW_8_0_30/src/
cmsenv

git cms-init
# for the following step you should have a GitHub's ssh key
git clone git@github.com:CMS-HSCP/SUSYBSMAnalysis-HSCP.git SUSYBSMAnalysis/HSCP 

# Pick the 2016 branch
pushd SUSYBSMAnalysis/HSCP/
git checkout -b Run2_2016
popd

scram b -j8
```

-------------------------------

## Run the code

### Step 1
```
cd SUSYBSMAnalysis/HSCP/test/AnalysisCode/
python Launch.py 1
```





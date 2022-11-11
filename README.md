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
cmsrel CMSSW_10_6_30
cd CMSSW_10_6_30/src/
cmsenv
```

For the following step you should have a ssh key associated to your GitHub account.
For more information, see [connecting-to-github-with-ssh-key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent).

```bash
git clone -b master git@github.com:CMS-HSCP/SUSYBSMAnalysis-HSCP.git SUSYBSMAnalysis 
```

To compile the code, run
```bash
cd SUSYBSMAnalysis
scram b -j
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

Don't forget to copy needed files:
```bash
cp HSCP/data/CorrFact*Pix*.txt .
cp HSCP/data/template*.root .
cp HSCP/data/MuonTimeOffset.txt .
```

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

<!--## Version numbering
// v19p0
// - change double to float
// - create fillDescription
// - intro ptErrOverPt vs ptErrOverPt2
// - change the order of preselection cuts
// - N-1 plots
// - Add two more cutflow histos, change boundary for ptErrOverPt2
// - Fix logic for new cutflow, fix the  change boundary for ptErrOverPt2
// - Make cuts into an array
// - Fix logic with not used variales
// - Change the cut flow order
// - Add Ih vs globalIas_ plot in preselection, change boundary for dxy/dz plots
// - Change dxy/dz cut default
// - Add plots for MiniIsol, MET, mT
// - Change MiniIsol definition, and plot range, move it to preselection
// - Change EoP to 0.8, then to 2.0 (essentially no cut)
// - Change to allTrackMCMatch
// - 18p3: PF matching to gentracks, change the binning of MiniIso histo
// - 18p4: fix for cutflowProbQfirst index, get rid of EoP cut
// - 18p5 change to new templates
// - 18p5: remove TK iso
// - 18p8: Add postPreselection plots
// - 19p0: One try with TOF
// - 19p1: Change mass binning, remove massT cut
// - 19p3: Simplify probQ cut, change mini-iso def
// - 19p4: Change mini-iso binning
// - 19p5: use charged iso in cutflow, dont cut away out of bound probs, only in preselection
// - 19p6: intro CutFlowEta and VsGenID
// - 19p7: intro NumEvents and HSCPCandidateType, for comparrison, put back EoP cut and TkIso cut (will remove in 19p8)
// - 19p8: - Cut on PF iso electrons, no cut on EoP and TkIso - Fixed N1_ plots, renamed BS_ to BefPreS_
// - 19p9: - Futher gen printouts, change back mass histo binning
// - 19p10: - Move sibling ID and angle to histos
// - 19p14: - Angles from the mother, other gen level plots
// - 19p15: - probQvsProbXY for possibly merged clusters, Change MiniIso to all, probQ vs globalIas_ correlation
// - 19p16: - add status check for gen particles, shift layer to make plots prettier
// - 19p17: - Add 2D genPT vs recoPT plot
// - 19p18: - Add 2D genPT vs recoPT plot as PostPreS and rename to BefPreS
// - 19p19: - Cut on probXY > 0.01, add the check on special cases in pixel CPE
// - 19p20: - Cut on probXY > 0.0, and cut on isPhoton
// - 19p21: - Cut on probXY > 0.01, for real this time
// - 19p22: - Cut on probXY > 0.0, loose NOPH>1
// - 19p23: - Add GenNumSibling plots, change the default IDs to 9999
// - 20p0: - Change EoP to use PF energy
// - 20p1: - Add check if secondaries are coming from pixel NI
// - 20p2: - Add RecoPFHT and RecoPFNumJets plots, add CutFlowPfType
// - 20p3: - Change the logic of CutFlowPfType and CutFlowEta plots,
//         - add BefPreS_GenPtVsGenMinPt, and BefPreS_GenPtVsdRMinGen
//         - change the logic, that the if the closest gen in not status=1 then it's not the match
// - 20p4: - Fix20p3, move the status check out of the OR
// - 20p5: - Add ErrorHisto, TriggerType, possible fix pfType plots by interoducing the ForIdx version
// - 20p6: - Further fix for pfType?
// - 20p7: - Add PostPreS_EoPVsPfType plot, cleanup gen print-outs, move them after the preS
// - 20p8: - Add not special in CPE and !pf_isPhoton to cutflow, Extended numJetPf to 30 jets
// - 20p9: - Fix for num of mothers, not cut on special in CPE, cut on EoP < 0.3, shift the integers with 0.5 for nicer plots
// - 21p0: - Cut on ProbXY > 0.001
//v22.1 Dylan
// - 21p1 add Regions used to validate the background estimate method
// - 21p2 - Fix bug in the miniIso definition
// - 21p3: - Cut if the minDr for them is > 0.1, change to no MET triggers
// - 21p9: - Change variable names
// - 22p0: - Exclude NumHits preselection cut, change pixel hits to 2, add lepton pt to miniIso
// - 22p1: - Minor technical changes
// - 22p2: - Change probQ to no use L1 when cutting on it
// - 22p3: - Fix N1 plots, that were buggy because of 22p0 (Exclude NumHits preselection cut)
// - 22p4: - Change NOM > 10, Eta < 1.2
// - 22p5: - Change Eta < 1.0
// - 22p6: - Include reverse cutflow
// - 22p7: - Include reverse cutflow, Variable vs globalIas_ plots ( I should do variable vs probQ too)
// - 22p8: - (probXYonTrack > 0.1) and a later point in the cutflow
// - 22p9: - (probXYonTrackNoLayer1 > 0.1) 
// - 23p0: - (probXYonTrackNoLayer1 > 0.01)
// - 23p4: - Add cluster-based probXY, probQ, size per layer plots
// - 23p5: - Fix the order of probs
// - 23p6: - Restore the default CutFlow from Dylan's test cutflow after Dylan version v25
// - 23p7: - Make the probs vs layers for data and signal too, (probXYonTrackNoLayer1 > 0.1
// - 23p9: - Move printouts for Morris' study to the preselection
// - 24p0: - CluSpecInCPEVsPixelLayer add all clusters, add pthat histo, gen enviroment ID plots
// - 24p1: - Change high globalIas_ to be globalIas_ > 0.6
// - 24p2: - Add nearest jet distance
// - 23p3: - Same as 24p2
// - 24p4: - Zoom in the dR jet plot
// - 24p5: - Fix definition for dRMinJet
// - 24p6: - Tighten MiniIso cut
// - 24p7: - NOMoNOH plot, MiniIso plot boundaries, add globalMinTrackProb variables, reverse cutflow code change
// - 24p8: - Tighten MiniIso cut to 0.02
// - 24p9: - Add distance to Calo jets
// - 25p1: - Add BefPreS_dRVsPtPfJet
// - 25p2: - Tighten GlobalMinFOVH to 0.9
// - 25p3: - CutFlowProbQ plot, match pt of gen candidate, tighten dRMinGen to 0.01
// - 25p4: - No cut on pt_err/pt
// - 25p5: - Add dRMinJet vs globalIas_ plots, loosen the cut on probXY
// - 25p6: - Cut on dRMinJet
// - 25p7: - Restrict track level pixel probs by their cluster level info
// - 25p8: - ProbQ with <.75 probs, dRVsPtPfJet with 20 GeV jets
// - 25p9: - ProbQ with <.9 probs
// - 26p0: - Cleaner cutflow
// - 26p1: - Restrict track level pixel probs by their cluster level info (specInCPE)
// - 26p2: - ProbQ with <.8 probs, cut on MassErr
// - 26p3: - ProbQ with <.8 probs and no SpansTwoRocs, some printouts for Morris, dRMinJet jet def change 
// - 26p5: - Remove MassErr cut
// - 26p6: - Remove calo jet requirements for EM fraction, cut on dRMinCaloJet
// - 26p7: - Fix out of bound probXY, remove some unused 3D histos, temp remove the cut on dRMinCaloJet, put back probXYonTrackNoLayer1 cut
// - 26p8: - Tighten cut on probXYonTrackNoLayer1 to 0.1
// - 26p9: - Dont cut on probXYonTrackNoLayer1, change to 1D template CPE (instead of CR)
// - 27p0: - Run with new CPE templates
// - 27p1: - Add new plot to check pt diff for PF and Calo jets, go back to probQ def w specInCPE, cut on dRMinCaloJet > 0.4
// - 27p2: - dont cut on dRMinCaloJet, high stat version
// - 27p3: - cut on probXY > 0.01, high stat version
// - 27p4: - ProbXY plots when globalIas_ > 0.6
// - 27p5: - CluProbXY plots when globalIas_ > 0.6, local angle plots when probXY less/more than minCut, lowBetaGamma plots for pixels and strips
// - 27p6: - probs with  && probQ < 0.8
// - 27p7: - Change histo boundary for strips
// - 27p8: - Rewrite computedEdx(), add PostPreS_closestPfJet*Fraction plots, change PF def back to >20 GeV jets, strips lowBetaGamma plots with layers
// - 27p9: - Change charges to e/um, intro genGammaBetaVsProbXYNoL1, for bad CPE default probXY to probXY = 0.009 add dRMinPfMet plot
// - 28p0: - PfMetPhi and PfMet plots, dPhi PfMet plots, protection for gen history with vertex, for bad CPE default probXY to probXY = 0, and dont use it
//         - BefPreS_CluNormChargeVsStripLayer_higherBetaGamma plot,
// - 28p1: - NormClu vs layer plots for diff status particles, modify the phi distribution
// - 28p2: - Skip the track if mom ID = cand ID and has 91 status
// - 28p3: - Skip the track if it has 91 status in the env
// - 28p4: - Dont skip, but increase binning for charge vs layer
// - 28p5: - Dont skip, add charge vs layer after preS for 91 statuses
// - 28p6: - Clean the logs, skip if it has 91 status in the env
// - 28p7: - add PostPreS_P, dont cut on mini-iso and see status 91
// - 28p8: - add back mini-iso, fix the trigInfo_ (not a global variable anymore)
// - 28p9: - add lowPt pt plots, fix some boundaries, fix trigInfo_ logic on return
// - 29p0: - Frozen preselection as agreed on Sept 8
// - 29p1: - PtErrOverPt a la Dylan, plus N1 plots to study it
// - 29p2: - TNOPH plots show the nonL1Pix hits, cut on ptErr/pt2 before PtErrOverPt a la Dylan
// - 29p3: - Dont cut on ptErrOverPt, add genTrack based iso plots, plots w PU bins
// - 29p4: - As 29p3 but bug fixed
// - 29p5: - Add HLT matching
// - 29p6: - Add cut on PFMiniIso
// - 29p7: - Address the question about trigger effs (temp commit)
// - 29p8: - Revert 29p7 changes, add cut on genTrack based variable cone size abs isolation
// - 29p9: - Event level matching of muon to HLT muon, add RecoHSCParticleType plots
// - 30p0: - Cut on rel PF mini iso then on TK mini iso
// - 30p1: - Cut on E/p
// - 30p2: - Add back ptErr/pt a la Dylan
// - 30p3: - Add a very loose cut on tProbQ (0.7)
// - 30p4: - Fix logic for trigger matching, change filter to final filter, go back to no isolation cuts
// - 30p5: - Cut on rel PF mini iso
// - 30p6: - Dont cut on PF, go back to 30p4 but no cut on the distance of the HLT and muons
// - 30p7: - Fix to not have nonGlobal but standalone muons as a match, cut on dR < 0.15
// - 30p8: - Add mini-Iso
// - 30p9: - Add TkIso, add E/p cut
// - 40p0: - Add probQ cut
// - 40p1: - Add ptErr/pT2 cut
// - 40p2: - Tighter probQ cut
// - 40p3: - Loosen the probQCut, add plot for hasFilled, add some ptErrOverPt2 plots
// - 40p4: - Bugfix to 40p3
// - 40p5: - Plot 1-probQ, change naming of Ias in histos
// - 40p6: - Add Ias,GenID Vs Dz,Dxy postPreS plots
// - 40p7: - Control region with pt < 55
// - 40p8: - Go back to standard wf
// - 40p9: - Extend the massT axis, rename histos to F and G, add PU systs for FvsG



//  
//v23 Dylan 
// - v23 fix clust infos
// - add Ih and globalIas_ Pixel only no BPIXL1
// - new step2 bkg estimate
// v24 Dylan
// - add miniIso with muon contribution
// - add miniIso in the tuple

# Support, Suggestions ?
> For any support/suggestions, mail to Tamas.Almos.VamiSPAMNOT@cern.ch

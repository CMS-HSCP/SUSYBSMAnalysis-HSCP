import sys, os
import FWCore.ParameterSet.Config as cms

isSignal = False
isBckg = False
isData = True
isSkimmedSample = False
GTAG = '94X_dataRun2_ReReco_EOY17_v2'  #This is the Global Tag 
OUTPUTFILE = 'HSCP.root'
LUMITOPROCESS = ''

#debug input files 
#this list is overwritten by CRAB
InputFileList = cms.untracked.vstring(
   'root://eoscms//eos/cms/store/data/Run2015B/DoubleMuon/RECO/PromptReco-v1/000/251/521/00000/5E8B4110-6429-E511-A64C-02163E0127DF.root',
)

#main EDM tuple cfg that depends on the above parameters
execfile( os.path.expandvars('${CMSSW_BASE}/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/HSCParticleProducer_cfg.py') )

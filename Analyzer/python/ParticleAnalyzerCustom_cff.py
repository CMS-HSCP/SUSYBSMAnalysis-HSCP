import FWCore.ParameterSet.Config as cms
import os

PATH_TO_DATA = "{}/src/SUSYBSMAnalysis/HSCP/data".format(os.getenv('CMSSW_BASE'))

analyzer = cms.EDAnalyzer('TriggerStudyAnalyzer'
    # collections
    ,hscpCollection   = cms.InputTag("HSCParticleProducer") 
    ,hscpIsoCollection = cms.InputTag("HSCPIsolation", "R03") #New format used for data since 17-07-2015, old: ("HSCPIsolation03")
    ,l1results   = cms.InputTag("gtStage2Digis","","RECO") #not sure?
    ,triggerResults = cms.InputTag("TriggerResults","","HLTX")
    ,pfMET = cms.InputTag("pfMet", "", "RECO")
    ,pfJet = cms.InputTag("ak4PFJetsCHS", "", "RECO")
    ,CaloMET = cms.InputTag("caloMet", "", "RECO")
    ,genParticleCollection = cms.InputTag("genParticlesSkimmed")
    #HLT triggers
    ,Trigger_Mu = cms.untracked.vstring("HLT_Mu45_eta2p1","HLT_Mu50_v") #"HLT_OldMu100_","HLT_TkMu100_","HLT_TkMu50_v"
    #,Trigger_MET  = cms.untracked.vstring("HLT_PFMET170_HBHECleaned_v","HLT_PFMET170_NoiseCleaned")
    ,Trigger_MET  = cms.untracked.vstring("HLT_PFMET120_PFMHT120_IDTight_v","HLT_PFHT500_PFMET100_PFMHT100_IDTight_v","HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v","HLT_MET105_IsoTrk50_v")
    # analysis parameters
    ,TypeMode        = cms.untracked.uint32(0) # 0:Tk only, 1:Tk+Muon, 2:Tk+TOF, 3:TOF onlypwd, 4:Q>1, 5:Q<1
    ,SampleType      = cms.untracked.uint32(2) # 0:Data, 1:Background, 2:Signal, 3:Signal Systematics
    ,SampleName      = cms.untracked.string("BaseName")
    ,Period          = cms.untracked.string("2016")

)

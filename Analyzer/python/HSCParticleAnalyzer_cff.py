import FWCore.ParameterSet.Config as cms
import os

PATH_TO_DATA = "{}/src/SUSYBSMAnalysis/HSCP/data".format(os.getenv('CMSSW_BASE'))

analyzer = cms.EDAnalyzer('Analyzer'
    # collections
    ,hscpCollection   = cms.InputTag("HSCParticleProducer")
    ,hscpIsoCollection = cms.InputTag("HSCPIsolation", "R03") #New format used for data since 17-07-2015, old: ("HSCPIsolation03")
    ,dedxCollection   = cms.InputTag("dedxHitInfo")
    ,muonTimeCollection    = cms.InputTag("muons","combined")
    ,muonDtTimeCollection  = cms.InputTag("muons","dt")
    ,muonCscTimeCollection = cms.InputTag("muons","csc")
    ,muonDtSegmentCollection  = cms.InputTag("dt4DSegments")
    ,muonCscSegmentCollection = cms.InputTag("cscSegments")
    ,offlinePrimaryVerticesCollection = cms.InputTag("offlinePrimaryVertices")
    ,refittedStandAloneMuonsCollection = cms.InputTag("refittedStandAloneMuons")
    ,offlineBeamSpotCollection = cms.InputTag("offlineBeamSpot")
    ,muonSegmentCollection = cms.InputTag("MuonSegmentProducer")
    ,muonCollection = cms.InputTag("muons")
    ,triggerResults = cms.InputTag("TriggerResults","","HLT")
    #HLT triggers
    ,Trigger_MET  = cms.untracked.vstring("HLT_PFMET120_PFMHT120_IDTight_v","HLT_PFMET170_NoiseCleaned")
    ,Trigger_Mu   = cms.untracked.vstring("HLT_Mu45_eta2p1","HLT_Mu50_v") #"HLT_OldMu100_","HLT_TkMu100_","HLT_TkMu50_v"
    # analysis parameters
    ,TypeMode        = cms.untracked.uint32(0) # 0:Tk only, 1:Tk+Muon, 2:Tk+TOF, 3:TOF onlypwd, 4:Q>1, 5:Q<1
    ,SampleType      = cms.untracked.uint32(0) # 0:Data, 1:MC, >=2:Signal
    ,BaseName        = cms.untracked.string("Data")
    # skip some distribution and trees
    ,SkipSelectionPlot   = cms.untracked.bool(False)
    # histogram bounds
    ,PtHistoUpperBound   = cms.untracked.double(4000)
    ,MassHistoUpperBound = cms.untracked.double(4000)
    ,MassNBins           = cms.untracked.uint32(400)
    ,IPbound             = cms.untracked.double(1.0)
    ,PredBins            = cms.untracked.uint32(0)
    ,EtaBins             = cms.untracked.uint32(60)
    ,dEdxS_UpLim         = cms.untracked.double(1.0)
    ,dEdxM_UpLim         = cms.untracked.double(30.0)
    ,DzRegions           = cms.untracked.uint32(6)
    ,GlobalMinPt         = cms.untracked.double(55.00)
    ,GlobalMinTOF        = cms.untracked.double(1.0)
    # scale factor and K & C
    #,DeDxSF_0        = cms.untracked.double(1.00000) #=1 if data
    #,DeDxSF_1        = cms.untracked.double(1.41822)
    #,DeDxK           = cms.untracked.double(2.275)
    #,DeDxC           = cms.untracked.double(3.675)
    # calibration
    ,enableDeDxCalibration = cms.untracked.bool(False)
    ,DeDxCalibration       = cms.untracked.string("{}/Data13TeVGains_v2.root".format(PATH_TO_DATA))
    ,DeDxTemplate          = cms.untracked.string("{}/Data13TeV16_dEdxTemplate.root".format(PATH_TO_DATA))
    ,Geometry              = cms.untracked.string("{}/CMS_GeomTree.root".format(PATH_TO_DATA))
    ,TimeOffset            = cms.untracked.string("{}/MuonTimeOffset.txt".format(PATH_TO_DATA))
)
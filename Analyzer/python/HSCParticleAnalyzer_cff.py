import FWCore.ParameterSet.Config as cms
import os

PATH_TO_DATA = "{}/src/SUSYBSMAnalysis/HSCP/data".format(os.getenv('CMSSW_BASE'))

analyzer = cms.EDAnalyzer('Analyzer'
    # collections
    ,hscpCollection   = cms.InputTag("HSCParticleProducer")
    ,dedxCollection   = cms.InputTag("dedxHitInfo")
    ,muonTimeCollection    = cms.InputTag("muons","combined")
    ,muonDtTimeCollection  = cms.InputTag("muons","dt")
    ,muonCscTimeCollection = cms.InputTag("muons","csc")
    ,muonDtSegmentCollection  = cms.InputTag("dt4DSegments")
    ,muonCscSegmentCollection = cms.InputTag("cscSegments")
    ,offlinePrimaryVerticesCollection = cms.InputTag("offlinePrimaryVertices")
    ,refittedStandAloneMuonsCollection = cms.InputTag("refittedStandAloneMuons")
    ,offlineBeamSpotCollection = cms.InputTag("offlineBeamSpot")
    ,hscpIsoCollection = cms.InputTag("HSCPIsolation", "R03")
    ,muonSegmentCollection = cms.InputTag("MuonSegmentProducer")
    # Configurations
    ,Debug           = cms.untracked.bool(False)
    ,TypeMode        = cms.untracked.uint32(0) # 0:Tk only, 1:Tk+Muon, 2:Tk+TOF, 3:TOF onlypwd, 4:Q>1, 5:Q<1
    ,SampleType      = cms.untracked.uint32(0) # 0:Data, 1:MC, >=2:Signal
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
    ,DeDxSF_0        = cms.untracked.double(1.00000)
    ,DeDxSF_1        = cms.untracked.double(1.41822)
    ,DeDxK           = cms.untracked.double(2.275)
    ,DeDxC           = cms.untracked.double(3.675)
    # calibration
    ,DeDxTemplate    = cms.untracked.string("{}/Data13TeV16_dEdxTemplate.root".format(PATH_TO_DATA))
    ,DeDxCalibration = cms.untracked.string("{}/Data13TeVGains_v2.root".format(PATH_TO_DATA))
    ,Geometry        = cms.untracked.string("{}/CMS_GeomTree.root".format(PATH_TO_DATA))
    ,TimeOffset      = cms.untracked.string("{}/MuonTimeOffset.txt".format(PATH_TO_DATA))
)
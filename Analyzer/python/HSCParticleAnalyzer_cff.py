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
    ,lumiScalers    = cms.InputTag("scalersRawToDigi")
    ,refittedStandAloneMuonsCollection = cms.InputTag("refittedStandAloneMuons")
    ,offlineBeamSpotCollection = cms.InputTag("offlineBeamSpot")
    ,muonSegmentCollection = cms.InputTag("MuonSegmentProducer")
    ,muonCollection = cms.InputTag("muons")
    ,triggerResults = cms.InputTag("TriggerResults","","HLT")
    ,pfMET = cms.InputTag("pfMet", "", "RECO")
    ,pfJet = cms.InputTag("ak4PFJetsCHS", "", "RECO")
    ,CaloMET = cms.InputTag("caloMet", "", "RECO")
    ,pileupInfo = cms.InputTag("addPileupInfo")
    ,genParticleCollection = cms.InputTag("genParticlesSkimmed")
    #HLT triggers
    ,Trigger_Mu = cms.untracked.vstring("HLT_Mu50_v")
    ,Trigger_MET  = cms.untracked.vstring("HLT_PFMET120_PFMHT120_IDTight_v","HLT_PFHT500_PFMET100_PFMHT100_IDTight_v","HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v","HLT_MET105_IsoTrk50_v")
    # analysis parameters
    ,TypeMode        = cms.untracked.uint32(0) # 0:Tk only, 1:Tk+Muon, 2:Tk+TOF, 3:TOF onlypwd, 4:Q>1, 5:Q<1
    ,SampleType      = cms.untracked.uint32(0) # 0:Data, 1:Background, 2:Signal, 3:Signal Systematics
    ,SampleName      = cms.untracked.string("BaseName")
    ,Period          = cms.untracked.string("2016")
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
    ,GlobalMaxPterr      = cms.untracked.double(0.25)
    ,GlobalMinPt         = cms.untracked.double(55.00)
    ,GlobalMinTOF        = cms.untracked.double(1.0)
    ,skipPixel           = cms.untracked.bool(True)
    ,useTemplateLayer    = cms.untracked.bool(False)
    # scale factor and K & C -- values obtained with harm-2 and 15% low values drop 
    ,DeDxSF_0        = cms.untracked.double(1.00000) #=1 if data
    ,DeDxSF_1        = cms.untracked.double(1.6107*1.0448500) #SF for run > 279479
    #,DeDxSF_1        = cms.untracked.double(1.41822)
    ,DeDxK           = cms.untracked.double(2.37) 
    ,DeDxC           = cms.untracked.double(2.93)
    ,FMIPX           = cms.untracked.double(4)
    ,saveTree        = cms.untracked.uint32(0) #0: tree not saved
    ,saveGenTree     = cms.untracked.uint32(0) #0: tree not saved
    # calibration
    ,enableDeDxCalibration = cms.untracked.bool(False)
    ,DeDxCalibration       = cms.untracked.string("{}/Data13TeVGains_v2.root".format(PATH_TO_DATA))
    ,DeDxTemplate          = cms.untracked.string("{}/dEdxTemplate_harm2_SO_in_noC_CCC_MG_2017B.root".format(PATH_TO_DATA))
    ,Geometry              = cms.untracked.string("{}/CMS_GeomTree.root".format(PATH_TO_DATA))
    ,TimeOffset            = cms.untracked.string("{}/MuonTimeOffset.txt".format(PATH_TO_DATA))
    ,pixelCPE              = cms.string("PixelCPEClusterRepair")
    ,trackProbQCut         = cms.untracked.double(1.0)
    ,debugLevel            = cms.untracked.uint32(0)
)

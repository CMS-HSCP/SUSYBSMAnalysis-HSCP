import FWCore.ParameterSet.Config as cms
import os

PATH_TO_DATA = "{}/src/SUSYBSMAnalysis/HSCP/data".format(os.getenv('CMSSW_BASE'))

analyzer = cms.EDAnalyzer('Analyzer'
    # collections
    #,tracks = cms.InputTag("generalTracks")
    ,hscpCollection   = cms.InputTag("HSCParticleProducer")
    ,dedxCollection   = cms.InputTag("dedxHitInfo")
    ,muonTimeCollection    = cms.InputTag("muons","combined")
    ,muonDtTimeCollection  = cms.InputTag("muons","dt")
    ,muonCscTimeCollection = cms.InputTag("muons","csc")
    ,muonDtSegmentCollection  = cms.InputTag("dt4DSegments")
    ,muonCscSegmentCollection = cms.InputTag("cscSegments")
    # Configurations
    ,Debug           = cms.untracked.bool(False)
    ,AddTree         = cms.untracked.bool(True)
    ,TypeMode        = cms.untracked.uint32(0) # 0:Tk only, 1:Tk+Muon, 2:Tk+TOF, 3:TOF onlypwd, 4:Q>1, 5:Q<1
    ,SampleType      = cms.untracked.uint32(0) # 0:Data, 1:MC, >=2:Signal
    ,SampleTxtFile   = cms.untracked.string("Analysis_Samples.txt")
    ,DeDxSF_0        = cms.untracked.double(1.00000)
    ,DeDxSF_1        = cms.untracked.double(1.41822)
    ,DeDxK           = cms.untracked.double(2.275)
    ,DeDxC           = cms.untracked.double(3.675)
    ,DeDxTemplate    = cms.untracked.string("{}/Data13TeV16_dEdxTemplate.root".format(PATH_TO_DATA))
    ,DeDxCalibration = cms.untracked.string("{}/Data13TeVGains_v2.root".format(PATH_TO_DATA))
    ,Geometry        = cms.untracked.string("{}/CMS_GeomTree.root".format(PATH_TO_DATA))
    ,TimeOffset      = cms.untracked.string("{}/MuonTimeOffset.txt".format(PATH_TO_DATA))
    #,dedxHitInfos = cms.InputTag("dedxHitInfo")
    #,miniTracks = cms.untracked.uint32(1000)
)
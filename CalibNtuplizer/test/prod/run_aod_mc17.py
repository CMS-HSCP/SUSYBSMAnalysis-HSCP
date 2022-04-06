import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = "106X_mc2017_realistic_v9" #pdmv 5 oct 2021

process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
         # aod 2017 w+1j
         'root://cmsxrootd.fnal.gov//store/mc/RunIISummer20UL17RECO/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_mc2017_realistic_v6-v1/00000/01EDA01D-B02E-754C-8882-A1906F60202C.root'
         # aod prod malgorzata
         #'root://cmsxrootd.fnal.gov//store/user/kazana/HSCP/MC2017/UL17_mc_PU/HSCPgluino_M_1000_TuneCP5_13TeV_pythia8/UL17_mc_PU/210326_181209/0000/RECO_1.root'
    )
)

process.stage = cms.EDAnalyzer('ntuple'
     , format_file       = cms.string('AOD')
     , isdata            = cms.bool(False)
     , year              = cms.untracked.int32(2017)
     , primaryVertexColl  = cms.InputTag('offlinePrimaryVertices')  #->AOD
     , isotracks             = cms.InputTag("isolatedTracks")
     , tracks             = cms.InputTag("generalTracks")
     , collectionHSCP     = cms.InputTag("HSCParticleProducer")
     , isoHSCP0            = cms.InputTag("HSCPIsolation","R005")
     , isoHSCP1            = cms.InputTag("HSCPIsolation","R01")
     , isoHSCP2            = cms.InputTag("HSCPIsolation","R03")
     , isoHSCP3            = cms.InputTag("HSCPIsolation","R05")
     , dedx               = cms.InputTag("dedxHitInfo")
     , MiniDedx           = cms.InputTag("isolatedTracks")
     , dEdxHitInfoPrescale = cms.InputTag("dedxHitInfo","prescale")
     , muons             = cms.InputTag("muons")
     , muonTOF            = cms.InputTag("muons","combined")
     , muonTDT            = cms.InputTag("muons","dt")
     , muonTCSC            = cms.InputTag("muons","csc")
     , printOut           = cms.untracked.int32(-1)
     , GenPart            = cms.InputTag("genParticles")
     , runOnGS            = cms.bool(False)
     , stripSimLinks      = cms.InputTag("simSiStripDigis")
     , ROUList            = cms.vstring(
                                      'TrackerHitsTIBLowTof',  # hard-scatter (prompt) collection, module label g4simHits
                                      'TrackerHitsTIBHighTof',
                                      'TrackerHitsTIDLowTof',
                                      'TrackerHitsTIDHighTof',
                                      'TrackerHitsTOBLowTof',
                                      'TrackerHitsTOBHighTof',
                                      'TrackerHitsTECLowTof',
                                      'TrackerHitsTECHighTof'
                                      )
    , lumiScalerTag     = cms.InputTag("scalersRawToDigi")
    , doRecomputeMuTim  = cms.bool(False)  # no recomputation for MC !! 
    , cscSegments     = cms.InputTag("cscSegments")
    , dt4DSegments    = cms.InputTag("dt4DSegments")
    , pfCand              = cms.InputTag("particleFlow", "", "RECO")
    , pfMET = cms.InputTag("pfMet", "", "RECO")
    , pfJet = cms.InputTag("ak4PFJetsCHS", "", "RECO")
    , CaloMET = cms.InputTag("caloMet", "", "RECO")
    , pixelCPE              = cms.string("PixelCPEClusterRepair")
    , trackProbQCut         = cms.untracked.double(1.0)
)



process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('nt_mc_aod.root')
 )

process.load("SUSYBSMAnalysis.HSCP.HSCParticleProducer_cff")
process.HSCParticleProducer.filter = cms.bool( False )  # no filter for MC !!

process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.HSCPTrigger = process.hltHighLevel.clone()
process.HSCPTrigger.TriggerResultsTag = cms.InputTag( "TriggerResults", "", "HLT" )
process.HSCPTrigger.andOr = cms.bool( True ) #OR
process.HSCPTrigger.throw = cms.bool( False )
process.HSCPTrigger.HLTPaths = [ #check triggers
#          "HLT_PFMET110_PFMHT110_IDTight_v*",  # sus-19-009 (2016) on MET_MHT
          "HLT_PFMET120_PFMHT120_IDTight_v*",  # sus-19-009 (2016-2018) 
          "HLT_PFMET170_NoiseCleaned_v*",      # sus-19-009 (2016) on MET, HSCP 2016
#          "HLT_PFMET170_JetIdCleaned_v*",      # sus-19-009 (2016)
#          "HLT_PFMET170_HBHECleaned_v*",       # sus-19-009 (2016)
#          "HLT_PFMET170_NotCleaned_v*",        # sus-19-009 (2016)
          "HLT_Mu50_v*",                       # fractionaly charged on SingleMu
          "HLT_OldMu100_*",                     # fractionaly charged on SingleMu
          "HLT_TkMu100_*",                      # fractionaly charged on SingleMu
          "HLT_Mu45_eta2p1*",                  # HSCP 2016
          "HLT_TkMu50_v*",                      # April2020 studies -> email 05/05
]

#process.p = cms.Path(process.HSCPTrigger + process.HSCParticleProducerSeq +process.stage)
process.p = cms.Path(process.HSCParticleProducerSeq +process.stage)  # no filter for MC!!
#process.p = cms.Path(process.dump)




## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.MessageLogger.cerr.default.limit = 10

import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = "94X_mc2017_realistic_v17"


#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
         # aod 2017
         'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17DRPremix/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/RECOSIMstep_94X_mc2017_realistic_v10-v1/00000/785EFE4D-43F1-E711-B3E1-0CC47A74527A.root'
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

import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = "102X_dataRun2_Prompt_v14"


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
         # file dataset=/ZeroBias/Run2018D-PromptReco-v2/AOD
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/273/00000/B7E2D71D-6AEE-1043-80AC-55165A154ADB.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/270/00000/8968E56A-37EA-7B4B-ABF6-D11F39169478.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/261/00000/32693582-388F-2D46-81F1-97393ECB5238.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/245/00000/BEB5B4B4-638E-8846-980B-9E4E9C366B8A.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/240/00000/556E2AA4-5FF3-B04D-B981-ADE7E38BBAA2.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/237/00000/56CC01C7-B07D-974F-9887-BDD436ABE2B4.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/236/00000/64ECA511-9B05-3E43-BDE4-930FD0C43727.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/235/00000/E24D80F7-0ED3-2D48-A8F6-937BA27D4B86.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/233/00000/052D7C9F-6C21-1847-BEDD-215DC974FBAE.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/232/00000/D49B9739-568F-644B-AB66-00248FF3D792.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/231/00000/81E7708B-0C20-6144-97C6-CA57599CB3B8.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/230/00000/3C6C9D7C-3EE0-1942-B9E8-BD5175CDA63A.root'
    )
)

process.stage = cms.EDAnalyzer('ntuple'
     , format_file       = cms.string('AOD')
     , isdata            = cms.bool(True)
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
)



process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('/opt/sbg/cms/ui3_data1/ccollard/HSCP/calib/nt_zerobias_aod.root')
 )

process.load("SUSYBSMAnalysis.HSCP.HSCParticleProducer_cff")

process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.HSCPTrigger = process.hltHighLevel.clone()
process.HSCPTrigger.TriggerResultsTag = cms.InputTag( "TriggerResults", "", "HLT" )
process.HSCPTrigger.andOr = cms.bool( True ) #OR
process.HSCPTrigger.throw = cms.bool( False )
process.HSCPTrigger.HLTPaths = [ #check triggers
#          "HLT_PFMET110_PFMHT110_IDTight_v*",  # sus-19-009 (2016) on MET_MHT
          "HLT_PFMET120_PFMHT120_IDTight_v*",  # sus-19-009 (2016-2018) 
#          "HLT_PFMET170_NoiseCleaned_v*",      # sus-19-009 (2016) on MET
#          "HLT_PFMET170_JetIdCleaned_v*",      # sus-19-009 (2016)
#          "HLT_PFMET170_HBHECleaned_v*",       # sus-19-009 (2016)
#          "HLT_PFMET170_NotCleaned_v*",        # sus-19-009 (2016)
          "HLT_Mu50_v*",                       # fractionaly charged on SingleMu
          "HLT_OldMu100_*",                     # fractionaly charged on SingleMu
          "HLT_TkMu100_*",                      # fractionaly charged on SingleMu
]

#process.p = cms.Path(process.HSCPTrigger + process.HSCParticleProducerSeq +process.stage)
process.p = cms.Path(process.HSCParticleProducerSeq +process.stage)
#process.p = cms.Path(process.dump)




## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.MessageLogger.cerr.default.limit = 10

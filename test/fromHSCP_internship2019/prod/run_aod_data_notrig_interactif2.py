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
         # file dataset=/ZeroBias/Run2018D-PromptReco-v2/AOD run=325022 site=T2_FI_HIP
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/FE93FE51-9A1A-8841-99BD-DCE50CC9262E.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/FDF621C8-3E80-B343-AB9B-05335677F2DC.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/FCAE89F8-8EA8-CE46-B5E6-2BE698FEF130.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/F86F78A9-455B-DE47-B7BE-FB42AD924C07.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/F3FD91CA-484C-9247-8B6B-B00F401A6901.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/F365CE5B-EABD-B344-BD8B-CAE1C8F51A24.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/F22CDE4A-CDC8-F649-A3EA-9F219578DD30.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/F2074423-E1C6-624D-BBF0-430274230C74.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/F175AAB8-5CFE-984F-9AE5-07C5E25377FC.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/F11E6507-AE46-6449-95DC-092831CF6F30.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/F0E3697D-04DB-074F-8A72-EF62324313D4.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/EFCFFA23-40FE-A645-A2CC-91ABA4A82F9C.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/EEEA3E2A-8F03-864A-91E0-6A03A25E1FED.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/EEE6F900-D973-3844-A4BB-564CE757AC4F.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/EE71973C-F897-8F49-BFBF-B6485C695DDE.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/EAAB40F7-09F6-7442-B866-C1D9D24B8A6E.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/E9408414-4C18-9D44-AABA-8DB642C21015.root',
         'root://cmsxrootd.fnal.gov//store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/325/022/00000/E8DE4F64-3FD8-A947-B574-D13AD40EF68E.root'
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
     fileName = cms.string('/opt/sbg/cms/ui3_data1/ccollard/HSCP/calib/nt_zerobias_aod_2.root')
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

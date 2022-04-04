import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
         # gluino 2016
         'root://cmsxrootd.fnal.gov//store/mc/Run3Summer19DR/HSCPgluino_M_1800_TuneCP5_14TeV_pythia8/AODSIM/HSCP_SIM_106X_mcRun3_2024_realistic_v4-v2/260000/E3F00FAB-DE78-7B43-AECB-A31B6BEE8BBB.root',
         'root://cmsxrootd.fnal.gov//store/mc/Run3Summer19DR/HSCPgluino_M_1800_TuneCP5_14TeV_pythia8/AODSIM/HSCP_SIM_106X_mcRun3_2024_realistic_v4-v2/260000/7C749E19-F085-E74C-8A75-21925E5FF646.root',
         'root://cmsxrootd.fnal.gov//store/mc/Run3Summer19DR/HSCPgluino_M_1800_TuneCP5_14TeV_pythia8/AODSIM/HSCP_SIM_106X_mcRun3_2024_realistic_v4-v2/260000/74997560-6AA0-E843-AB54-9F220CE24BC1.root',
         'root://cmsxrootd.fnal.gov//store/mc/Run3Summer19DR/HSCPgluino_M_1800_TuneCP5_14TeV_pythia8/AODSIM/HSCP_SIM_106X_mcRun3_2024_realistic_v4-v2/260000/21912A8F-B46A-DA43-AD7C-756C606A6565.root'
    )
)

process.stage = cms.EDAnalyzer('ntuple'
     , tracks             = cms.InputTag("generalTracks")
     , dedx               = cms.InputTag("dedxHitInfo")
     , dEdxHitInfoPrescale = cms.InputTag("dedxHitInfo","prescale")
     , muons             = cms.InputTag("muons")
     , printOut           = cms.untracked.int32(-1)
     , GenPart            = cms.InputTag("prunedGenParticles")
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
)

process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('/opt/sbg/cms/ui3_data1/ccollard/HSCP/gluino1800_run3/test_gluino_2024.root')
 )

process.p = cms.Path(process.stage)
#process.p = cms.Path(process.dump)


## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.MessageLogger.cerr.default.limit = 10

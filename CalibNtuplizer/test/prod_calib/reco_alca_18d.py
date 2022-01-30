import FWCore.ParameterSet.Config as cms

process = cms.Process("DEDXUNCSKIM")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')

process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')


process.options   = cms.untracked.PSet(
      wantSummary = cms.untracked.bool(True),
      SkipEvent = cms.untracked.vstring('ProductNotFound'),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(200)
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:root://xrootd-cms.infn.it//store/data/Run2018D/ZeroBias/ALCARECO/SiStripCalMinBias-PromptReco-v2/000/325/022/00000/01238A92-7BE9-424E-860B-F0178FB44485.root'
    ),
   inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")
)

#process.GlobalTag.globaltag = GTAG
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "102X_dataRun2_v11", '')


process.load('Alignment.CommonAlignmentProducer.AlignmentTrackSelector_cfi')
process.tracksForDeDx = process.AlignmentTrackSelector.clone(
    src = 'generalTracks',
    filter = True,
    applyBasicCuts = True,
    ptMin = 0.3,
    nHitMin = 6,
    chi2nMax = 10.,
)
process.tracksForDeDx.src = cms.InputTag('ALCARECOSiStripCalMinBias') #for SiStripCalMinBias ALCARECO format

process.load('RecoVertex.BeamSpotProducer.BeamSpot_cff')
from RecoTracker.TrackProducer.TrackRefitter_cfi import *
process.RefitterForDeDx = TrackRefitter.clone(
      src = cms.InputTag("tracksForDeDx"),
      NavigationSchool = cms.string("")                                   
)

from RecoTracker.DeDx.dedxEstimators_cff import *
process.dedxHitInfo = dedxHitInfo.clone()
process.dedxHitInfo.tracks=cms.InputTag("RefitterForDeDx")
process.dedxHitInfo.trajectoryTrackAssociation = cms.InputTag("RefitterForDeDx")
process.dedxHitInfo.minTrackPt = cms.double(0.0)

#make the pool output
process.Out = cms.OutputModule("PoolOutputModule",
     outputCommands = cms.untracked.vstring(
         "drop *",
         "keep EventAux_*_*_*",
         "keep LumiSummary_*_*_*",
         "keep *_RefitterForDeDx_*_DEDXUNCSKIM",
         "keep *_dedxHitInfo_*_DEDXUNCSKIM",
    ),
    fileName = cms.untracked.string("dEdxSkim.root"),
    SelectEvents = cms.untracked.PSet(
       SelectEvents = cms.vstring('*')
    ),
)


process.stage = cms.EDAnalyzer('calib_ntuple'
#     , primaryVertexColl  = cms.InputTag('offlinePrimaryVertices')  #->AOD
     , tracks             = cms.InputTag("RefitterForDeDx")
     , dedx               = cms.InputTag("dedxHitInfo")
     , dEdxHitInfoPrescale = cms.InputTag("dedxHitInfo","prescale")
     , printOut           = cms.untracked.int32(-1)
)


# Output service 
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('calib_ntuple.root')
)

#schedule the sequence
process.p = cms.Path(process.tracksForDeDx * process.offlineBeamSpot * process.RefitterForDeDx * process.dedxHitInfo * process.stage)
#process.p = cms.Path(process.tracksForDeDx * process.offlineBeamSpot * process.RefitterForDeDx * process.dedxHitInfo )
process.endPath1 = cms.EndPath(process.Out)
process.schedule = cms.Schedule(process.p, process.endPath1)


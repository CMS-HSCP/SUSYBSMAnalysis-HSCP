import sys, os
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')

# defaults
options.outputFile = 'Histos.root'
options.maxEvents = -1 # -1 means all events

options.register('SAMPLE', 'isSignal',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Sample Type. Use: isSignal or isBckg or isData"
)
options.register('isSkimmedSample', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "is sample Skimmed? True or False"
)
options.register('LUMITOPROCESS', 'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Lumi to process"
)
options.parseArguments()


process = cms.Process("HSCPAnalysis")

#diventano var parsing
#The following parameters need to be provided
#isSignal, isBckg, isData, isSkimmedSample, GTAG, InputFileList
isSignal = True
isBckg = False
isData = False
isSkimmedSample = False

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("SUSYBSMAnalysis.Analyzer.HSCParticleAnalyzer_cff")
process.load("SUSYBSMAnalysis.HSCP.HSCParticleProducer_cff")

process.options   = cms.untracked.PSet(
      #wantSummary = cms.untracked.bool(True),
      #SkipEvent = cms.untracked.vstring('ProductNotFound'),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring("file:HSCP_Gluino_Mass1800_AOD_1.root"),
   inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = "106X_upgrade2018_realistic_v11_L1v1"


process.genParticlesSkimmed = cms.EDFilter("GenParticleSelector",
  filter = cms.bool(False),
  src = cms.InputTag("genParticles"),
  cut = cms.string('pt > 5.0'),
  stableOnly = cms.bool(True)
)

########################################################################
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")

process.trackExtraRekeyer = cms.EDProducer("TrackExtraRekeyer",
  src = cms.InputTag("generalTracks"),
  association = cms.InputTag("muonReducedTrackExtras"),
)

import RecoTracker.TrackProducer.TrackRefitter_cfi
process.myRefittedTracks = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone()
process.myRefittedTracks.src= 'trackExtraRekeyer'
process.myRefittedTracks.NavigationSchool = ''
process.myRefittedTracks.Fitter = 'FlexibleKFFittingSmoother'
process.myRefittedTracks.TrajectoryInEvent = True
process.myRefittedTracks.TTRHBuilder = "WithAngleAndTemplate"

# run the EDAnalyzer
### set your configirattion here (default: python/HSCParticleAnalyzer_cff.py)
process.analyzer.TypeMode=0
process.analyzer.SampleType=2
process.analyzer.saveTree=0 #all saved
process.analyzer.saveGenTree=0
process.analyzer.DeDxCalibration="Data13TeVGains_v2.root"
process.analyzer.DeDxTemplate="dEdxTemplate_harm2_SO_in_noC_CCC_MG_2017B.root"
process.analyzer.Geometry="CMS_GeomTree.root"
process.analyzer.TimeOffset="MuonTimeOffset.txt"
process.analyzer.trajInputLabel    = cms.untracked.InputTag('myRefittedTracks')
process.analyzer.triggerResults = cms.InputTag("TriggerResults","","RECO")
process.analyzer.probQCut = cms.untracked.double(0.1)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string(options.outputFile)
)

process.TrackRefitter_step = cms.Path(
  process.offlineBeamSpot*
  process.trackExtraRekeyer*
  process.myRefittedTracks
)

process.analyzer_step = cms.Path(
#  process.HSCParticleProducerSeq*process.analyzer
   process.genParticlesSkimmed*process.HSCParticleProducerSeq*process.analyzer
)


#schedule the sequence
process.schedule = cms.Schedule(
    process.TrackRefitter_step,
    process.analyzer_step,
)


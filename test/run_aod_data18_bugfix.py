import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
from Configuration.Eras.Modifier_pf_badHcalMitigation_cff import pf_badHcalMitigation

process = cms.Process("HSCPAnalysis",Run2_2018,pf_badHcalMitigation)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')


## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.options = cms.untracked.PSet(
#    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('MismatchedInputFiles'),
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = "106X_dataRun2_v27"

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# /SingleMuon/Run2018*-12Nov2019_UL2018*/AOD
# /SingleMuon/Run2018D-12Nov2019_UL2018-v4/AOD
#root://cms-xrd-global.cern.ch//store/data/Run2018D/SingleMuon/AOD/12Nov2019_UL2018-v4/100000/63BF07E8-4AAA-DC4D-99DE-48B52FCB0C06.root
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
         'root://cms-xrd-global.cern.ch//store/data/Run2018D/SingleMuon/AOD/12Nov2019_UL2018-v4/100000/63BF07E8-4AAA-DC4D-99DE-48B52FCB0C06.root'
    ),
        secondaryFileNames = cms.untracked.vstring(
         'root://cms-xrd-global.cern.ch//store/data/Run2018D/SingleMuon/ALCARECO/SiPixelCalSingleMuon-ForPixelALCARECO_UL2018-v1/250000/EE812551-E563-5B44-B6AC-F8190A191ABF.root',
    ),
    eventsToProcess = cms.untracked.VEventRange('324237:37111027')
)

#---------------------- Refitter -----------------------
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")

process.MeasurementTrackerEvent.pixelClusterProducer = 'ALCARECOSiPixelCalSingleMuon::RECO'
process.MeasurementTrackerEvent.stripClusterProducer = 'ALCARECOSiPixelCalSingleMuon::RECO'
process.MeasurementTrackerEvent.inactivePixelDetectorLabels = cms.VInputTag()
process.MeasurementTrackerEvent.inactiveStripDetectorLabels = cms.VInputTag()

#process.TrackRefitter.src = cms.InputTag('globalMuons::RECO')
process.TrackRefitter.src = 'ALCARECOSiPixelCalSingleMuon::RECO'
#process.TrackRefitter.src = 'generalTracks'
process.TrackRefitter.TrajectoryInEvent = True
process.TrackRefitter.NavigationSchool = ''
process.TrackRefitter.TTRHBuilder = "WithAngleAndTemplate"

# CMS Path
process.TrackRefitter_step = cms.Path(
  process.offlineBeamSpot*
  process.MeasurementTrackerEvent*
  process.TrackRefitter
)

process.schedule = cms.Schedule(
    process.TrackRefitter_step
    #,process.HSCParticleProducer_step,
    #process.Ntuple_step,
    #process.endjob_step
)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)



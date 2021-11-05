import sys, os
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')

# defaults
options.outputFile = 'Histos.root'
options.maxEvents = -1 # -1 means all events

options.register('GTAG', '106X_dataRun2_v20',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Global Tag"
)
options.register('SAMPLE', 'isData',
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
#isSignal = True
#isBckg = False
#isData = False
#isSkimmedSample = False
#GTAG = 'START72_V1::All'

## print configuration:
print('\nCMSSW version : {}'.format(os.environ['CMSSW_VERSION']))
print('Global Tag    : {}'.format(options.GTAG))
if options.SAMPLE=='isData':
   print('Lumi File     : {}'.format(options.LUMITOPROCESS))
print('Sample Type   : {}'.format(options.SAMPLE))
print('is skimmed    : {}'.format(options.isSkimmedSample))
print('Output File   : {}'.format(options.outputFile))
print('Input Files   : {}\n'.format(options.inputFiles))

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')

process.options   = cms.untracked.PSet(
      wantSummary = cms.untracked.bool(True),
      SkipEvent = cms.untracked.vstring('ProductNotFound'),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring("file:/afs/cern.ch/work/t/tvami/public/HSCP/NewFrameWork/2018_SingleMuon_AOD_100.root"),
   inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")
)
#if(options.SAMPLE=='isSignal'): process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck") 


#The duplicateCheckMode works only if we submit with Condor - not with Crab - checks process history, run number, lumi number
process.source.duplicateCheckMode = cms.untracked.string("checkAllFilesOpened")


from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.GTAG, '')

process.HSCPTuplePath = cms.Path() 
#J: should I add something like:  process.HSCPTuplePath += process.source.duplicateCheckMode ?

########################################################################
#Run the Skim sequence if necessary
if(not options.isSkimmedSample):
   process.nEventsBefSkim  = cms.EDProducer("EventCountProducer")

   process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
   process.HSCPTrigger = process.hltHighLevel.clone()
   process.HSCPTrigger.TriggerResultsTag = cms.InputTag( "TriggerResults", "", "HLT" )
   process.HSCPTrigger.andOr = cms.bool( True ) #OR
   process.HSCPTrigger.throw = cms.bool( False )
   if(options.SAMPLE=='isData'):
      process.HSCPTrigger.HLTPaths = [ #check triggers
          "HLT_PFMET120_PFMHT120_IDTight_v*",
          "HLT_Mu50_v*",
          "HLT_PFHT500_PFMET100_PFMHT100_IDTight_v*",
          "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v*",
          "HLT_MET105_IsoTrk50_v*",
      ]
   elif(options.SAMPLE=='isBckg'):
      #to be updated to Run2 Triggers, in the meanwhile keep all of them to study trigger efficiency
      process.HSCPTrigger.HLTPaths = [
          "HLT_PFMET120_PFMHT120_IDTight_v*", #check triggers  
          "HLT_Mu50_v*",
          "HLT_PFHT500_PFMET100_PFMHT100_IDTight_v*",
          "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v*",
          "HLT_MET105_IsoTrk50_v*",
      ]
   else:
      #do not apply trigger filter on signal
      process.HSCPTrigger.HLTPaths = ["*"]  
   
   process.HSCPTuplePath += process.nEventsBefSkim + process.HSCPTrigger 

########################################################################

#Run the HSCP EDM-tuple Sequence on skimmed sample
process.nEventsBefEDM   = cms.EDProducer("EventCountProducer")
process.load("SUSYBSMAnalysis.HSCP.HSCParticleProducer_cff") 
process.HSCPTuplePath += process.nEventsBefEDM + process.HSCParticleProducerSeq

########################################################################  
# Only for MC samples, save skimmed genParticles

if(options.SAMPLE=='isSignal' or options.SAMPLE=='isBckg'):
   process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
   process.genParticlesSkimmed = cms.EDFilter("GenParticleSelector",
        filter = cms.bool(False),
        src = cms.InputTag("genParticles"),
        cut = cms.string('pt > 5.0'),
        stableOnly = cms.bool(True)
   )

   process.HSCPTuplePath += process.genParticlesSkimmed

########################################################################

#make the pool output
process.Out = cms.OutputModule("PoolOutputModule",
     outputCommands = cms.untracked.vstring(
         "drop *",
         "keep EventAux_*_*_*",
         "keep LumiSummary_*_*_*",
         "keep edmMergeableCounter_*_*_*",
         "keep GenRunInfoProduct_*_*_*",
         "keep GenEventInfoProduct_generator_*_*",
         "keep *_genParticlesSkimmed_*_*",
         "keep *_genParticlePlusGeant_*_*",
         "keep *_offlinePrimaryVertices_*_*",
         "keep recoTracks_generalTracks_*_*",
         "keep recoTracks_standAloneMuons_*_*",
         "keep recoTrackExtras_standAloneMuons_*_*",
         "keep TrackingRecHitsOwned_standAloneMuons_*_*",
         "keep recoTracks_globalMuons_*_*",  
         "keep recoTrackExtras_globalMuons_*_*",
         "keep recoMuons_muons_*_*",
         "keep recoMuonTimeExtraedmValueMap_muons_*_*",
         "keep edmTriggerResults_TriggerResults_*_*",
         "keep *_ak4PFJetsCHS__*", 
         "keep recoPFMETs_pfMet__*",     
         "keep *_HSCParticleProducer_*_*",
         "keep *_HSCPIsolation*_*_*",
         "keep *_dedxHitInfo*_*_*",
         "keep triggerTriggerEvent_hltTriggerSummaryAOD_*_*",
         "keep *_offlineBeamSpot_*_*",
         "keep *_MuonSegmentProducer_*_*",
         "keep *_g4SimHits_StoppedParticles*_*",
         "keep PileupSummaryInfos_addPileupInfo_*_*",
         "keep *_dt4DSegments__*",  
         "keep *_cscSegments__*",  
         "keep *_scalersRawToDigi_*_*", 
         "keep *_caloMet_*_*",
    ),
    fileName = cms.untracked.string(options.outputFile),
    SelectEvents = cms.untracked.PSet(
       SelectEvents = cms.vstring('*')
    ),
)

if(options.SAMPLE=='isData' and len(options.LUMITOPROCESS)>0):
   import FWCore.PythonUtilities.LumiList as LumiList
   #process.source.lumisToProcess = LumiList.LumiList(filename = options.LUMITOPROCESS).getVLuminosityBlockRange()
   #process.source.lumisToProcess = LumiList.LumiList(url = https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt).getVLuminosityBlockRange()

if(options.SAMPLE=='isBckg' or options.SAMPLE=='isData'):
   process.Out.SelectEvents.SelectEvents =  cms.vstring('HSCPTuplePath')  #take just the skimmed ones
   process.Out.outputCommands.extend(["drop triggerTriggerEvent_hltTriggerSummaryAOD_*_*"])
else:
   process.Out.SelectEvents = cms.untracked.PSet()


########################################################################
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")

process.trackExtraRekeyer = cms.EDProducer("TrackExtraRekeyer",
                                         src = cms.InputTag("generalTracks"),
                                         association = cms.InputTag("muonReducedTrackExtras"),
#                                         association = cms.InputTag("slimmedMuonTrackExtras"),
                                         )

import RecoTracker.TrackProducer.TrackRefitter_cfi
process.myRefittedTracks = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone()
process.myRefittedTracks.src= 'trackExtraRekeyer'
process.myRefittedTracks.NavigationSchool = ''
process.myRefittedTracks.Fitter = 'FlexibleKFFittingSmoother'
process.myRefittedTracks.TrajectoryInEvent = True
process.myRefittedTracks.TTRHBuilder = "WithAngleAndTemplate"

# run the EDAnalyzer



process.load("SUSYBSMAnalysis.Analyzer.HSCParticleAnalyzer_cff")
### set your configirattion here (default: python/HSCParticleAnalyzer_cff.py)
#process.analyzer.SampleTxtFile=options.sampleTxtFile
process.analyzer.TypeMode=0
process.analyzer.SampleType=0
process.analyzer.saveTree=0 #all saved
process.analyzer.saveGenTree=0
process.analyzer.DeDxCalibration="Data13TeVGains_v2.root"
process.analyzer.DeDxTemplate="dEdxTemplate_harm2_SO_in_noC_CCC_MG_2017B.root"
process.analyzer.Geometry="CMS_GeomTree.root"
process.analyzer.TimeOffset="MuonTimeOffset.txt"
process.analyzer.trajInputLabel    = cms.untracked.InputTag('myRefittedTracks')

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string(options.outputFile)
                                   )
###process.analyzer.OutputFile = 'Data_2017_UL'

process.analysis = cms.Path(process.analyzer)
#process.p = cms.Path(process.analyzer*process.dump)
#process.analyzer = cms.EDAnalyzer('Analyzer')

#process.prodAnalysis = cms.Path(process.HSCPTuplePath*process.analyzer)

process.TrackRefitter_step = cms.Path(
  process.offlineBeamSpot*
#  process.MeasurementTrackerEvent*
  process.trackExtraRekeyer*
  process.myRefittedTracks
#  process.TrackRefitter
)

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.endjob_step = cms.EndPath(process.endOfProcess)

process.HSCPTuplePath += process.analyzer

########################################################################

process.tsk = cms.Task()
for mod in process.producers_().itervalues():
    process.tsk.add(mod)
for mod in process.filters_().itervalues():
    process.tsk.add(mod)

#schedule the sequence
process.endPath1 = cms.EndPath(process.Out)
process.schedule = cms.Schedule(process.TrackRefitter_step, process.HSCPTuplePath, process.endjob_step)
#process.schedule = cms.Schedule(process.prodAnalysis, process.endjob_step)
#process.schedule = cms.Schedule(process.HSCPTuplePath, process.analysis, process.endjob_step)
#process.schedule = cms.Schedule(process.HSCPTuplePath, process.endPath1)


# process task? in 9.4

import sys, os
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('python')

options.register('GTAG', '94X_dataRun2_ReReco_EOY17_v2',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Global Tag"
)
options.register('OUTPUTFILE', 'HSCP.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output File"
)

options.register('LUMITOPROCESS', '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Lumi to process"
)
options.parseArguments()


process = cms.Process("HSCPAnalysis")

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
print '------------'
print 'inputFiles'
print options.inputFiles
print '------------'
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options.inputFiles),
                            #fileNames =  options.inputFiles,
                            inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")
                            )

#The duplicateCheckMode works only if we submit with Condor - not with Crab - checks process history, run number, lumi number
process.source.duplicateCheckMode = cms.untracked.string("checkAllFilesOpened")


from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.GTAG, '')

process.HSCPTuplePath = cms.Path() 

########################################################################
#Run the Skim sequence if necessary

process.nEventsBefSkim  = cms.EDProducer("EventCountProducer")

process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.HSCPTrigger = process.hltHighLevel.clone()
process.HSCPTrigger.TriggerResultsTag = cms.InputTag( "TriggerResults", "", "HLT" )
process.HSCPTrigger.andOr = cms.bool( True ) #OR
process.HSCPTrigger.throw = cms.bool( False )
process.HSCPTrigger.HLTPaths = [ #check triggers
   "HLT_Mu50_v*",
   ]
process.HSCPTuplePath += process.nEventsBefSkim + process.HSCPTrigger 

########################################################################

#Run the HSCP EDM-tuple Sequence on skimmed sample
process.nEventsBefEDM   = cms.EDProducer("EventCountProducer")
process.load("SUSYBSMAnalysis.HSCP.HSCParticleProducer_cff") 
process.HSCPTuplePath += process.nEventsBefEDM + process.HSCParticleProducerSeq


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
    ),
    fileName = cms.untracked.string(options.OUTPUTFILE),
    SelectEvents = cms.untracked.PSet(
       SelectEvents = cms.vstring('*')
    ),
)


import FWCore.PythonUtilities.LumiList as LumiList
if options.LUMITOPROCESS: 
    process.source.lumisToProcess = LumiList.LumiList(filename = options.LUMITOPROCESS).getVLuminosityBlockRange()

process.Out.SelectEvents.SelectEvents =  cms.vstring('HSCPTuplePath')  #take just the skimmed ones
process.Out.outputCommands.extend(["drop triggerTriggerEvent_hltTriggerSummaryAOD_*_*"])

########################################################################

process.tsk = cms.Task()
for mod in process.producers_().itervalues():
    process.tsk.add(mod)
for mod in process.filters_().itervalues():
    process.tsk.add(mod)

#schedule the sequence
process.endPath1 = cms.EndPath(process.Out)
process.schedule = cms.Schedule(process.HSCPTuplePath, process.endPath1)


# process task? in 9.4

import sys, os
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

options = VarParsing('analysis')

# defaults
options.outputFile = 'Histos.root'
# -1 means all events
options.maxEvents = -1

#options.register('GTAG', '106X_upgrade2018_realistic_v11_L1v1',
options.register('GTAG', '106X_upgrade2018_realistic_v11BasedCandidateTmp_2022_08_09_01_32_34',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Global Tag"
)
options.register('SAMPLE', 'isBckg',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Sample Type. Use: isSignal or isBckg or isData"
)
options.register('YEAR', '2018',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Year. Use: 2017 or 2018"
)
options.register('isSkimmedSample', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "is sample Skimmed? True or False"
)
#options.register('LUMITOPROCESS', '',
options.register('LUMITOPROCESS', 'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Lumi to process"
)
options.parseArguments()

process = cms.Process("HSCPAnalysis")
#process = cms.Process("HSCPAnalysis",Run2_2018)

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
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.Services_cff')

process.options   = cms.untracked.PSet(
      wantSummary = cms.untracked.bool(True),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source("PoolSource",
   #fileNames = cms.untracked.vstring("/store/mc/RunIISummer20UL18RECO/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/230000/064A8795-8468-3849-B543-BDD6287EE510.root"),
   fileNames = cms.untracked.vstring("/store/mc/RunIISummer20UL18RECO/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/280005/D8AB7663-12E6-6247-BF03-0F24B7D7D4C6.root "),
   inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.GTAG, '')

process.HSCPTuplePath = cms.Path() 

########################################################################
#Run the Skim sequence if necessary
if(not options.isSkimmedSample):
   process.nEventsBefSkim  = cms.EDProducer("EventCountProducer")

   process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
   process.HSCPTrigger = process.hltHighLevel.clone()
   process.HSCPTrigger.TriggerResultsTag = cms.InputTag( "TriggerResults", "", "HLT" )
   process.HSCPTrigger.andOr = cms.bool( True ) #OR
   process.HSCPTrigger.throw = cms.bool( False )  
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
#        stableOnly = cms.bool(True)
   )

   process.HSCPTuplePath += process.genParticlesSkimmed

########################################################################
# electron VID
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
electron_id_config = cms.PSet(electron_ids = cms.vstring([                   
                    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff',
                    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff',
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff', 
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
                    ]))  


                 
switchOnVIDElectronIdProducer(process,DataFormat.AOD)
for idmod in electron_id_config.electron_ids.value():
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.HSCPTuplePath += process.egmGsfElectronIDSequence

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
   process.source.lumisToProcess = LumiList.LumiList(filename = options.LUMITOPROCESS).getVLuminosityBlockRange()
   #process.source.lumisToProcess = LumiList.LumiList(url = https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt).getVLuminosityBlockRange()

if(options.SAMPLE=='isBckg' or options.SAMPLE=='isData'):
   process.Out.SelectEvents.SelectEvents =  cms.vstring('HSCPTuplePath')  #take just the skimmed ones
   process.Out.outputCommands.extend(["drop triggerTriggerEvent_hltTriggerSummaryAOD_*_*"])
else:
   process.Out.SelectEvents = cms.untracked.PSet()


########################################################################

# run the EDAnalyzer

if options.SAMPLE=='isData' :
   SampleType = 0
   if options.YEAR=='2017' :
       K = 2.30
       C = 3.17
       SF0 = 1.0
       SF1 = 1.0325
       IasTemplate = "template_2017C.root" 
   
   if options.YEAR=='2018' :
       K = 2.27
       C = 3.16
       SF0 = 1.0
       SF1 = 1.0817
       IasTemplate = "template_2017C.root" #FIXME template 2018?
    #HSCP_minPt = 55

elif options.SAMPLE=='isBckg':
   SampleType = 1
   if options.YEAR=='2017' :
       K = 2.26
       C = 3.22
       SF0 = 1.0079
       SF1 = 1.0875
       IasTemplate = "templateMC.root"
    
   if options.YEAR=='2018' :
       K = 2.27
       C = 3.22
       SF0 = 1.0047
       SF1 = 1.1429
       IasTemplate = "templateMC.root"

else :
   SampleType = 2
   if options.YEAR=='2017' :
       K = 2.26
       C = 3.22
       SF0 = 1.0079
       SF1 = 1.0875
       IasTemplate = "templateMC.root"
    
   if options.YEAR=='2018' :
       K = 2.27
       C = 3.22
       SF0 = 1.0047
       SF1 = 1.1429
       IasTemplate = "templateMC.root"

process.load("SUSYBSMAnalysis.Analyzer.HSCParticleAnalyzer_cfi")
process.HSCParticleAnalyzer.TypeMode = 0 # 0: Tracker only
process.HSCParticleAnalyzer.SampleType = SampleType 
process.HSCParticleAnalyzer.SaveTree = 6 #6 is all saved, 0 is none
process.HSCParticleAnalyzer.SaveGenTree = 0
process.HSCParticleAnalyzer.DeDxTemplate=IasTemplate
process.HSCParticleAnalyzer.TimeOffset="MuonTimeOffset.txt"
process.HSCParticleAnalyzer.Period = "2018"
process.HSCParticleAnalyzer.DebugLevel = 0 
process.HSCParticleAnalyzer.DeDxK = K
process.HSCParticleAnalyzer.DeDxC = C
process.HSCParticleAnalyzer.DeDxSF_0 = SF0
process.HSCParticleAnalyzer.DeDxSF_1 = SF1
process.HSCParticleAnalyzer.GlobalMinIh = C

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string(options.outputFile)
                                   )

process.analysis = cms.Path(process.HSCParticleAnalyzer)

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.endjob_step = cms.EndPath(process.endOfProcess)

process.HSCPTuplePath += process.HSCParticleAnalyzer

########################################################################

process.tsk = cms.Task()
for mod in process.producers_().itervalues():
    process.tsk.add(mod)
for mod in process.filters_().itervalues():
    process.tsk.add(mod)

#schedule the sequence
process.endPath1 = cms.EndPath(process.Out)
process.schedule = cms.Schedule(process.HSCPTuplePath, process.endjob_step)


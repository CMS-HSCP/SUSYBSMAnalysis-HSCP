import sys, os
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')

# defaults
options.outputFile = 'Histos.root'
options.maxEvents = 10000 # -1 means all events

# was 106X_dataRun2_v20
options.register('GTAG', '106X_dataRun2_v36',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Global Tag"
)
options.register('SAMPLE', 'isData',
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
process.load("SUSYBSMAnalysis.Analyzer.metFilters_cff")

process.options   = cms.untracked.PSet(
      wantSummary = cms.untracked.bool(True),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source("PoolSource",
   # fileNames = cms.untracked.vstring('file:/eos/uscms/store/user/lpchscp/noreplica/008C9591-D7BF-FF41-B484-AE207413C07F.root'),
   # fileNames = cms.untracked.vstring("/store/data/Run2017B/SingleMuon/AOD/15Feb2022_UL2017-v1/2820000/0014B98B-5C06-A140-82C4-38E4C1BE6366.root"),
   # fileNames = cms.untracked.vstring("/store/data/Run2017C/SingleMuon/AOD/15Feb2022_UL2017-v1/2560000/8B124DEA-526B-394E-8BB6-B67B1FF5EC4B.root"), # this one seems to be on tape now...
   # fileNames = cms.untracked.vstring("/store/data/Run2017C/SingleMuon/AOD/15Feb2022_UL2017-v1/2820000/006B1C4B-D30C-014D-A85A-4A4F46591DA3.root"), #file 1/4 on disk: didn't work
   # fileNames = cms.untracked.vstring("/store/data/Run2017C/SingleMuon/AOD/15Feb2022_UL2017-v1/2820000/04D68268-64FD-3C4A-8B0D-B16BA6095A43.root"), #file 2/4 on disk: seems to work!
   # fileNames = cms.untracked.vstring("/store/data/Run2017C/SingleMuon/AOD/15Feb2022_UL2017-v1/2820000/022CE3B9-9141-4E4C-8B54-432F7F940ED5.root"), #file 3/4 on disk:
   # fileNames = cms.untracked.vstring("/store/data/Run2017C/SingleMuon/AOD/15Feb2022_UL2017-v1/2820000/02CE93D9-96B1-EF47-8364-A9D2FD5AF7C9.root"), #file 4/4 on disk:
   # fileNames = cms.untracked.vstring("/store/data/Run2017D/SingleMuon/AOD/15Feb2022_UL2017-v1/2560000/004E92B3-AC0F-6947-A6D6-A2D2FEFB095D.root"),
   # fileNames = cms.untracked.vstring("/store/data/Run2017E/SingleMuon/AOD/15Feb2022_UL2017-v1/2560000/09081324-7FFE-AD43-892A-67B617FEF398.root"),
   # fileNames = cms.untracked.vstring("/store/data/Run2017F/SingleMuon/AOD/15Feb2022_UL2017-v1/2560000/01C1F2EC-6368-0240-B53F-E3DE2453A607.root"), # this one not working even though it's on disk
   # fileNames = cms.untracked.vstring("/store/data/Run2017A/MET/AOD/09Aug2019_UL2017_rsb-v1/120000/006BDDD8-AF7E-EE4B-ACD2-28C9B5A2972F.root"),
   # fileNames = cms.untracked.vstring("/store/data/Run2017B/MET/AOD/09Aug2019_UL2017_rsb-v1/00000/000FC5BA-7C2F-D543-8D50-85B037E34CE1.root"),
   # fileNames = cms.untracked.vstring("/store/data/Run2017C/MET/AOD/09Aug2019_UL2017_rsb-v1/10000/008C9591-D7BF-FF41-B484-AE207413C07F.root"), # now moved to tape
    fileNames = cms.untracked.vstring("/store/data/Run2017D/MET/AOD/09Aug2019_UL2017_rsb-v1/120000/005E3EF6-8DEC-4841-A1E3-3DBC6321FDDD.root"),
                            
   inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")
)

#The duplicateCheckMode works only if we submit with Condor - not with Crab - checks process history, run number, lumi number
process.source.duplicateCheckMode = cms.untracked.string("checkAllFilesOpened")

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



process.HSCPTuplePath += process.nEventsBefEDM + process.HSCParticleProducerSeq + process.metFilters

########################################################################  
# Only for MC samples, save skimmed genParticles

if(options.SAMPLE=='isSignal' or options.SAMPLE=='isBckg'):
   process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
   process.genParticlesSkimmed = cms.EDFilter("GenParticleSelector",
        filter = cms.bool(False),
        src = cms.InputTag("genParticles"),
        cut = cms.string('pt > 0.0'),
        stableOnly = cms.bool(True)
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
for mod in process.producers_().values():
    process.tsk.add(mod)
for mod in process.filters_().values():
    process.tsk.add(mod)

#schedule the sequence
process.endPath1 = cms.EndPath(process.Out)

process.schedule = cms.Schedule(process.HSCPTuplePath,  process.endjob_step)



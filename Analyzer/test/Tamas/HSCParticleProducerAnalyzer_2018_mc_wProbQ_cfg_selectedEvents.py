import sys, os
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')

# defaults
options.outputFile = 'Histos.root'
# -1 means all events
options.maxEvents = -1

options.register('GTAG', '106X_upgrade2018_realistic_v11_L1v1',
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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source("PoolSource",
   #fileNames = cms.untracked.vstring("/store/mc/RunIISummer20UL18RECO/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/230000/064A8795-8468-3849-B543-BDD6287EE510.root"),
  # fileNames = cms.untracked.vstring("/store/mc/RunIISummer20UL18RECO/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/280005/D8AB7663-12E6-6247-BF03-0F24B7D7D4C6.root "),
#   fileNames = cms.untracked.vstring("file:F8A9F740-F226-D443-A132-45CBC706B908.root"),
#   fileNames = cms.untracked.vstring("file:D49DD2CE-E848-9E42-9CB0-AE4E6C60280A.root"),
  # fileNames = cms.untracked.vstring("file:589444D4-2FA5-2D42-B952-F719D35D1EF4.root"),
  # fileNames = cms.untracked.vstring("file:3A0707EB-0FC4-E74B-9491-27C51157FB89.root"),
 #  fileNames = cms.untracked.vstring("file:B406DEC6-85A0-E342-A677-281FC853ABEE.root"),
 #  fileNames = cms.untracked.vstring("/store/mc/RunIISummer20UL18RECO/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/100000/215F439C-73BF-8D47-AE10-681B747F764F.root"),
   fileNames = cms.untracked.vstring("/store/mc/RunIISummer20UL18RECO/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/2430000/E09ACB33-2178-7346-9B8F-1B2E37A01299.root"),
   inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")
)

#process.source.eventsToProcess = cms.untracked.VEventRange('1:29057:317155112')
#process.source.eventsToProcess = cms.untracked.VEventRange('1:30965:337976289')
#process.source.eventsToProcess = cms.untracked.VEventRange('1:34335:374762778')
#process.source.eventsToProcess = cms.untracked.VEventRange('1:9798:91439758')
#process.source.eventsToProcess = cms.untracked.VEventRange('1:1707:15925987')
#process.source.eventsToProcess = cms.untracked.VEventRange('1:21841:165597238')
process.source.eventsToProcess = cms.untracked.VEventRange('1:30134:441662388')



from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.GTAG, '')

process.HSCPTuplePath = cms.Path() 

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


process.load("SUSYBSMAnalysis.Analyzer.HSCParticleAnalyzer_cff")
### set your configirattion here (default: python/HSCParticleAnalyzer_cff.py)
#process.analyzer.SampleTxtFile=options.sampleTxtFile
process.analyzer.TypeMode = 0 # 0: Tracker only
process.analyzer.SampleType = SampleType 
process.analyzer.SaveTree = 0 #6 is all saved, 0 is none
process.analyzer.SaveGenTree = 0
process.analyzer.DeDxTemplate=IasTemplate
process.analyzer.TimeOffset="MuonTimeOffset.txt"
process.analyzer.TrackProbQCut = 1.0
process.analyzer.Period = "2018"
process.analyzer.DebugLevel = 10 
process.analyzer.DeDxK = K
process.analyzer.DeDxC = C
process.analyzer.DeDxSF_0 = SF0
process.analyzer.DeDxSF_1 = SF1
process.analyzer.GlobalMinIh = C

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string(options.outputFile)
                                   )

process.analysis = cms.Path(process.analyzer)

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.endjob_step = cms.EndPath(process.endOfProcess)

process.HSCPTuplePath += process.analyzer

#schedule the sequence
process.schedule = cms.Schedule(process.HSCPTuplePath, process.endjob_step)


import os
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')

# defaults
options.outputFile = 'AnalysisTuple.root'
options.maxEvents = -1 # -1 means all events

options.register('TYPE', 0, 
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "0:Tk only, 1:Tk+Muon, 2:Tk+TOF, 3:TOF onlypwd, 4:Q>1, 5:Q<1"
)
options.register('SAMPLE', 'isData',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Sample Type. Use: isData or isBckg or isSignal or isSignalSyst" 
)
options.register('NAME', 'Data',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of the tree (pattern). e.g Gluino_13TeV_M2400"
)
options.register('ERA', '2017',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "ERA. e.g 2017A"
)

options.parseArguments()

process = cms.Process("Analyzer")

## print configuration:
print('\n')
print('CMSSW version : {}'.format(os.environ['CMSSW_VERSION']))
print('Analysis Type : {}'.format(options.TYPE))
print('Sample Type   : {}'.format(options.SAMPLE))
print('Sample Name   : {}'.format(options.NAME))
print('Sample Era    : {}'.format(options.ERA))
print('Output File   : {}'.format(options.outputFile))
print('Input Files   : {}\n'.format(options.inputFiles))

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO' #"DEBUG"
process.MessageLogger.categories.append('Analyzer')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options.inputFiles)
)

process.load("SUSYBSMAnalysis.Analyzer.HSCParticleAnalyzer_cff")
### set your configirattion here (default: python/HSCParticleAnalyzer_cff.py)

process.analyzer.TypeMode=options.TYPE

if options.SAMPLE=='isData'        : process.analyzer.SampleType=0
elif options.SAMPLE=='isBckg'      : process.analyzer.SampleType=1
elif options.SAMPLE=='isSignal'    : process.analyzer.SampleType=2
elif options.SAMPLE=='isSignalSyst': process.analyzer.SampleType=3
else: exit('Unkown SAMPLE')

process.analyzer.SampleName=options.NAME
process.analyzer.Period=options.ERA

process.analyzer.saveTree=1
process.analyzer.saveGenTree=1

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string(options.outputFile)
                                   )

process.p = cms.Path(process.analyzer)
#process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.p = cms.Path(process.analyzer*process.dump)

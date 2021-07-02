import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')

# defaults
options.outputFile = 'AnalysisTuple.root'
options.maxEvents = -1 # -1 means all events

#options.register('sampleTxtFile', 'Analysis_Samples.txt',
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.string,
#    "Sample text file"
#)
options.parseArguments()

process = cms.Process("Analyzer")

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
#process.analyzer.SampleTxtFile=options.sampleTxtFile
process.analyzer.TypeMode=0
process.analyzer.SampleType=0
process.analyzer.saveTree=6 #all saved
process.analyzer.saveGenTree=0

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string(options.outputFile)
                                   )
###process.analyzer.OutputFile = 'Data_2017_UL'

process.p = cms.Path(process.analyzer)
#process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.p = cms.Path(process.analyzer*process.dump)

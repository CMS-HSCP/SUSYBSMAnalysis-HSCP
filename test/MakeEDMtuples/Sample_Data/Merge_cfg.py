import FWCore.ParameterSet.Config as cms
process = cms.Process("MergeHLT")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.load("FWCore.MessageService.MessageLogger_cfi")
from SUSYBSMAnalysis.HSCP.HSCPVersion_cff import *

process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring( *(
XXX_INPUT_XXX
   ) )
)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.HSCPHLTDuplicate = cms.EDFilter("HSCPHLTFilter",
   RemoveDuplicates = cms.bool(True),
   TriggerProcess   = cms.string("HLT"),
   MuonTrigger1Mask    = cms.int32(0),  #Activated
   PFMetTriggerMask    = cms.int32(0),  #Activated
   L2MuMETTriggerMask  = cms.int32(0),  #Activated
)
process.DuplicateFilter = cms.Path(process.HSCPHLTDuplicate   )

process.Out = cms.OutputModule("PoolOutputModule",
     outputCommands = cms.untracked.vstring(
	 "keep *"
    ),
    fileName = cms.untracked.string('XXX_SAVEPATH_XXX'),
    SelectEvents = cms.untracked.PSet(
       SelectEvents = cms.vstring('DuplicateFilter')
    ),
)
process.endPath = cms.EndPath(process.Out)
process.schedule = cms.Schedule(process.DuplicateFilter, process.endPath)

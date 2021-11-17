import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
"file:/eos/cms/store/group/phys_muon/cericeci/7B3D3FE8-0912-164D-967D-45CE1DE4923D.root",
"file:/eos/cms/store/group/phys_muon/cericeci/D358B712-A437-F942-BABD-A3EBC4D68083.root",
"file:/eos/cms/store/group/phys_muon/cericeci/8A3751B5-890A-3D44-BFC0-B3685C8E3344.root",
"file:/eos/cms/store/group/phys_muon/cericeci/2E09E6E0-4E34-D349-9EA1-272EB69769B3.root",
"file:/eos/cms/store/group/phys_muon/cericeci/9027F6EA-74C2-8543-8B6F-23B6B7A5A9D2.root",
"file:/eos/cms/store/group/phys_muon/cericeci/945D2092-9A33-1740-BCAF-F48CD21F946D.root",
"file:/eos/cms/store/group/phys_muon/cericeci/E6436B9B-55BB-0B41-B388-A410B70C9DEB.root",
"file:/eos/cms/store/group/phys_muon/cericeci/7FCDBC49-AEFA-F34B-AAD4-5EE5A5D74959.root",

        #'file:/eos/cms/store/group/phys_muon/cericeci/HSCP_MINIAOD.root'
        #'file:/eos/cms/store/group/phys_muon/cericeci/temp/AOD_DYMM_10_6.root'
        #'file:/eos/cms/store/group/phys_muon/cericeci/temp/MINI_DYMM_10_6.root'
        #'/store/relval/CMSSW_3_2_1/RelValSingleMuPt100/GEN-SIM-RECO/MC_31X_V3-v1/0006/DC15F12B-9477-DE11-B1E0-000423D98C20.root',
    )
)

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
from Configuration.StandardSequences.Reconstruction_cff import *

process.prefer("GlobalTag")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v9', '')


process.maxEvents = cms.untracked.PSet(
	input=cms.untracked.int32(-1)
)


from UserCode.MuonTiming.muonTiming_mini_cfi import *
process.muontiming = muontiming_mini
process.p   = cms.Path(process.muontiming)
process.out = cms.OutputModule("PoolOutputModule", 
                               fileName = cms.untracked.string("muon_timing_mini.root"),
                               outputCommands=cms.untracked.vstring(
                                   'drop *',
                                   #"keep *_*muon*_*_*",
                                   "keep recoMuonTimeExtraedmValueMap_*_*_*",
                               )

)
print process
process.schedule = cms.Schedule(process.p)	
process.output_step = cms.EndPath(process.out)
process.schedule.extend([process.output_step])

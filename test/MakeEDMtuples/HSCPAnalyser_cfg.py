# To Modify ..copied from HNL analysis

from os import path as path
import FWCore.ParameterSet.Config as cms

debugLevel    = -1 

isMC_         = True
isMCSignal_   = False

jec_tag_DATA  = 'JetCorrectorParametersCollection_Summer16_23Sep2016AllV4_DATA_AK4PFchs'
jec_tag_MC    = 'JetCorrectorParametersCollection_Summer16_23Sep2016V4_MC_AK4PFchs'
jec_file_DATA = 'sqlite_file:Summer16_23Sep2016AllV4_DATA.db'
jec_file_MC   = 'sqlite_file:Summer16_23Sep2016V4_MC.db'
jec_tag       = jec_tag_MC if isMC_ else jec_tag_DATA
jec_file      = jec_file_MC if isMC_ else jec_file_DATA
algorithm     = "AK4PFchs"

GT_MC   = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
GT_DATA = '80X_dataRun2_2016SeptRepro_v7'
GT      =  GT_MC if isMC_ else GT_DATA

process = cms.Process("AnalysisProc")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 1

import FWCore.PythonUtilities.LumiList as LumiList
LumiList.LumiList().getVLuminosityBlockRange()

#from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, GT)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource", 
                            fileNames =  cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/80000/4E597432-24BE-E611-ACBB-00266CFFBFC0.root'
))
process.TFileService = cms.Service("TFileService", fileName = cms.string("Analysis_output.root"))
process.load('HNL.DisplacedAdaptiveVertexFinder.displacedInclusiveVertexing_cff')

from HNL.HeavyNeutralLeptonAnalysis.ele_Sequence_cff import addElectronSequence

addElectronSequence(process)

process.load("CondCore.CondDB.CondDB_cfi")

process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")

process.jetCorrFactors = process.updatedPatJetCorrFactors.clone(
    src = cms.InputTag("slimmedJets"),
    levels = ['L1FastJet', 
              'L2Relative', 
              'L3Absolute',
              'L2L3Residual'],
    payload = 'AK4PFchs' ) 

process.slimmedJetsJEC = process.updatedPatJets.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("jetCorrFactors"))
    )


#process.ak4PFJetsCHSL1FastL2L3 = process.ak4PFCHSJetsL1.clone(correctors = ['ak4PFL1FastL2L3Corrector'])
#process.correctJets           = cms.Sequence(process.ak4PFL1FastL2L3CorrectorChain * process.ak4PFJetsCHSL1FastL2L3)

process.HeavyNeutralLepton = cms.EDAnalyzer('HeavyNeutralLeptonAnalysis',
                                            debugLevel            = cms.int32(debugLevel),
                                            isMC                  = cms.bool(isMC_),
                                            isMCSignal            = cms.bool(isMCSignal_),
                                            muons                 = cms.InputTag("muons"),
                                            offlinePrimaryVertices= cms.InputTag("offlinePrimaryVertices"),
                                            stAlMuonTracks        = cms.InputTag("refittedStandAloneMuons"),
                                            beamSpot              = cms.InputTag("offlineBeamSpot"),
                                            isoHSCP               = cms.InputTag("HSCPIsolation","","R03"),
                                            genParticles          = cms.InputTag(""), #quale uso?
                                            collectionHSCP        = cms.InputTag("HSCParticleProducer"),
                                            dedxHit               = cms.InputTag("dedxHitInfo"),
                                            muonTOF               = cms.InputTag("muons","",""), #somwehere in the include
                                            muonTOFDT             = cms.InputTag("muons","",""),
                                            muonTOFCSC            = cms.InputTag("muons","",""),
                                            segCollectionCSC      = cms.InputTag("cscSegments"), 
                                            segCollectionDT       = cms.InputTag("dt4DSegments"),
                                            triggerResults        = cms.InputTag(""),

                                            
                                            )

process.p = cms.Path(
    process.displacedInclusiveVertexing 
    *process.ele_Sequence
    *process.jetCorrFactors
    *process.slimmedJetsJEC
    *process.HeavyNeutralLepton
    )

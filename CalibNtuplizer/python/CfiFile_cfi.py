import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('ntuple'
#     ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
      , tracks             = cms.InputTag("generalTracks")
      , dedx               = cms.InputTag("dedxHitInfo")
      , printOut = cms.untracked.int32(0)
      , GenPart            = cms.InputTag("prunedGenParticles")
      , runOnGS            = cms.bool(False)
      , stripSimLinks      = cms.InputTag("simSiStripDigis")
      , ROUList            = cms.vstring(
                                      'TrackerHitsTIBLowTof',  # hard-scatter (prompt) collection, module label g4simHits
                                      'TrackerHitsTIBHighTof',
                                      'TrackerHitsTIDLowTof',
                                      'TrackerHitsTIDHighTof',
                                      'TrackerHitsTOBLowTof',
                                      'TrackerHitsTOBHighTof',
                                      'TrackerHitsTECLowTof',
                                      'TrackerHitsTECHighTof'
                                      )

)

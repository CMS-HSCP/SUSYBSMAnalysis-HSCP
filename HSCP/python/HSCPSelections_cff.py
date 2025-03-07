
import FWCore.ParameterSet.Config as cms

HSCPSelectionEmpty = cms.PSet(
      cms.PSet(
         onlyConsiderTrack        = cms.bool(False),
         onlyConsiderMuon         = cms.bool(False),
         onlyConsiderMuonSTA      = cms.bool(False),
         onlyConsiderMuonGB       = cms.bool(False),
         onlyConsiderMuonTK       = cms.bool(False),
         onlyConsiderMTMuon       = cms.bool(False),
         onlyConsiderRpc          = cms.bool(False),
         onlyConsiderEcal         = cms.bool(False),

         minTrackHits             = cms.int32 (-1),
         minTrackP                = cms.double(-1),
         minTrackPt               = cms.double(-1),

         minDedx                  = cms.double(-1),

         minMuonP                 = cms.double(-1),
         minMuonPt                = cms.double(-1),
         minMTMuonPt              = cms.double(-1),
         minSAMuonPt              = cms.double(-1),

         maxMuTimeDtBeta          = cms.double(-1),
         minMuTimeDtNdof          = cms.double(-1),
         maxMuTimeCscBeta         = cms.double(-1),
         minMuTimeCscNdof         = cms.double(-1),
         maxMuTimeCombinedBeta    = cms.double(-1),
         minMuTimeCombinedNdof    = cms.double(-1),

         maxBetaRpc               = cms.double(-1),
         maxBetaEcal              = cms.double(-1),
      ),
)


HSCPSelectionDefault = HSCPSelectionEmpty.clone()
HSCPSelectionDefault.minTrackHits             = cms.int32(2)
HSCPSelectionDefault.minTrackPt               = cms.double(45.0)
HSCPSelectionDefault.minMuonPt                = cms.double(5.0)

HSCPSelectionMET = HSCPSelectionEmpty.clone()
HSCPSelectionMET.minTrackHits                 = cms.int32(2)
HSCPSelectionMET.minTrackPt                   = cms.double(45.0)
HSCPSelectionMET.minMuonPt                    = cms.double(5.0)

HSCPSelectionHighdEdx = HSCPSelectionDefault.clone()
HSCPSelectionHighdEdx.onlyConsiderTrack       = cms.bool(True)
HSCPSelectionHighdEdx.minDedxEstimator1       = cms.double(3.5)

HSCPSelectionHighTOF = HSCPSelectionDefault.clone()
HSCPSelectionHighTOF.onlyConsiderMuon         = cms.bool(True)
HSCPSelectionHighTOF.maxMuTimeDtBeta          = cms.double(0.9)

HSCPSelectionMTMuonOnly = HSCPSelectionEmpty.clone()
HSCPSelectionMTMuonOnly.onlyConsiderMTMuon    = cms.bool(True)
HSCPSelectionMTMuonOnly.minMTMuonPt           = cms.double(70.0)


HSCPSelectionSAMuonOnly = HSCPSelectionEmpty.clone()
HSCPSelectionSAMuonOnly.onlyConsiderMuonSTA    = cms.bool(True)
HSCPSelectionSAMuonOnly.minSAMuonPt           = cms.double(70.0)

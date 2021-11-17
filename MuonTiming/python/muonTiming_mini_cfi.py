import FWCore.ParameterSet.Config as cms

from RecoMuon.MuonIdentification.MuonTimingFiller_cfi import *

muontiming_mini = cms.EDProducer('MuonTimingProducer_Mini',
  TimingFillerBlock,
  MuonCollection = cms.InputTag("slimmedMuons"),
)

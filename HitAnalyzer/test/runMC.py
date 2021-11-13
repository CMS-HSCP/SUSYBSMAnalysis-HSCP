import FWCore.ParameterSet.Config as cms

f=open("list.txt", "r")
files =[]
files.append([x.strip() for x in f.readlines()])

process = cms.Process("Demo")



### standard MessageLoggerConfiguration
process.load("FWCore.MessageService.MessageLogger_cfi")

### Standard Configurations

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')




### Track refitter specific stuff
process.load("RecoTracker.TrackProducer.TrackRefitters_cff") #the correct one




process.TrackRefitter.NavigationSchool = ''
process.TrackRefitter.Fitter = 'FlexibleKFFittingSmoother'

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use#
    fileNames = cms.untracked.vstring(
        #files[:][0]
        'file:54159F2D-02ED-E711-ADB2-0242AC130002.root'
    )
)
process.muFilter = cms.EDFilter("MuonFilter",

)

process.demo = cms.EDAnalyzer('HitAnalyzer'
)


process.p = cms.Path(process.TrackRefitter+process.muFilter+process.demo)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v6', '')






#jsonFile='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_27Jan2016ReReco_Collisions15_ZeroTesla_25ns_JSON_MuonPhys.txt'



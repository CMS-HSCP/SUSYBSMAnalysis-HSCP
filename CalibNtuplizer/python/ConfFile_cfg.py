import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:/opt/sbg/cms/ui3_data1/ccollard/HSCP/test_aod_mc/FEAE43F8-6D9C-E711-A831-3417EBE645A9.root'  # AODSIM
#        'file:/opt/sbg/cms/ui3_data1/ccollard/HSCP/test_aod_mc/0E8CD6D5-EC97-E711-A564-28924A33AF26.root'   # GEN-SIM-RECO
#        'file:root://xrootd-cms.infn.it//store/user/jpriscia/HSCP_2017/SingleMuon/2017Apr25/180427_073624/0000/HSCP_1.root'  #2017 data
#         '/store/mc/RunIIFall17DRPremix/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/94X_mc2017_realistic_v11-v1/90000/B8C3405B-DD34-E811-938A-009C02AAB554.root'
#        'file:/opt/sbg/cms/ui3_data1/ccollard/HSCP/ttbar/step3.root'   # GEN-SIM-RECO
#1
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/E6F831BC-43D2-E811-AE09-0CC47AD99176.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/E420CC7E-43D2-E811-8A04-48D539F38676.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/D8633B09-C3D1-E811-9ACB-0025905C2CA4.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/B6FA8401-43D2-E811-A41F-0025905C5476.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/A4BD6428-C3D1-E811-8717-F02FA768CF56.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/9AA80E25-BED1-E811-A4A4-48FD8EE73A79.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/607E1735-C3D1-E811-BE2F-48FD8EE739A9.root' #AODSIM MinBias
#2
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/5E0C834C-B4D1-E811-A038-0CC47AF9B2CA.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/5C24C04E-B4D1-E811-8CD3-0025904B7C26.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/58F5702B-B9D1-E811-9DC4-0CC47AFB7F5C.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/50648F41-B9D1-E811-ABED-0025905C542C.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/4A2EC479-64D2-E811-A1C0-0025905C3DF6.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/4891CDBD-B6D1-E811-B837-0025905C3D3E.root'
#3
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/3A76C9BD-B6D1-E811-B975-0025905C54FE.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/380943E8-C2D1-E811-ACE3-0CC47AA992B4.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/34CB0E8F-64D2-E811-82FC-48FD8EE73ABB.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/60000/30688C64-64D2-E811-AC26-0CC47AD98A92.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/F6F5825F-D7AB-E811-B596-00010100195C.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/F2CB6CB2-64AB-E811-AEEA-24BE05CEDCD1.root'
#4
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/F2719A60-D2AC-E811-8735-0001010022C3.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/F09F7B60-2EAD-E811-9B18-F01FAFE5CA86.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/EA2F076D-4EAC-E811-8962-000101003817.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/D835713D-4CAD-E811-8693-5065F3817281.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/D4AD5A94-19AE-E811-A1F3-901B0E5427B6.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/MinBias_TuneMBR_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/D2854275-12AD-E811-B8E8-0001010024D9.root' 
#         gluino
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/HSCPgluino_M-2400_TuneCP5_13TeV-pythia8/AODSIM/PU2017_HSCP1_94X_mc2017_realistic_v11-v2/110000/2A1E4FBE-BE74-E811-AF2E-509A4C8449D1.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/HSCPgluino_M-2400_TuneCP5_13TeV-pythia8/AODSIM/PU2017_HSCP1_94X_mc2017_realistic_v11-v2/100000/E2D56EB9-D46F-E811-AB47-90E2BAD1BDF0.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/HSCPgluino_M-2400_TuneCP5_13TeV-pythia8/AODSIM/PU2017_HSCP1_94X_mc2017_realistic_v11-v2/100000/DC9D25F0-1F70-E811-96FC-064A2A000168.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/HSCPgluino_M-2400_TuneCP5_13TeV-pythia8/AODSIM/PU2017_HSCP1_94X_mc2017_realistic_v11-v2/100000/DC74B110-F972-E811-9A34-009C02AABEB8.root',
#         'root://cmsxrootd.fnal.gov///store/mc/RunIIFall17DRPremix/HSCPgluino_M-2400_TuneCP5_13TeV-pythia8/AODSIM/PU2017_HSCP1_94X_mc2017_realistic_v11-v2/100000/DC20A0B2-F872-E811-8225-34E6D7E05F01.root'
#         nouveau ttbar semilept
         '/store/mc/RunIIFall17DRPremix/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/94X_mc2017_realistic_v11-v1/00001/1080AC42-DA25-E811-9587-7CD30AC03054.root',
         #'/store/mc/RunIIFall17DRPremix/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/94X_mc2017_realistic_v11-v1/00001/40E47FB4-F625-E811-A892-0026B9277A0B.root',
         '/store/mc/RunIIFall17DRPremix/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/94X_mc2017_realistic_v11-v1/00001/4A1D4643-C226-E811-8CBC-BC305B390A80.root',
         '/store/mc/RunIIFall17DRPremix/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/94X_mc2017_realistic_v11-v1/00001/4AC3CB50-C625-E811-B340-BC305B390AA7.root',
         '/store/mc/RunIIFall17DRPremix/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/94X_mc2017_realistic_v11-v1/00001/4EEC88FE-EC25-E811-83DF-0023AEEEB6FA.root',
         '/store/mc/RunIIFall17DRPremix/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/94X_mc2017_realistic_v11-v1/00001/5C56AB98-1C26-E811-A005-0026B92779BD.root',
         '/store/mc/RunIIFall17DRPremix/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/94X_mc2017_realistic_v11-v1/00001/68E0EF6B-0F26-E811-B628-7CD30AC03054.root',
         '/store/mc/RunIIFall17DRPremix/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/94X_mc2017_realistic_v11-v1/00001/78FD9F66-FB25-E811-8A78-28924A35056E.root'
    )
)

process.stage = cms.EDAnalyzer('ntuple'
     , tracks             = cms.InputTag("generalTracks")
     , dedx               = cms.InputTag("dedxHitInfo")
     , printOut           = cms.untracked.int32(-1)
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

process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService",
#     fileName = cms.string('/opt/sbg/cms/ui3_data1/ccollard/HSCP/ttbar/ntuple_ttbar_aod.root')
#     fileName = cms.string('test_ttbar_aod.root')
#     fileName = cms.string('test_minbias_aod.root')
#     fileName = cms.string('/opt/sbg/cms/ui3_data1/ccollard/HSCP/minbias_check/test_minbias_aod4.root')
#     fileName = cms.string('/opt/sbg/cms/ui3_data1/ccollard/HSCP/gluino2400_aod/test_gluino_aod1.root')
     fileName = cms.string('/opt/sbg/cms/ui3_data1/ccollard/HSCP/ttbar/test_ttbar_aod2.root')
 )

process.p = cms.Path(process.stage)
#process.p = cms.Path(process.dump)


## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.MessageLogger.cerr.default.limit = 10

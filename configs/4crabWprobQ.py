from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Analysis_SingleMuon_UL2018C_wProbQ'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HSCParticleProducerAnalyzer_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.JobType.inputFiles = ['SUSYBSMAnalysis/HSCP/data/Data13TeVGains_v2.root','SUSYBSMAnalysis/HSCP/data/dEdxTemplate_harm2_SO_in_noC_CCC_MG_2017B.root','SUSYBSMAnalysis/HSCP/data/CMS_GeomTree.root','SUSYBSMAnalysis/HSCP/data/MuonTimeOffset.txt']

config.Data.inputDataset = '/SingleMuon/tvami-crab_PrivateReAOD_2018_SingleMuon_AOD_v1-f31fbc775312df6a6c2c9ecba57d4788/USER'
config.Data.inputDBS = 'phys03'
    #config.Data.splitting = 'Automatic'
    #config.Data.splitting = 'LumiBased'
    #config.Data.unitsPerJob = 1 #20
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

    #config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
    #config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
    #config.Data.lumiMask = '/afs/cern.ch/work/d/dapparu/private/thesis/hscp/CMSSW_10_6_20/src/SUSYBSMAnalysis/crab_projects/crab_EDMProd_MET_UL2017B_297047_297484/results/processedLumis.json'

    #config.Data.runRange = '299368-300284'
    #config.Data.runRange = '297047-306462'

config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.Data.outLFNDirBase = '/store/user/tvami/HSCP'

config.Site.storageSite = 'T2_HU_Budapest'

from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Analysis_SingleMuon_UL2017CEra_wProbQ_newMethod_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HSCParticleProducerAnalyzer_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.JobType.inputFiles = ['Data13TeVGains_v2.root','dEdxTemplate_harm2_SO_in_noC_CCC_MG_2017B.root','CMS_GeomTree.root','MuonTimeOffset.txt','Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt']


config.Data.inputDataset = '/SingleMuon/Run2017C-09Aug2019_UL2017-v1/AOD'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'Automatic'
    #config.Data.splitting = 'LumiBased'
    #config.Data.unitsPerJob = 1 #20
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 1
config.Data.publication = True 
config.Data.outputDatasetTag = config.General.requestName
config.Data.outLFNDirBase = '/store/user/tvami/HSCP'
config.Site.storageSite = 'T2_HU_Budapest'

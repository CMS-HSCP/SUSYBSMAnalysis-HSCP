from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Analysis_SingleMuon_UL20177C_wProbQ_NewMethod_v2'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HSCParticleProducerAnalyzer_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.JobType.inputFiles = ['Data13TeVGains_v2.root','dEdxTemplate_harm2_SO_in_noC_CCC_MG_2017B.root','CMS_GeomTree.root','MuonTimeOffset.txt']


config.Data.inputDataset = '/SingleMuon/tvami-crab_HLT2AOD_2017_SingleMuon_v4-7cee3b736accffaaa03458294e357e5d/USER'
#config.Data.inputDataset = '/SingleMuon/tvami-crab_PrivateReAOD_2018_SingleMuon_AOD_v1-f31fbc775312df6a6c2c9ecba57d4788/USER'
config.Data.inputDBS = 'phys03'
    #config.Data.splitting = 'Automatic'
    #config.Data.splitting = 'LumiBased'
    #config.Data.unitsPerJob = 1 #20
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.publication = True 
config.Data.outputDatasetTag = config.General.requestName
config.Data.outLFNDirBase = '/store/user/tvami/HSCP'
config.Site.storageSite = 'T2_HU_Budapest'

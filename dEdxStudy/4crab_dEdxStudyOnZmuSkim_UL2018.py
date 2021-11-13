from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Analysis_SingleMuon_UL2018C_dEdxSkimmer_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'dEdxSkimmer.py'
config.JobType.outputFiles = ['dEdxSkim.root']
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 3000
config.JobType.maxMemoryMB = 4000

config.Data.inputDataset = '/SingleMuon/Run2018C-ZMu-12Nov2019_UL2018-v2/RAW-RECO'
#config.Data.inputDBS = 'phys03'
#config.Data.splitting = 'Automatic'
#config.Data.splitting = 'LumiBased'
#config.Data.unitsPerJob = 1 #20
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.Data.outLFNDirBase = '/store/user/tvami/HSCP'

config.Site.storageSite = 'T2_HU_Budapest'

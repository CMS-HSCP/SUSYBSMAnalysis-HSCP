from CRABClient.UserUtilities import config
config = config()

config.section_('General')
config.General.requestName = 'Analysis_HSCPgluino_M_1800_CodeV4p6_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HSCParticleProducerAnalyzer_SignalMC_cfg.py'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.maxJobRuntimeMin = 3000
config.JobType.maxMemoryMB = 3000
config.JobType.inputFiles = ['dEdxTemplate_harm2_SO_in_noC_CCC_MG_2017B.root','CMS_GeomTree.root','MuonTimeOffset.txt']

config.section_('Data')
config.Data.inputDataset = '/HSCPgluino_M_1800/tvami-crab_PrivateHSCP_2018_Gluino_Mass1800_AOD_NoPU_v1-1363d434c2f54354c4b307a55f861103/USER'
config.Data.inputDBS = 'phys03'
#config.Data.splitting = 'Automatic'
config.Data.splitting = 'LumiBased'
    #config.Data.unitsPerJob = 1 #20
#config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 50
config.Data.publication = True 
config.Data.outputDatasetTag = config.General.requestName
config.Data.outLFNDirBase = '/store/user/tvami/HSCP'
config.Data.ignoreLocality = True

config.section_('Site')
config.Site.whitelist = ['T2_DE_DESY','T2_FR_IPHC','T2_CH_CERN','T2_IT_Bari','T1_IT_*','T2_US_*']
config.Site.storageSite = 'T2_HU_Budapest'


from CRABClient.UserUtilities import config
config = config()

config.section_('General')
#config.General.requestName = 'Analysis_2018_MyHSCPgluino_NoRew_NoCust_M_1800_woProbQ_CodeV27p0_v1'
config.General.requestName = 'Analysis_2018_MyHSCPgluino_NoRew_M_1800_woProbQ_CodeV27p5_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HSCParticleProducerAnalyzer_2018_MySignalMC_woProbQ_cfg.py'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.maxJobRuntimeMin = 3000
config.JobType.maxMemoryMB = 3500
config.JobType.inputFiles = ['templateMC.root','MuonTimeOffset.txt','Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt']

config.section_('Data')
config.Data.inputDBS = 'phys03'
#config.Data.inputDataset = '/HSCPgluino_M_1800/tvami-crab_PrivateHSCP_2018_Gluino_Mass1800_DIGI2AOD_NoPU_v3-0bfbad649d32c05924b5bfd4b5874292/USER'
config.Data.inputDataset = '/HSCPgluino_M_1800/tvami-crab_PrivateHSCP_2018_Gluino_Mass1800_DIGI2AOD_NoPU_NoRew_v1-e88ed05559facf15902f2fe2f6141103/USER'
#config.Data.inputDataset = '/HSCPgluino_M_1800/tvami-crab_PrivateHSCP_2018_Gluino_Mass1800_DIGI2AOD_NoPU_NoRewGT_v1-a007e86beae864f0107e60af996b4558/USER'
#config.Data.inputDataset = '/HSCPgluino_M_1800/tvami-crab_PrivateHSCP_2018_Gluino_Mass1800_DIGI2AOD_NoPU_v3-0bfbad649d32c05924b5bfd4b5874292/USER'
config.Data.splitting = 'LumiBased'
#config.Data.unitsPerJob = 1 #20
#config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 50
#config.Data.totalUnits = config.Data.unitsPerJob * 1000
config.Data.publication = True
config.Data.outputDatasetTag = config.General.requestName
config.Data.outLFNDirBase = '/store/user/tvami/HSCP'
config.Data.ignoreLocality = True
config.Data.runRange = '0'

config.section_('Site')
config.Site.whitelist = ['T2_DE_DESY','T2_FR_IPHC','T2_CH_CERN','T2_IT_Bari','T1_IT_*','T2_US_*']
config.Site.storageSite = 'T2_HU_Budapest'
  

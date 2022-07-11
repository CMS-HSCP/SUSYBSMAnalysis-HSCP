
from CRABClient.UserUtilities import config
config = config()

config.section_('General')
config.General.requestName = 'Analysis_ROVIDMINTA_wProbQ_CodeVVERZIO_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HSCParticleProducerAnalyzer_data_wProbQ_cfg.py'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.maxJobRuntimeMin = 3000
config.JobType.maxMemoryMB = 3500
config.JobType.inputFiles = ['SUSYBSMAnalysis/HSCP/data/CorrFact2017PixL1.txt','SUSYBSMAnalysis/HSCP/data/CorrFact2017PixL2.txt','SUSYBSMAnalysis/HSCP/data/CorrFact2017PixL3.txt','SUSYBSMAnalysis/HSCP/data/CorrFact2017PixL4.txt','SUSYBSMAnalysis/HSCP/data/CorrFact2017PixR1.txt','SUSYBSMAnalysis/HSCP/data/CorrFact2017PixR2.txt','SUSYBSMAnalysis/HSCP/data/template_2017C.root','MuonTimeOffset.txt','Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt']

config.section_('Data')
config.Data.inputDataset = 'MINTA'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 50
config.Data.publication = True
config.Data.outputDatasetTag = config.General.requestName
config.Data.outLFNDirBase = '/store/user/tvami/HSCP'
config.Data.ignoreLocality = True

config.section_('Site')
config.Site.whitelist = ['T2_DE_DESY','T2_FR_IPHC','T2_CH_CERN','T2_IT_Bari','T1_IT_*','T2_US_*']
config.Site.storageSite = 'T2_HU_Budapest'
#config.Site.storageSite = 'T3_US_FNALLPC'
  

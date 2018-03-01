from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.requestName = 'UCLouvain_process'
config.General.transferOutputs = True
config.section_('JobType')
config.JobType.psetName = '../HSCParticleProducer_Data_cfg.py'
config.JobType.pluginName = 'Analysis'
#config.JobType.inputFiles = ['Data7TeV_Deco_SiStripDeDxMip_3D_Rcd.db']
#config.JobType.pyCfgParams = ['globalTag=GR_P_V42B::All']
config.JobType.outputFiles = ['HSCP.root']
config.section_('Data')
config.Data.inputDataset = '/SingleMuon/Run2015B-PromptReco-v1/AOD'
config.Data.publication = False
config.Data.unitsPerJob = 6 # No. of lumisections each job works through
config.Data.splitting = 'LumiBased'
config.Data.inputDBS = 'global'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/DCSOnly/json_DCSONLY_Run2015B.txt'
#config.Data.runRange = '100000-900000'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_BE_UCL'

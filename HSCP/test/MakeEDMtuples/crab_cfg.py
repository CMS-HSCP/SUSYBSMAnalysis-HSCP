from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'HSCP_data_filtering_METF' #o qualunque cosa tu voglia, meglio specificare il dataset
config.General.workArea = 'crab_jobs' 
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/home/ucl/cp3/jpriscia/CMSSW_9_4_3/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/HSCParticleProducerMET2017_cfg.py' #metti il full path, e` piu` sicuro
#config.JobType.pyCfgParams = ['GTAG=94X_dataRun2_ReReco_EOY17_v2']

config.Data.inputDataset = '/MET/Run2017F-17Nov2017-v1/AOD' #il dataset
#config.Data.inputDBS = 'global' #solo per dataset privati
config.Data.splitting = 'LumiBased' #'FileBased' for MC
config.Data.unitsPerJob = 15 #lumisections per job (dati) file per job (MC)
config.Data.outLFNDirBase = '/store/user/%s/HSCP_2017/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = '2017Apr25' #qualcosa per ricordarti chi sono
config.Data.lumiMask = '/home/ucl/cp3/jpriscia/CMSSW_9_4_3/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/Sample_Data/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt' #solo per dati, full path al json file

#config.Site.whitelist   = ['T2_BE_UCL']
config.Site.storageSite = 'T2_BE_UCL'

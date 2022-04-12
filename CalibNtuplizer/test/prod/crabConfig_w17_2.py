from WMCore.Configuration import Configuration
# More details here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial

config = Configuration()

config.section_('General')
config.General.requestName = 'w17_2'
#config.General.workArea = '.'
#config.General.instance = 'private'
config.General.transferOutputs = True
config.General.transferLogs = False
#config.General.serverUrl = 'You need to set the CRAB3 server URL'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_aod_mc17.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ['minbias_template_uncorr_iter1.root','minbias_template_corr_iter1.root']

config.section_('Data')
#config.Data.inputDataset = '/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM'
config.Data.inputDataset = '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM'
config.Data.allowNonValidInputDataset = True # FIXME
#config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'LumiBased'
#config.Data.unitsPerJob = 25  
config.Data.splitting = 'Automatic'
#config.Data.splitting = 'FileBased'  # special case for memory issues
#config.Data.unitsPerJob = 100
#config.Data.lumiMask = 'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'
#config.Data.runRange = '315257,315259,315264,315265,315267,315270,315339,315357,315361,315363,315365,315366,315420,315489,315490'
config.Data.publication = False
config.Data.outputDatasetTag = 'Prod2'
config.Data.outLFNDirBase = '/store/user/ccollard/HSCP/prodOct2021_CMSSW_10_6_27/'
#config.Data.ignoreLocality = True  # to be used only if  use the whitelist

config.section_('Site')
config.Site.storageSite = 'T2_FR_IPHC'
#config.Site.storageSite = 'T2_CH_CERN'
#config.Site.blacklist = ['T2_IT_Legnaro']
#config.Site.whitelist = ['T2_FR_IPHC']

config.section_('User')


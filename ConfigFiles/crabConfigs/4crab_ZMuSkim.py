import os
import glob

from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'SingleMuonZMuSkim_2018C_BasedHSCPNtuple_v1'

config.section_('JobType')
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_ZMuSkim_data18_Tamas.py'
config.JobType.inputFiles = ['minbias_template_corr_iter1.root','minbias_template_uncorr_iter1.root']
config.JobType.outputFiles = ['nt_data_ZMuSkim.root']
config.JobType.disableAutomaticOutputCollection = True
config.JobType.maxMemoryMB = 5000
# config.JobType.numCores = 8

config.section_('Data')
#config.Data.inputDBS = 'phys03'
config.Data.inputDataset = '/SingleMuon/Run2018C-ZMu-12Nov2019_UL2018-v2/RAW-RECO'
config.Data.outLFNDirBase = '/store/user/tvami/HSCP/'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 30 
config.Data.publication = True

config.section_('Site')
config.Site.storageSite = 'T2_HU_Budapest'
#config.Site.storageSite = 'T3_US_FNALLPC'

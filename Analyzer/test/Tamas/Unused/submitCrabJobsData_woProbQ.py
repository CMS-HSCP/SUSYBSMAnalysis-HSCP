import sys, os, time, re
import numpy as np
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog codeVersion")
(opt,args) = parser.parse_args()

datasetList = [
#"/SingleMuon/Run2017C-09Aug2019_UL2017-v1/AOD",
#"/MET/Run2017C-09Aug2019_UL2017_rsb-v1/AOD",
#"/SingleMuon/Run2018A-15Feb2022_UL2018-v1/AOD",
#"/SingleMuon/Run2018B-15Feb2022_UL2018-v1/AOD",
"/SingleMuon/Run2018C-15Feb2022_UL2018-v1/AOD",
#"/SingleMuon/Run2018D-15Feb2022_UL2018-v1/AOD",
#"/MET/Run2018C-15Feb2022_UL2018-v1/AOD",
]

codeVersion = sys.argv[1]
#just the number, like 18p2

if not os.path.exists("submittedConfigs"): os.makedirs("submittedConfigs")

if not os.path.exists("4crab_Template_Data.py"):
  TEMPLATE = '''
from CRABClient.UserUtilities import config
config = config()

config.section_('General')
config.General.requestName = 'Analysis_2018_ROVIDMINTA_CodeVVERZIO_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HSCParticleProducerAnalyzer_data_cfg.py'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.maxJobRuntimeMin = 3000
config.JobType.maxMemoryMB = 3500
config.JobType.inputFiles = ['SUSYBSMAnalysis/HSCP/data/CorrFact2018PixL1.txt','SUSYBSMAnalysis/HSCP/data/CorrFact2018PixL2.txt','SUSYBSMAnalysis/HSCP/data/CorrFact2018PixL3.txt','SUSYBSMAnalysis/HSCP/data/CorrFact2018PixL4.txt','SUSYBSMAnalysis/HSCP/data/CorrFact2018PixR1.txt','SUSYBSMAnalysis/HSCP/data/CorrFact2018PixR2.txt','SUSYBSMAnalysis/HSCP/data/template_2017C.root','MuonTimeOffset.txt','Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt']

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
config.Site.blacklist = ['T2_US_Nebraska']
config.Site.storageSite = 'T2_HU_Budapest'
#config.Site.storageSite = 'T3_US_FNALLPC'
  '''

  with open("4crab_Template_Data.py", "w") as text_file:
      text_file.write(TEMPLATE)

for i in datasetList:
  print("Submit for sample "+i)
  os.system("cp 4crab_Template_Data.py 4crab_toSubmit_Data.py")
  replaceVERZIO = "sed -i 's/VERZIO/"+codeVersion+"/g' 4crab_toSubmit_Data.py"
  os.system(replaceVERZIO)
  shortSampleName = i[1:(i.find('-'))].replace("/","_")
  replaceROVIDMINTA = "sed -i 's/ROVIDMINTA/"+shortSampleName+"/g' 4crab_toSubmit_Data.py"
  os.system(replaceROVIDMINTA)
  replaceMINTA = "sed -i 's/MINTA/"+i.replace("/","\/")+"/g' 4crab_toSubmit_Data.py"
  os.system(replaceMINTA)
  os.system("crab submit -c 4crab_toSubmit_Data.py")

  os.system("mv 4crab_toSubmit_Data.py submittedConfigs/.")



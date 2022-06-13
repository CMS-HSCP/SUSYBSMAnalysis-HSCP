import sys, os, time, re
import numpy as np
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog codeVersion")
(opt,args) = parser.parse_args()

datasetList = [
"/HSCPgluino_M-1000_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/HSCPgluino_M-100_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/HSCPgluino_M-1400_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/HSCPgluino_M-1600_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/HSCPgluino_M-1800_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/HSCPgluino_M-2000_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/HSCPgluino_M-200_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/HSCPgluino_M-2200_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/HSCPgluino_M-2400_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/HSCPgluino_M-2600_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/HSCPgluino_M-400_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/HSCPgluino_M-500_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/HSCPgluino_M-800_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
]

codeVersion = sys.argv[1]
#just the number, like 18p2

didVoms = input("Push enter if you alread did voms-proxy-init -rfc -voms cms -valid 192:00 otherwise say no and do it\n")
if(didVoms):
 sys.exit()

if not os.path.exists("submittedConfigs"): os.makedirs("submittedConfigs")

if not os.path.exists("4crab_Template_wProbQ.py"):
  TEMPLATE = '''
from CRABClient.UserUtilities import config
config = config()

config.section_('General')
config.General.requestName = 'Analysis_2018_ROVIDMINTA_wProbQ_CodeVVERZIO_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.instance = 'preprod'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HSCParticleProducerAnalyzer_2018_mc_wProbQ_cfg.py'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.maxJobRuntimeMin = 3000
config.JobType.maxMemoryMB = 3500
config.JobType.inputFiles = ['dEdxTemplate_harm2_SO_in_noC_CCC_MG_2017B.root','CMS_GeomTree.root','MuonTimeOffset.txt','Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt']

config.section_('Data')
config.Data.inputDataset = 'MINTA'
#config.Data.inputDBS = 'phys03'
#config.Data.splitting = 'Automatic'
config.Data.splitting = 'LumiBased'
    #config.Data.unitsPerJob = 1 #20
#config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 50
config.Data.totalUnits = config.Data.unitsPerJob * 1000
config.Data.publication = True
config.Data.outputDatasetTag = config.General.requestName
config.Data.outLFNDirBase = '/store/user/tvami/HSCP'
config.Data.ignoreLocality = True
config.Data.runRange = '0'

config.section_('Site')
config.Site.whitelist = ['T2_DE_DESY','T2_FR_IPHC','T2_CH_CERN','T2_IT_Bari','T1_IT_*','T2_US_*']
config.Site.storageSite = 'T2_HU_Budapest'
  '''

  with open("4crab_Template_woProbQ.py", "w") as text_file:
      text_file.write(TEMPLATE)

for i in datasetList:
  print("Submit for sample "+i)
  os.system("cp 4crab_Template_wProbQ.py 4crab_toSubmit_wProbQ.py")
  replaceVERZIO = "sed -i 's/VERZIO/"+codeVersion+"/g' 4crab_toSubmit_wProbQ.py"
  os.system(replaceVERZIO)
  shortSampleName = i[1:(i.find('TuneCP5'))-1]
  replaceROVIDMINTA = "sed -i 's/ROVIDMINTA/"+shortSampleName+"/g' 4crab_toSubmit_wProbQ.py"
  os.system(replaceROVIDMINTA)
  replaceMINTA = "sed -i 's/MINTA/"+i.replace("/","\/")+"/g' 4crab_toSubmit_wProbQ.py"
  os.system(replaceMINTA)
  os.system("crab submit -c 4crab_toSubmit_wProbQ.py")
  os.system("mv 4crab_toSubmit_wProbQ.py submittedConfigs/.")



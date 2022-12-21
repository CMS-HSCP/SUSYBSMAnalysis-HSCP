import sys, os, time, re
import numpy as np
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog codeVersion")
(opt,args) = parser.parse_args()

datasetList = [
"/SingleMuon/Run2017B-15Feb2022_UL2017-v1/AOD",
"/SingleMuon/Run2017C-15Feb2022_UL2017-v1/AOD",
"/SingleMuon/Run2017D-15Feb2022_UL2017-v1/AOD",
"/SingleMuon/Run2017E-15Feb2022_UL2017-v1/AOD",
"/SingleMuon/Run2017F-15Feb2022_UL2017-v1/AOD",
#"/SingleMuon/Run2017G-15Feb2022_UL2017-v1/AOD", # 5 TeV
#"/SingleMuon/Run2017H-15Feb2022_UL2017-v1/AOD", # lowPU
"/SingleMuon/Run2018A-15Feb2022_UL2018-v1/AOD",
"/SingleMuon/Run2018B-15Feb2022_UL2018-v1/AOD", 
"/SingleMuon/Run2018C-15Feb2022_UL2018-v1/AOD",
"/SingleMuon/Run2018D-15Feb2022_UL2018-v1/AOD"
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
config.General.requestName = 'Analysis_ROVIDMINTA_CodeVVERZIO_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HSCParticleProducerAnalyzer_master_cfg.py'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.maxJobRuntimeMin = 3000
config.JobType.maxMemoryMB = 4000
config.JobType.inputFiles = ['SUSYBSMAnalysis/HSCP/data/CorrFactEVADPixL1.txt','SUSYBSMAnalysis/HSCP/data/CorrFactEVADPixL2.txt','SUSYBSMAnalysis/HSCP/data/CorrFactEVADPixL3.txt','SUSYBSMAnalysis/HSCP/data/CorrFactEVADPixL4.txt','SUSYBSMAnalysis/HSCP/data/CorrFactEVADPixR1.txt','SUSYBSMAnalysis/HSCP/data/CorrFactEVADPixR2.txt','SUSYBSMAnalysis/HSCP/data/GiTEMPLATE','MuonTimeOffset.txt']
config.JobType.pyCfgParams = ['GTAG=106X_dataRun2_v36', 'SAMPLE=isData', 'YEAR=EVAD', 'ERA=IDOSZAK']

config.section_('Data')
config.Data.inputDataset = 'MINTA'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 50
#MASZK2018config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
#MASZK2017config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
config.Data.publication = True
config.Data.outputDatasetTag = config.General.requestName
config.Data.outLFNDirBase = '/store/user/tvami/HSCP'
config.Data.ignoreLocality = True

config.section_('Site')
config.Site.whitelist = ['T2_DE_DESY','T2_CH_CERN','T2_IT_Bari','T1_IT_*','T2_US_*', 'T3_US_FNALLPC','T2_HU_Budapest','T2_FR_*', 'T2_UK_London_IC']
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
  if ("2017" in i) :
    EVAD = "2017"
    replaceMASZK = "sed -i 's/#MASZK2017//g' 4crab_toSubmit_Data.py"
  elif ("2018" in i) :
    EVAD = "2018"
    replaceMASZK = "sed -i 's/#MASZK2018//g' 4crab_toSubmit_Data.py"
  replaceEVAD  = "sed -i 's/EVAD/"+EVAD+"/g' 4crab_toSubmit_Data.py"
  os.system(replaceEVAD)
  os.system(replaceMASZK)
  if ("A-" in i) :
  	IDOSZAK = "A"
  elif ("B-" in i) :
  	IDOSZAK = "B"
  elif ("C-" in i) :
  	IDOSZAK = "C"
  elif ("D-" in i) :
  	IDOSZAK = "D"
  elif ("E-" in i) :
  	IDOSZAK = "E"
  elif ("F-" in i) :
  	IDOSZAK = "F"
  elif ("G-" in i) :
  	IDOSZAK = "G"
  elif ("H-" in i) :
  	IDOSZAK = "H"
  if ("AODSIM" in i ):
        IDOSZAK = "MC"
  replaceIDOSZAK  = "sed -i 's/IDOSZAK/"+IDOSZAK+"/g' 4crab_toSubmit_Data.py"
  os.system(replaceIDOSZAK)
  	
  replaceGiTEMPLATE = "sed -i 's/GiTEMPLATE/template_"+EVAD+IDOSZAK+"_v2.root/g' 4crab_toSubmit_Data.py"
  os.system(replaceGiTEMPLATE)

  os.system("crab submit -c 4crab_toSubmit_Data.py")
  os.system("mv 4crab_toSubmit_Data.py submittedConfigs/"+codeVersion+"_"+EVAD+IDOSZAK+".py")


os.system("rm 4crab_Template_Data.py")

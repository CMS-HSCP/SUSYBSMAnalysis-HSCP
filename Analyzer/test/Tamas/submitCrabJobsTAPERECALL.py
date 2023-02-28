import sys, os, time, re
import numpy as np
from threading import Thread
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog codeVersion")
(opt,args) = parser.parse_args()

datasetList = [
#"/SingleMuon/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/AOD",
#"/SingleMuon/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/AOD",
#"/SingleMuon/Run2016C-21Feb2020_UL2016_HIPM-v1/AOD",
#"/SingleMuon/Run2016D-21Feb2020_UL2016_HIPM-v1/AOD",
#"/SingleMuon/Run2016E-21Feb2020_UL2016_HIPM-v1/AOD",
#"/SingleMuon/Run2016F-21Feb2020_UL2016_HIPM-v1/AOD",
#"/SingleMuon/Run2016F-21Feb2020_UL2016-v1/AOD",
#"/SingleMuon/Run2016G-21Feb2020_UL2016-v1/AOD",
#"/SingleMuon/Run2016H-21Feb2020_UL2016-v1/AOD",
"/SingleMuon/Run2017B-15Feb2022_UL2017-v1/AOD",
"/SingleMuon/Run2017C-15Feb2022_UL2017-v1/AOD",
"/SingleMuon/Run2017D-15Feb2022_UL2017-v1/AOD",
"/SingleMuon/Run2017E-15Feb2022_UL2017-v1/AOD",
"/SingleMuon/Run2017F-15Feb2022_UL2017-v1/AOD",
"/SingleMuon/Run2018A-15Feb2022_UL2018-v1/AOD",
"/SingleMuon/Run2018B-15Feb2022_UL2018-v1/AOD", 
"/SingleMuon/Run2018C-15Feb2022_UL2018-v1/AOD",
"/SingleMuon/Run2018D-15Feb2022_UL2018-v1/AOD"
"/QCD_Pt-50To80_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM",
"/QCD_Pt-80To120_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM",
"/QCD_Pt-120To170_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM",
"/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM",
"/ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/ZToMuMu_M-120To200_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/ZToMuMu_M-200To400_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/ZToMuMu_M-400To800_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/ZToMuMu_M-800To1400_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/ZToMuMu_M-1400To2300_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM",
"/ZToMuMu_M-2300To3500_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/ZToMuMu_M-3500To4500_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM",
"/ZToMuMu_M-4500To6000_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
"/ZToMuMu_M-6000ToInf_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM",
]

codeVersion = sys.argv[1]
#just the number, like 18p2

if not os.path.exists("submittedConfigs"): os.makedirs("submittedConfigs")

if not os.path.exists("4crab_TemplateMTtape.py"):
  TEMPLATE = '''
from CRABClient.UserUtilities import config
config = config()

config.section_('General')
config.General.requestName = 'Analysis_ROVIDMINTA_CodeVVERZIO_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HSCParticleProducerAnalyzer_master_cfg_TAPERECALL.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 300
config.JobType.maxMemoryMB = 2500
config.JobType.inputFiles = ['SUSYBSMAnalysis/HSCP/data/CorrFactEVADPixL1.txt','SUSYBSMAnalysis/HSCP/data/CorrFactEVADPixL2.txt','SUSYBSMAnalysis/HSCP/data/CorrFactEVADPixL3.txt','SUSYBSMAnalysis/HSCP/data/CorrFactEVADPixL4.txt','SUSYBSMAnalysis/HSCP/data/CorrFactEVADPixR1.txt','SUSYBSMAnalysis/HSCP/data/CorrFactEVADPixR2.txt','SUSYBSMAnalysis/HSCP/data/GiTEMPLATE','MuonTimeOffset.txt']
config.JobType.pyCfgParams = ['GTAG=106X_dataRun2_v36', 'SAMPLE=isData', 'YEAR=EVAD', 'ERA=IDOSZAK']

config.section_('Data')
config.Data.inputDataset = 'MINTA'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 500
config.Data.outputDatasetTag = config.General.requestName
config.Data.outLFNDirBase = '/store/user/tvami/HSCP'
config.Data.ignoreLocality = True
config.Data.partialDataset = False 
config.Data.publication = False

config.section_('Site')
config.Site.whitelist = ['T2_DE_DESY','T2_CH_CERN','T2_IT_Bari','T1_IT_*','T2_US_*', 'T3_US_FNALLPC','T2_HU_Budapest','T2_FR_*', 'T2_UK_London_IC']
config.Site.storageSite = 'T2_HU_Budapest'
#config.Site.storageSite = 'T3_US_FNALLPC'
  '''

  with open("4crab_TemplateMTtape.py", "w") as text_file:
      text_file.write(TEMPLATE)

def task(i):
  os.system("cp 4crab_TemplateMTtape.py 4crab_toSubmit"+str(i.replace("/","_"))+".py")
  replaceVERZIO = "sed -i 's/VERZIO/"+codeVersion+"/g' 4crab_toSubmit"+str(i.replace("/","_"))+".py"
  os.system(replaceVERZIO)
  shortSampleName = i[1:(i.find('-'))].replace("/","_") if ("SingleMuon" in i) else  i[1:(i.find('TuneCP5'))].replace("/","_")
  replaceROVIDMINTA = "sed -i 's/ROVIDMINTA/"+shortSampleName+"/g' 4crab_toSubmit"+str(i.replace("/","_"))+".py"
  os.system(replaceROVIDMINTA)
  replaceMINTA = "sed -i 's/MINTA/"+i.replace("/","\/")+"/g' 4crab_toSubmit"+str(i.replace("/","_"))+".py"
  os.system(replaceMINTA)
  if ("2016" in i) :
    EVAD = "2016"
    replaceMASZK = "sed -i 's/#MASZK2016//g' 4crab_toSubmit"+str(i.replace("/","_"))+".py"
  if ("2017" in i) :
    EVAD = "2017"
    replaceMASZK = "sed -i 's/#MASZK2017//g' 4crab_toSubmit"+str(i.replace("/","_"))+".py"
  elif ("2018" in i) :
    EVAD = "2018"
    replaceMASZK = "sed -i 's/#MASZK2018//g' 4crab_toSubmit"+str(i.replace("/","_"))+".py"
  replaceEVAD  = "sed -i 's/EVAD/"+EVAD+"/g' 4crab_toSubmit"+str(i.replace("/","_"))+".py"
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
  replaceIDOSZAK  = "sed -i 's/IDOSZAK/"+IDOSZAK+"/g' 4crab_toSubmit"+str(i.replace("/","_"))+".py"
  os.system(replaceIDOSZAK)
  	
  replaceGiTEMPLATE = "sed -i 's/GiTEMPLATE/template_2018MC_v4.root/g' 4crab_toSubmit"+str(i.replace("/","_"))+".py"
  os.system(replaceGiTEMPLATE)

  os.system("crab submit -c 4crab_toSubmit"+str(i.replace("/","_"))+".py")
  os.system("mv 4crab_toSubmit"+str(i.replace("/","_"))+".py submittedConfigs/"+codeVersion+"_"+EVAD+IDOSZAK+".py")

for dataset in datasetList:
  t = Thread(target=task, args=(dataset,))
  t.start()


os.system("rm 4crab_TemplateMTtape.py")
os.system("rm *pyc")

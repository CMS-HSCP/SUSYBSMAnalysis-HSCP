import sys, os, time, re
import numpy as np
from threading import Thread
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog codeVersion")
(opt,args) = parser.parse_args()

datasetList = [
"/QCD_Pt-50To80_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM",
"/QCD_Pt-80To120_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM",
"/QCD_Pt-120To170_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM",
"/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM",
"/ZToMuMu_M-50To120_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/ZToMuMu_M-120To200_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/ZToMuMu_M-200To400_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/ZToMuMu_M-400To800_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/ZToMuMu_M-800To1400_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/ZToMuMu_M-1400To2300_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM",
"/ZToMuMu_M-2300To3500_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/ZToMuMu_M-3500To4500_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM",
"/ZToMuMu_M-4500To6000_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
"/ZToMuMu_M-6000ToInf_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM",
]

codeVersion = sys.argv[1]
#just the number, like 18p2

if not os.path.exists("submittedConfigs"): os.makedirs("submittedConfigs")

if not os.path.exists("4crab_Bkg_TemplateMT.py"):
  TEMPLATE = '''
from CRABClient.UserUtilities import config
config = config()

config.section_('General')
config.General.requestName = 'Analysis_2017_ROVIDMINTA_CodeVVERZIO_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HSCParticleProducerAnalyzer_2017_mc_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 3000
config.JobType.maxMemoryMB = 4000
config.JobType.inputFiles = ['SUSYBSMAnalysis/HSCP/data/template_2017MC_v5.root','SUSYBSMAnalysis/HSCP/data/MuonTimeOffset.txt']
config.JobType.pyCfgParams = ['GTAG=106X_mc2017_realistic_v10', 'SAMPLE=isBckg', 'YEAR=2017', 'ERA=MC']

config.section_('Data')
config.Data.inputDataset = 'MINTA'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 75
config.Data.outputDatasetTag = config.General.requestName
config.Data.outLFNDirBase = '/store/user/sawkim/HSCP'
config.Data.ignoreLocality = True
config.Data.partialDataset = True
config.Data.publication = False

config.section_('Site')
config.Site.whitelist = ['T2_DE_DESY','T2_CH_CERN','T2_IT_Bari','T1_IT_*','T2_US_*', 'T3_US_FNALLPC','T2_HU_Budapest','T2_FR_*', 'T2_UK_London_IC']
config.Site.storageSite = 'T2_HU_Budapest'
  '''

  with open("4crab_Bkg_TemplateMT.py", "w") as text_file:
      text_file.write(TEMPLATE)

def task(i):
  #print("Submit for sample "+i)
  os.system("cp 4crab_Bkg_TemplateMT.py 4crab_Bkg_toSubmit"+str(i.replace("/","_"))+".py")
  replaceVERZIO = "sed -i 's/VERZIO/"+codeVersion+"/g' 4crab_Bkg_toSubmit"+str(i.replace("/","_"))+".py"
  os.system(replaceVERZIO)
  shortSampleName = i[1:(i.find('TuneCP5'))-1]
  replaceROVIDMINTA = "sed -i 's/ROVIDMINTA/"+shortSampleName+"/g' 4crab_Bkg_toSubmit"+str(i.replace("/","_"))+".py"
  os.system(replaceROVIDMINTA)
  replaceMINTA = "sed -i 's/MINTA/"+i.replace("/","\/")+"/g' 4crab_Bkg_toSubmit"+str(i.replace("/","_"))+".py"
  os.system(replaceMINTA)
  os.system("crab submit -c 4crab_Bkg_toSubmit"+str(i.replace("/","_"))+".py")
  os.system("mv 4crab_Bkg_toSubmit"+str(i.replace("/","_"))+".py submittedConfigs/.")
  

for dataset in datasetList:
  t = Thread(target=task, args=(dataset,))
  t.start()

os.system("rm 4crab_Bkg_TemplateMT.py")
os.system("rm *pyc")

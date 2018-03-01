#!/usr/bin/env python

import urllib
import string
import os,sys,time
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor
import glob
import commands
import json
import collections

# 278816

PrimaryDatasetsToMatch = "/*/Run2016*-PromptReco-v*/AOD"
GlobalTag              = "80X_dataRun2_Prompt_v8"
StorageSite            = 'T2_BE_UCL'
AllLumisFile           = "FinalBatch.json"
StorageDir             = "/storage/data/cms/store/user/jozobec"
LumisPerJob            = 5

MuonTrigger1Mask       = 0
PFMetTriggerMask       = 0
L2MuMETTriggerMask     = 0

TransferDirectlyToStorage = True

LumiFiles = []
FilesToPush = []

def initProxy():
      print "You are going to run on a sample over grid using either CRAB or the AAA protocol, it is therefore needed to initialize your grid certificate"
      os.system('mkdir -p ~/x509_user_proxy; voms-proxy-init --voms cms -valid 192:00 --out ~/x509_user_proxy/x509_proxy')# shamelessly stolen from Loic

def getRunNum(lumiFileName):
   return lumiFileName.split('_')[1].split('.')[0]

def getRunList(lumiFileName):
   runs = []
   globalLumiFile = open(lumiFileName)
   runList = json.load(globalLumiFile).items()
   runList.sort()
   for entry in runList:
      runs.append("%s" % str(entry[0]))
   return runs

def splitLumiFileByRuns (lumiFileName):
   globalLumiFile = open(lumiFileName)
   runList = json.load(globalLumiFile).items()
   runList.sort()
   for run in runList:
      fileName = "Run_%s.json" % str(run[0])
      LumiFiles.append(fileName)
      f = open (fileName, 'w')
      f.write("{\"%s\": %s}" % (str(run[0]), str(run[1])))
      f.close()

def findDatasetsForRun (lumiFileName, dataset):
   runNumber = getRunNum (lumiFileName)
   return os.popen("das_client --limit 0 --query 'dataset dataset=%s run=%s' | grep -E '(/SingleMuon/|/DoubleMuon/)|(/MET/)'" % (dataset, runNumber)).read().split()

def createCrabConfigFile (lumiFileName, dataset):
   appendString = "MET"
   if dataset.find("SingleMu") != -1:
      appendString = "SingleMuon"
   elif dataset.find("DoubleMu") != -1:
      appendString = "DoubleMuon"

   runNumber    = getRunNum(lumiFileName)
   confFileName = "crab3_%s_%s.py" % (runNumber, appendString)
   FilesToPush.append(confFileName)
   f = open(confFileName, "w")
   f.write("from WMCore.Configuration import Configuration\n")
   f.write("config = Configuration()\n")
   f.write("config.section_('General')\n")
   f.write("config.General.requestName = 'Run%s_%s'\n" % (runNumber, appendString))
   f.write("config.General.transferOutputs = True\n")
   f.write("config.section_('JobType')\n")
   f.write("config.JobType.psetName = '../HSCParticleProducer_Data_cfg.py'\n")
   f.write("config.JobType.pluginName = 'Analysis'\n")
   f.write("config.JobType.pyCfgParams = ['globalTag=%s']\n" % GlobalTag)
   f.write("config.JobType.outputFiles = ['HSCP.root']\n")
   f.write("config.section_('Data')\n")
   f.write("config.Data.inputDataset = '%s'\n" % dataset)
   f.write("config.Data.publication = False\n")
   f.write("config.Data.unitsPerJob = %i\n" % LumisPerJob)
   f.write("config.Data.splitting = 'LumiBased'\n")
   f.write("config.Data.inputDBS = 'global'\n")
   f.write("config.Data.lumiMask = '%s'\n" % lumiFileName)
   f.write("config.section_('User')\n")
   f.write("config.section_('Site')\n")
   f.write("config.Site.storageSite = '%s'\n" % StorageSite)
   f.close()

def createToMergeList (paths):
   listOfFiles = []
   for path in paths:
      listOfFiles += os.popen("find %s -maxdepth 1 -name \"HSCP*.root\" -size +1024c" % path).read().split()
   listOfFiles.sort()
   toMergeList=""
   for i in range(0, len(listOfFiles)):
      toMergeList += "      \"file:"+listOfFiles[i]+"\",\n"
   LaunchOnCondor.Jobs_Inputs.append(toMergeList)

def createMergeConfigTemplate (templateName):
   f = open (templateName, 'w')
   f.write("import FWCore.ParameterSet.Config as cms\n")
   f.write("process = cms.Process(\"MergeHLT\")\n\n")
   f.write("process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )\n")
   f.write("process.load(\"FWCore.MessageService.MessageLogger_cfi\")\n")
   f.write("from SUSYBSMAnalysis.HSCP.HSCPVersion_cff import *\n\n")
   f.write("process.MessageLogger.cerr.FwkReport.reportEvery = 5000\n")
   f.write("process.source = cms.Source(\"PoolSource\",\n")
   f.write("   fileNames = cms.untracked.vstring( *(\n\n")
   f.write("XXX_INPUT_XXX\n")
   f.write("   ) )\n")
   f.write(")\n\n")
   f.write("process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )\n\n")
   f.write("process.HSCPHLTDuplicate = cms.EDFilter(\"HSCPHLTFilter\",\n")
   f.write("   RemoveDuplicates = cms.bool(True),\n")
   f.write("   TriggerProcess   = cms.string(\"HLT\"),\n")
   f.write("   MuonTrigger1Mask    = cms.int32(%i),\n" % MuonTrigger1Mask)
   f.write("   PFMetTriggerMask    = cms.int32(%i),\n" % PFMetTriggerMask)
   f.write("   L2MuMETTriggerMask  = cms.int32(%i),\n" % L2MuMETTriggerMask)
   f.write(")\n")
   f.write("process.DuplicateFilter = cms.Path(process.HSCPHLTDuplicate   )\n\n")
   f.write("process.Out = cms.OutputModule(\"PoolOutputModule\",\n")
   f.write("     outputCommands = cms.untracked.vstring(\n")
   f.write("	 \"keep *\"\n")
   f.write("    ),\n")
   f.write("    fileName = cms.untracked.string('XXX_SAVEPATH_XXX'),\n")
   f.write("    SelectEvents = cms.untracked.PSet(\n")
   f.write("       SelectEvents = cms.vstring('DuplicateFilter')\n")
   f.write("    ),\n")
   f.write(")\n")
   f.write("process.endPath = cms.EndPath(process.Out)\n")
   f.write("process.schedule = cms.Schedule(process.DuplicateFilter, process.endPath)\n")
   f.close()



if sys.argv[1] == '1':
   print "Splitting lumis ..."
   splitLumiFileByRuns (AllLumisFile)
   print "Done!\nWriting crab3 config files ...\n============================"
   for lumiFile in LumiFiles:
      matchedDatasets = findDatasetsForRun (lumiFile, PrimaryDatasetsToMatch)
      matchedDatasets.sort()
      print "----------------------------"
      for dataset in matchedDatasets:
         createCrabConfigFile (lumiFile, dataset)
         print "Run %s: %s ready!" % (getRunNum(lumiFile), dataset)

   print "Done!\nSubbmitting jobs on crab!\n============================"
   for confFile in FilesToPush:
      print "Submitting %s ..." % confFile
      os.system("crab submit -c %s" % confFile)
   


if sys.argv[1] == '2':
   print "Merging EDM files ..."
   FarmDirectory = "MERGECrab"
   MergeTemplateName = "Merge_Template_cfg.py"
   EndPath = ""

   if TransferDirectlyToStorage:
      print "Grid certificate is needed for the final lcg-cp command ..."
      initProxy()
      EndPath = "%s/HSCP2016" % StorageDir
   else:
      EndPath = "%s/%s/outputs" % (os.getcwd(), FarmDirectory)

   if not os.path.isdir(EndPath):
      os.system("mkdir -p %s" % EndPath)
   runs = getRunList(AllLumisFile)
   createMergeConfigTemplate(MergeTemplateName)
   LaunchOnCondor.SendCluster_Create(FarmDirectory, "HSCPEdmMerge")
   for run in runs:
      paths = ["%s/DoubleMuon/crab_Run%s_DoubleMuon/*/0000/" % (StorageDir, run),
               "%s/MET/crab_Run%s_MET/*/0000/" % (StorageDir, run),
               "%s/SingleMuon/crab_Run%s_SingleMuon/*/0000/" % (StorageDir, run)]
      createToMergeList(paths)

      LaunchOnCondor.Jobs_InitCmds   = ['export HOME=%s' % os.environ['HOME'], 'export X509_USER_PROXY=$HOME/x509_user_proxy/x509_proxy']
      LaunchOnCondor.Jobs_FinalCmds  = ['rm -f %s/Run2016_%s.root' % (EndPath, run)]


      if TransferDirectlyToStorage:
         LaunchOnCondor.Jobs_FinalCmds += ["lcg-cp -v -n 10 -D srmv2 -b file://${PWD}/Run2016_%s.root srm://ingrid-se02.cism.ucl.ac.be:8444/srm/managerv2\?SFN=%s/Run2016_%s.root && rm -f Run2016_%s.root" % (run, EndPath, run, run)] # if you do not use zsh, change '\?' to '?'
      else:
         LaunchOnCondor.Jobs_FinalCmds += ["mv Run2016_%s.root %s" % (run, EndPath)]

      LaunchOnCondor.Jobs_Queue = '8nh'
      LaunchOnCondor.SendCluster_Push(["CMSSW", MergeTemplateName, 'XXX_SAVEPATH_XXX', 'Run2016_%s.root' % run])
   LaunchOnCondor.SendCluster_Submit()
   os.system("rm -f %s" % MergeTemplateName)



if sys.argv[1] == '3':
   print "Validating the content and estimating the total integrated luminosity..."
   EndPath = ""
   if TransferDirectlyToStorage:
      EndPath = "%s/HSCP2016" % StorageDir
   else:
      FarmDirectory = "MERGECrab"
      EndPath = "%s/%s/outputs" % (os.getcwd(), FarmDirectory)
   runs = getRunList(AllLumisFile)
   f = open ("Analysis_Samples_tmp.txt", 'w')
   for run in runs:
      if os.path.isfile("%s/Run2016_%s.root" % (EndPath, run)):
         f.write("\"CMSSW_8_0\", 0, \"Data13TeV16\", \"%s/Run2016_%s.root\", \"Data\", \"NoPU\", 0, 0.0000000000E+00, 0, 0.000, 0.000, 0.000\n" % (EndPath, run))
   f.close()

#   listOfProcessedFiles = os.popen('find %s -maxdepth 1 -name "Run2016_*.root" | sort -V' % EndPath).read().split()
#   for line in listOfProcessedFiles:
#      f.write("\"CMSSW_8_0\", 0, \"Data13TeV16\", \"%s\", \"Data\", \"NoPU\", 0, 0.0000000000E+00, 0, 0.000, 0.000, 0.000\n" % line)
#   f.close()



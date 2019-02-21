#!/usr/bin/env python

import urllib
import string
import os,sys,time
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor  
import glob
import commands
import json
import collections # kind of map

goodLumis = {}
def LoadJson(jsonPath):
   jsonFile = open(jsonPath,'r')
   runList=json.load(jsonFile,encoding='utf-8').items()
   runList.sort()
   for run in runList :
      goodLumis[int(run[0])] = []
      for lumi in run[1] : goodLumis[int(run[0])].append(lumi)

def IsGoodRun(R):
   if(R in goodLumis): return True
   return False

def IsGoodLumi(R, L):
   for lumi in goodLumis[R]:
      if(L>=lumi[0] and L<=lumi[1]): return True
   return False


def IsFileWithGoodLumi(F):  
   #check if the file contains at least one good lumi section using DAS_CLIENT --> commented out because VERY SLOW!
   #print 'das_client.py --limit=0 --query "lumi file='+f+' |  grep lumi.run_number,lumi.number"'      
   #containsGoodLumi = False
   #for run in commands.getstatusoutput('das_client.py --limit=0 --query "run file='+f+'"')[1].replace('[','').replace(']','').split(','):
   #   if(not IsGoodRun(int(run))):continue
   #   for lumi in commands.getstatusoutput('das_client.py --limit=0 --query "lumi file='+f+'"')[1].replace('[','').replace(']','').split(','):
   #      if(IsGoodLumi(run, lumi)):return True

   #FASTER technique only based on run number and file name parsing
   run = int(F.split('/')[8])*1000+int(F.split('/')[9])
   if(IsGoodRun(run)):return True
   return False

   
def getChunksFromList(MyList, n):
  return [MyList[x:x+n] for x in range(0, len(MyList), n)]

JSON = 'Calib.json'
#JSON = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/DCSOnly/json_DCSONLY_Run2015B.txt'#/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/DCSOnly/json_DCSONLY.txt'   
LOCALTIER   = 'T2_CH_CERN'
#DATASETMASK = '/StreamExpress/Run2015B-SiStripCalMinBias-Express-v1/ALCARECO'
#DATASETMASK = '/MinimumBias/Run2015B-SiStripCalMinBias-PromptReco-v1/ALCARECO'
#DATASETMASK = '/StreamExpress/Run2015B-SiStripCalMinBias-Express-v1/ALCARECO'
#DATASETMASK = '/ZeroBias/Run2015B-PromptReco-v1/RECO'
DATASETMASK = ['/ZeroBias/Run2016F-SiStripCalMinBias-18Apr2017-v1/ALCARECO', '/ZeroBias/Run2016G-SiStripCalMinBias-18Apr2017-v1/ALCARECO']
EndPath     = "/storage/data/cms/store/user/jpriscia/dEdxCalib"
ISLOCAL     = False
TransferDirectlyToStorage = False
LoadJson(JSON)

def initProxy():
   if(not os.path.isfile(os.path.expanduser('~/x509_user_proxy/x509_proxy')) or ((time.time() - os.path.getmtime(os.path.expanduser('~/x509_user_proxy/x509_proxy')))>600)):
      print "You are going to run on a sample over grid using either CRAB or the AAA protocol, it is therefore needed to initialize your grid certificate"
      os.system('mkdir -p ~/x509_user_proxy; voms-proxy-init --voms cms -valid 192:00 --out ~/x509_user_proxy/x509_proxy')#all must be done in the same command to avoid environement problems.  Note that the first sourcing is only needed in Louvain



def filesFromDataset(dataset):
   ISLOCAL=False
   command_out = commands.getstatusoutput('das_client.py --limit=0 --query "site dataset='+dataset+' | grep site.name,site.dataset_fraction"')
   for site in command_out[1].split('\n'):
      if(LOCALTIER in site and '100.00%' in site): 
         ISLOCAL=True
         break

   Files = []
   command_out = commands.getstatusoutput('das_client.py --limit=0 --query "file dataset='+dataset+'"')
   for f in command_out[1].split():
      if(not IsFileWithGoodLumi(f)):continue
      if(ISLOCAL): Files += [f]
      else       : Files += ['root://cms-xrd-global.cern.ch/' + f]
   return Files

def filesFromDataset2(dataset):
   Files = []
   Runs  = ['278018', '278308', '279931', '280385']
   for run in Runs:
      output = os.popen('das_client --limit=0 --query \'file run=%s dataset=%s\'' % (run, dataset)).read().split('\n')
      for f in output:
         if len(f)<5: continue
         Files.append([run, 'root://cms-xrd-global.cern.ch//' + f])
   return Files


#get the list of sample to process from das and datasetmask query
print("Initialize your grid proxy in case you need to access remote samples\n")
initProxy()

#   command_out = commands.getstatusoutput('das_client --limit=0 --query "dataset='+DATASETMASK+'"')
datasetList = DATASETMASK

#get the list of samples to process from a local file
#datasetList= open('DatasetList','r')
JobName = "DEDXSKIMMER"
FarmDirectory = "FARM_EDM"
LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
#LaunchOnCondor.Jobs_Queue = '8nh'
LaunchOnCondor.subTool = 'condor'


if not TransferDirectlyToStorage:
   os.system("mkdir -p /eos/user/j/jpriscia/out");
else:
   os.system('mkdir -p %s/{278018,278308,279931,280385}' % EndPath)
for DATASET in datasetList :
   DATASET = DATASET.replace('\n','')
   FILELIST = filesFromDataset2(DATASET)
   print DATASET + " --> " + str(FILELIST)

   LaunchOnCondor.Jobs_InitCmds = []
#      if(not ISLOCAL):LaunchOnCondor.Jobs_InitCmds = ['export X509_USER_PROXY=~/x509_user_proxy/x509_proxy; voms-proxy-init --noregen;']
   if(not ISLOCAL):LaunchOnCondor.Jobs_InitCmds = ['export X509_USER_PROXY=~/x509_user_proxy/x509_proxy']

   for inFileList in getChunksFromList(FILELIST,1):
      os.system("cp dEdxSkimmer_Template_cfg.py dEdxSkimmer_cfg.py")
      f = open("dEdxSkimmer_cfg.py", "a")
      f.write("\n")
      f.write("process.Out.fileName = cms.untracked.string('dEdxSkim.root')\n")
      f.write("\n")
      for inFile in inFileList:
          f.write("process.source.fileNames.extend(['"+inFile[1]+"'])\n")
      f.write("\n")
      f.write("#import PhysicsTools.PythonAnalysis.LumiList as LumiList\n")
      f.write("#process.source.lumisToProcess = LumiList.LumiList(filename = '"+JSON+"').getVLuminosityBlockRange()")
      f.write("\n")

      if("/ALCARECO" in DATASET):
         f.write("\n")
         f.write("process.tracksForDeDx.src = cms.InputTag('ALCARECOSiStripCalMinBias') #for SiStripCalMinBias ALCARECO format\n")
         f.write("\n")
      f.close()   
      if not TransferDirectlyToStorage:
         LaunchOnCondor.Jobs_FinalCmds = ["cp dEdxSkim.root " + "/eos/user/j/jpriscia/out/dEdxSkim_%s_%i.root && rm dEdxSkim.root" % (inFile[0], LaunchOnCondor.Jobs_Count)]
      else:
         LaunchOnCondor.Jobs_FinalCmds = ["lcg-cp -v -n 10 -D srmv2 -b file://${PWD}/dEdxSkim.root srm://ingrid-se02.cism.ucl.ac.be:8444/srm/managerv2\?SFN=%s/%s/dEdxSkim_%s_%i.root && rm -f dEdxSkim.root" % (EndPath, inFile[0], inFile[0], LaunchOnCondor.Jobs_Count)] # if you do not use zsh, change '\?' to '?'
      LaunchOnCondor.SendCluster_Push  (["CMSSW", "dEdxSkimmer_cfg.py" ])
      os.system("rm -f dEdxSkimmer_cfg.py")

LaunchOnCondor.SendCluster_Submit()


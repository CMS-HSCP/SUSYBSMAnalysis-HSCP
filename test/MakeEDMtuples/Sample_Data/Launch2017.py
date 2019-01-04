
#!/usr/bin/env python 

import urllib
import string
import os,sys,time
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor
import glob
import commands
import json
import collections
from pdb import set_trace 
from FWCore.ParameterSet.VarParsing import VarParsing
import argparse
import itertools


parser = argparse.ArgumentParser()

parser.add_argument('STEP', type=int)
parser.add_argument('GTAG', type=str)
parser.add_argument('LUMITOPROCESS', type=str)
parser.add_argument('SAMPLE', type=str,choices=['isData','isSignal','isBkg'])
parser.add_argument('isSkimmedSample', type=bool)
parser.add_argument('--OUTPUTFILE', type=str,default='HSCP.root')
args = parser.parse_args()


#JSON = 'test_JSON.txt'
JSON = args.LUMITOPROCESS
OUTPUTFILE = args.OUTPUTFILE
LOCALTIER   = 'T2_BE_UCL'
#DATASETMASKS = ['/MET/Run2017*-17Nov2017-v1/AOD', '/SingleMuon/Run2017*-17Nov2017-v1/AOD']#, '/DoubleMuon/Run2017*-17Nov2017-v1/AOD']
ISLOCAL     = False
PrimaryDatasetsToMatch = "/*/Run2017*-17Nov2017-v1/AOD"


goodLumis = {}
def LoadJson(jsonPath):
   jsonFile = open(jsonPath,'r')
   runList=json.load(jsonFile,encoding='utf-8').items()
   runList.sort()
   for run in runList :
      goodLumis[int(run[0])] = []
      runList=[range(i[0], i[-1]+1) for i in run[1]]
      mergedList = list(itertools.chain.from_iterable(runList))
      goodLumis[int(run[0])]=[list(l) for l in zip(mergedList, mergedList)]
      #print [list(l) for l in zip(mergedList, mergedList)]
      #for lumi in run[1] : goodLumis[int(run[0])].append(lumi)

LumiFiles=[]
def splitLumiFileByRuns (goodLumis):
   # globalLumiFile = open(lumiFileName)
   #runList = json.load(globalLumiFile).items()
   #runList.sort()
   #for run in runList:
   #   goodLumis[int(run[0])] = []
   #   runList=[range(i[0], i[-1]+1) for i in run[1]]
   #   mergedList = list(itertools.chain.from_iterable(runList))
   #   goodLumis[int(run[0])].append([list(l) for l in zip(mergedList, mergedList)])
   #   print goodLumis[int(run[0])]
   #set_trace()
   for key,value in goodLumis.iteritems() :
      lumi_per_job = [value[i:i + 5] for i in range(0, len(value), 5)]
      for lumi_block in lumi_per_job:
         fileName = "Run_%s_%s.json" %(str(key),str(lumi_block[0][0]))
         LumiFiles.append(fileName)
         #lumi_block_formatted = [lumi_block[i:i + 1] for i in range(0, len(lumi_block), 1)]
         f = open (fileName, 'w')
         f.write("{\"%s\": %s}" % (str(key), str(lumi_block)))
         f.close()

def getRunNum(lumiFileName):
   return lumiFileName.split('_')[1].split('.')[0]
def getLumiNum(lumiFileName):
   return lumiFileName.split('_')[2].split('.')[0]

def findDatasetsForRun (lumiFileName, dataset):
   runNumber = getRunNum (lumiFileName)
   return os.popen("dasgoclient --limit 0 --query 'dataset dataset=%s run=%s' | grep -E '(/SingleMuon/)|(/MET/)'" % (dataset, runNumber)).read().split()

def filesFromDaset(lumi_dict,dataset):
    ISLOCAL=False
    command_out_1 = commands.getstatusoutput('dasgoclient --limit=0 --query "site dataset='+dataset+' | grep site.name,site.dataset_fraction"')
    for site in command_out_1[1].split('\n'):
        if(LOCALTIER in site and '100.00%' in site):
            ISLOCAL=True
            break
    command_out = commands.getstatusoutput("dasgoclient --limit 0 --query 'file dataset="+dataset+" run=%s'"% (runN))
   
    for f in command_out[1].split():
       value = eval(commands.getstatusoutput("dasgoclient --limit 0 --query 'lumi file="+f+"'")[1].replace(']\n[',','))
       if(ISLOCAL and LOCALTIER=='T2_CH_CERN'):  f = "root://eoscms//eos/cms"+f
       elif(ISLOCAL):    f = f
       elif 'T1_US_FNAL_Disk   100.00%' in command_out_1[1]: f = 'root://cmsxrootd.fnal.gov/' + f
       else:   f = 'root://cms-xrd-global.cern.ch/' + f
       print f
       lumi_dict[f]= value


def findFilesPerLumi(lumi_dict,lumiFile):
    jsonFile = open(lumiFile,'r')
    runList=json.load(jsonFile,encoding='utf-8').items()
    runList.sort()
    lumiSectInFile = list(itertools.chain.from_iterable(runList[0][1]))
    #set_trace()
    lumiSectInFile1, lumiSectInFile2 = zip(*runList[0][1])
    lumiSectInFile =list(lumiSectInFile1)
    for key, value in lumi_dict.iteritems():
       #print set(value), set(lumiSectInFile)
       if not set(lumiSectInFile).isdisjoint(set(value)):
          FilesByRun.append(key)
          #if(ISLOCAL and LOCALTIER=='T2_CH_CERN'):
          #   FilesByRun.append("root://eoscms//eos/cms"+key)
          #elif(ISLOCAL):     FilesByRun.append(f)
          #else:     FilesByRun.append('root://cms-xrd-global.cern.ch/' + key)

################################
##Proxy initialization##
def initProxy():
    print "You are going to run on a sample over grid using either CRAB or the AAA protocol, it is therefore needed to initialize your grid certificate"
    os.system('mkdir -p ~/x509_user_proxy; voms-proxy-init --voms cms -valid 192:00 --out ~/x509_user_proxy/x509_proxy')
    os.system('voms-proxy-init --voms cms -valid 192:00')

################################

LoadJson(JSON)
if args.STEP==0:
   os.system("rm -f *.d *.so *.pcm *.root *.pyc")
   os.system("rm -rdf FARM_TEMP out pictures")

if args.STEP==1:
    print 'proxy initialization..\n==============='
    initProxy()
    print 'splitting json by run....\n=============='
    splitLumiFileByRuns(goodLumis)
    
    os.system("mkdir -p out")
    JobName = "HSCPEdmProd"
    FarmDirectory = "FARM"
    LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
    LaunchOnCondor.Jobs_Queue = '8nh'
    LaunchOnCondor.Jobs_InitCmds = []
    oldRunNumber=0
    lumi_dict={}
    for lumiFile in LumiFiles:
        
      matchedDatasets = findDatasetsForRun(lumiFile, PrimaryDatasetsToMatch)
      matchedDatasets.sort()
      runN = getRunNum(lumiFile)
      lumiN = getLumiNum(lumiFile)
      print matchedDatasets

      
      if runN !=oldRunNumber:
         lumi_dict={}
         FilesByRun=[]
         print "new Run Number", runN
         for dataset in matchedDatasets:
            filesFromDaset(lumi_dict,dataset)
         findFilesPerLumi(lumi_dict,lumiFile)
         oldRunNumber=runN

      else:
         print "old Run Number", runN
         #print lumi_dict
         FilesByRun=[]
         findFilesPerLumi(lumi_dict,lumiFile)


      if(not ISLOCAL):LaunchOnCondor.Jobs_InitCmds = ['export HOME=%s' % os.environ['HOME'], 'export X509_USER_PROXY=~/x509_user_proxy/x509_proxy; voms-proxy-init --noregen;','set -o errexit']
      os.system("mkdir -p out/"+str(runN)+"_"+str(lumiN))
      LaunchOnCondor.Jobs_FinalCmds = ["mv %s %s/out/%s_%s/HSCP_%s_%s.root" % (OUTPUTFILE, os.getcwd(), runN,lumiN, runN,lumiN),"touch %s/out/%s_%s/HSCP_%s_%s.txt"% (os.getcwd(), runN,lumiN, runN,lumiN)]
      LaunchOnCondor.SendCluster_Push(
         ["CMSSW", ["HSCParticleProducer2017_cfg.py"]], 
         "GTAG=%s OUTPUTFILE=%s SAMPLE=%s isSkimmedSample=%s LUMITOPROCESS=%s inputFiles=%s" %(args.GTAG, OUTPUTFILE, args.SAMPLE, args.isSkimmedSample,os.getcwd()+'/'+lumiFile,','.join(FilesByRun)),
         index='%s_%s' % (runN,lumiN)
         )
      print '-------------'
      print lumiFile
      print FilesByRun
      print len(FilesByRun)
      print '-------------'
    LaunchOnCondor.SendCluster_Submit()

else: print "choose step 0 or 1"

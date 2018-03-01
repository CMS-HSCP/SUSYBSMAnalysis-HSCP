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
      for lumi in run[1] : goodLumis[int(run[0])].append(lumi)

LumiFiles=[]
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

def getRunNum(lumiFileName):
   return lumiFileName.split('_')[1].split('.')[0]

def findDatasetsForRun (lumiFileName, dataset):
   runNumber = getRunNum (lumiFileName)
   return os.popen("dasgoclient --limit 0 --query 'dataset dataset=%s run=%s' | grep -E '(/SingleMuon/|/DoubleMuon/)|(/MET/)'" % (dataset, runNumber)).read().split()


def filesFromDaset(dataset):
    ISLOCAL=False
    command_out = commands.getstatusoutput('dasgoclient --limit=0 --query "site dataset='+dataset+' | grep site.name,site.dataset_fraction"')
    for site in command_out[1].split('\n'):
        if(LOCALTIER in site and '100.00%' in site):
            ISLOCAL=True
            break
    command_out = commands.getstatusoutput("dasgoclient --limit 0 --query 'file dataset="+dataset+" run=%s'"% (runN))

    for f in command_out[1].split():
        if(ISLOCAL and LOCALTIER=='T2_CH_CERN'): FilesByRun.append("root://eoscms//eos/cms"+f)
        elif(ISLOCAL):                           FilesByRun.append(f)
        else:                                    FilesByRun.append('root://cms-xrd-global.cern.ch/' + f)
            

################################
##Proxy initialization##
def initProxy():
    print "You are going to run on a sample over grid using either CRAB or the AAA protocol, it is therefore needed to initialize your grid certificate"
    #os.system('mkdir -p ~/x509_user_proxy; voms-proxy-init --voms cms -valid 192:00 --out ~/x509_user_proxy/x509_proxy')
    os.system('voms-proxy-init --voms cms -valid 192:00')

################################

LoadJson(JSON)

if args.STEP==0:
   os.system("rm -f *.d *.so *.pcm *.root *.pyc")
   os.system("rm -rdf FARM_TEMP out pictures")

if args.STEP==1:
    print 'proxy initialization..\n==============='
    #initProxy()
    print 'splitting json by run....\n=============='
    splitLumiFileByRuns(JSON)
    
    os.system("mkdir -p out")
    JobName = "HSCPEdmProd"
    FarmDirectory = "FARM"
    LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
    LaunchOnCondor.Jobs_Queue = '8nh'
    LaunchOnCondor.Jobs_InitCmds = []
    for lumiFile in LumiFiles:
        
      matchedDatasets = findDatasetsForRun(lumiFile, PrimaryDatasetsToMatch)
      matchedDatasets.sort()
      runN = getRunNum(lumiFile)
      print matchedDatasets

      FilesByRun=[]
      for dataset in matchedDatasets:
          filesFromDaset(dataset)

      if(not ISLOCAL):LaunchOnCondor.Jobs_InitCmds = ['export HOME=%s' % os.environ['HOME'], 'export X509_USER_PROXY=~/x509_user_proxy/x509_proxy; voms-proxy-init --noregen;']
      os.system("mkdir -p out/"+str(runN));
      LaunchOnCondor.Jobs_FinalCmds = ["mv %s %s/out/%s/HSCP_%s.root" % (OUTPUTFILE, os.getcwd(), runN, runN)]
      LaunchOnCondor.SendCluster_Push(
         ["CMSSW", ["HSCParticleProducer2017_cfg.py"]], 
         "GTAG=%s OUTPUTFILE=%s SAMPLE=%s isSkimmedSample=%s LUMITOPROCESS=%s inputFiles=%s" %(args.GTAG, OUTPUTFILE, args.SAMPLE, args.isSkimmedSample,os.getcwd()+lumiFile,FilesByRun),
         index='%s_' % runN
         )
      print '-------------'
      print lumiFile
      print FilesByRun
      print '-------------'
    LaunchOnCondor.SendCluster_Submit()

else: print "choose step 0 or 1"

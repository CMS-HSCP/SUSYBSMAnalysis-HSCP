#!/usr/bin/env python

import urllib
import string
import os
import sys
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor  
import glob
import commands

LOCALTIER   = 'T2_BE_UCL'
#DATASETMASK = '/*HSCP*/RunIISpring15DR74-Asympt25ns*/AODSIM'
DATASETMASK = '/*HSCP*/RunIISpring15DR74-Asympt25ns*/GEN-SIM-RECO'
ISLOCAL     = False

def initProxy():
#   if(not os.path.isfile(os.path.expanduser('~/x509_user_proxy/x509_proxy')) or ((time.time() - os.path.getmtime(os.path.expanduser('~/x509_user_proxy/x509_proxy')))>600 and  int(commands.getstatusoutput('(export X509_USER_PROXY=~/x509_user_proxy/x509_proxy;voms-proxy-init --noregen;voms-proxy-info -all) | grep timeleft | tail -n 1')[1].split(':')[2])<8 )):
      print "You are going to run on a sample over grid using either CRAB or the AAA protocol, it is therefore needed to initialize your grid certificate"
      os.system('mkdir -p ~/x509_user_proxy; voms-proxy-init --voms cms -valid 192:00 --out ~/x509_user_proxy/x509_proxy')#all must be done in the same command to avoid environement problems.  Note that the first sourcing is only needed in Louvain


def nameFromDataset(line):
   Name =''
   if('gluino'   in line): Name+= 'Gluino'
   if('stop'     in line): Name+= 'Stop'
   if('ppstau'   in line): Name+= 'PPStau'
   if('gmstau'   in line): Name+= 'GMStau'
   if('mchamp'   in line): Name+= 'DY'

   if('7TeV'  in line): Name+= '_7TeV'
   if('8TeV'  in line): Name+= '_8TeV'
   if('13TeV' in line): Name+= '_13TeV'
   if('14TeV' in line): Name+= '_14TeV'

   MassStr = '_M' + (line.split('_M-')[1]).split('_')[0]
   Name+=MassStr

   if('onlyneutral'  in line): Name+= 'N'

   QStr = ''
   if('mchamp' in line): QStr = '_Q' + (line.split('mchamp')[1]).split('_')[0]
   Name+=QStr

   return Name

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
      if(ISLOCAL): Files += [f]
      else       : Files += ['root://cms-xrd-global.cern.ch/' + f]
   return Files
  



#get the list of sample to process from das and datasetmask query
print("Initialize your grid proxy in case you need to access remote samples\n")
#initProxy()

command_out = commands.getstatusoutput('das_client.py --limit=0 --query "dataset='+DATASETMASK+'"')
datasetList = command_out[1].split()

#get the list of samples to process from a local file
#datasetList= open('DatasetList','r')
for DATASET in datasetList :
   DATASET = DATASET.replace('\n','')
   NAME = nameFromDataset(DATASET)
   FILELIST = filesFromDataset(DATASET)
   print DATASET + " --> " + NAME + " --> " + str(FILELIST)

   JobName = NAME
   FarmDirectory = "FARM_EDM"
   LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
   LaunchOnCondor.Jobs_InitCmds = []
   if(not ISLOCAL):LaunchOnCondor.Jobs_InitCmds = ['export X509_USER_PROXY=~/x509_user_proxy/x509_proxy; voms-proxy-init --noregen;']

   f = open("HSCPEDM_cfg.py", "w")
   f.write("import sys, os\n")
   f.write("import FWCore.ParameterSet.Config as cms\n")
   f.write("\n")
   if('HSCP' in DATASET):
         f.write("isSignal = True\n")
         f.write("isBckg = False\n")
   else:
         f.write("isSignal = False\n")
         f.write("isBckg = True\n")
   f.write("isData = False\n")
   f.write("isSkimmedSample = False\n")
   f.write("GTAG = 'MCRUN2_74_V9'\n")
   f.write("OUTPUTFILE = 'XXX_OUTPUT_XXX.root'\n")
   f.write("InputFileList = cms.untracked.vstring()\n")
   f.write("\n")
   for inFile in FILELIST:
      f.write("InputFileList.extend(['"+inFile+"'])\n")
   f.write("\n")
   f.write("#main EDM tuple cfg that depends on the above parameters\n")
   f.write("execfile( os.path.expandvars('${CMSSW_BASE}/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/HSCParticleProducer_cfg.py') )\n")
   f.close()

   LaunchOnCondor.SendCluster_Push  (["CMSSW", "HSCPEDM_cfg.py" ])
   LaunchOnCondor.SendCluster_Submit()

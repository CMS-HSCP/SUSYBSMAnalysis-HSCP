#!/usr/bin/env python
import urllib
import string
import os
import sys
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor
import glob
import time

LaunchOnCondor.Jobs_InitCmds       = ['ulimit -c 0;']  #disable production of core dump in case of job crash

UseRemoteSamples                   = True
RemoteStorageDir                   = '/storage/data/cms/store/user/jozobec/HSCP2016/'
#RemoteStorageDir                   = '/store/group/phys_exotica/hscp/'


def initProxy():
   if(not os.path.isfile(os.path.expanduser('~/private/x509_proxy')) or ((time.time() - os.path.getmtime(os.path.expanduser('~/private/x509_proxy')))>600)):
      print "You are going to run on a sample over grid using either CRAB or the AAA protocol, it is therefore needed to initialize your grid certificate"
      os.system('mkdir -p ~/private; voms-proxy-init --voms cms -valid 192:00 --out ~/private/x509_proxy')#all must be done in the same command to avoid environement problems.  Note that the first sourcing is only needed in Louvain


if sys.argv[1]=='1':
   if UseRemoteSamples:
      initProxy()
   print("compile the Stability code")
   os.system("sh " + os.getcwd() + "/StabilityCheck.sh ") #just compile

   print 'STABILITY'
   FarmDirectory = "FARM"
   JobName = "HSCPStability"
   LaunchOnCondor.Jobs_RunHere = 0
   LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
   LaunchOnCondor.subTool = 'condor'

   #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/StabilityCheck.C", '"ANALYSE"'])
   #LaunchOnCondor.SendCluster_Push(["BASH", "sh " + os.getcwd()+"/StabilityCheck.sh " + os.getcwd()+"/pictures"])


   #NJobs = 500
   #for Job in range(0,NJobs) :
   #      LaunchOnCondor.SendCluster_Push(["BASH", "sh " + os.getcwd()+"/StabilityCheck.sh " + os.getcwd()+"/pictures " + str(Job) +" " + str(NJobs)])
   #LaunchOnCondor.SendCluster_Submit()

   cwd = '%s/src/SUSYBSMAnalysis/HSCP/test/UsefulScripts/StabilityCheck' % os.environ['CMSSW_BASE']
   f= open('%s/../../AnalysisCode/Analysis_Samples.txt' % cwd,'r')
#   f= open('Analysis_Samples_tmp.txt','r')
   index = -1
   for line in f :
      index+=1           
      if(line.startswith('#')):continue
      vals=line.split(',')
      if(int(vals[1])==2):continue
      if int(vals[1]) == 0 and str(vals[2]).find("Data13TeV16") == -1: continue
      if(UseRemoteSamples and int(vals[1])==0 and vals[3].find('2016')):
         LaunchOnCondor.Jobs_InitCmds = ['ulimit -c 0',
                                         'export HOME=%s' % os.environ['HOME'],
                                         'export X509_USER_PROXY=%s/x509_user_proxy/x509_proxy; voms-proxy-init --noregen;' % os.environ['HOME'],
                                         'export REMOTESTORAGEPATH='+RemoteStorageDir.replace('/storage/data/cms/store/', '/store/')]
      else: LaunchOnCondor.Jobs_InitCmds = ['ulimit -c 0']
      LaunchOnCondor.SendCluster_Push(["BASH", "sh " + os.getcwd()+"/StabilityCheck.sh " + os.getcwd()+"/pictures " + str(index) +" " + str(1)])
   f.close()
   LaunchOnCondor.SendCluster_Submit()

elif sys.argv[1]=='2':
   os.system('find pictures/Histos_*.root  -type f -size +1024c | xargs hadd -f  pictures/Histos.root')

if sys.argv[1]=='3':
   os.system('sh MakePlot.sh')


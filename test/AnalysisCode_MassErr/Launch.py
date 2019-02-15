#!/usr/bin/env python

import urllib
import string
import os
import sys
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor 
import glob
import commands
import time

LaunchOnCondor.Jobs_InitCmds       = ['ulimit -c 0;']  #disable production of core dump in case of job crash
LaunchOnCondor.subTool = 'condor'
#LaunchOnCondor.Jobs_Queue = '8nh'

UseRemoteSamples          = True 
RemoteStorageDir          = '/storage/data/cms/store/user/jozobec/HSCP2016/'
RemoteServer              = 'cms-xrd-global.cern.ch'
#RemoteStorageDir          = '/store/group/phys_exotica/hscp/'

#the vector below contains the "TypeMode" of the analyses that should be run
AnalysesToRun = [0,2]#,4]#,3,5]

CMSSW_74X_VERSION = '7_4_16'
base80X = os.environ['CMSSW_BASE']
base74X = ""
dirs = base80X.split('/')
for i in range(0, len(dirs)-1):
   base74X+=dirs[i]+'/'
base74X+='CMSSW_'+CMSSW_74X_VERSION

if len(sys.argv)==1:       
        print "Please pass in argument a number between 1 and 5"
        print "  1  - Submit the Core of the (TkOnly+TkTOF) Analysis                 --> submitting %d x #SignalPoints jobs" % len(AnalysesToRun)
        print "  2  - Merge all output files and estimate backgrounds                --> submitting %d                 jobs"  % len(AnalysesToRun)
        print "  3  - Run the control plot macro                                     --> interactive"
        print "  4o - Run the Optimization macro based on best Exp Limit             --> submitting %d x #SignalPoints jobs (OPTIONAL)"  % len(AnalysesToRun)
        print "  4  - compute the limits from the cuts set in Analysis_Cuts.txt file --> submitting %d x #SignalPoints jobs (cut must edited by hand in Analysis_Cuts.txt)"  % len(AnalysesToRun)
        print "  5  - Run the exclusion plot macro                                   --> interactive"
        sys.exit()


CMSSW_VERSION = os.getenv('CMSSW_VERSION','CMSSW_VERSION')
if CMSSW_VERSION == 'CMSSW_VERSION':
  print 'please setup your CMSSW environement'
  sys.exit(0)

if("CMSSW_8" in CMSSW_VERSION): CMSSW_VERSION = CMSSW_VERSION + " CMSSW_7_4"

def skipSamples(type, name):
   if(name.find("AMSB")!=-1 and type!=0):return True; #only consider AMSB in TkOnly analysis

   if(type==3):
      if(name.find("Gluino")==-1 and name.find("Stop")==-1 and name.find("Stau")==-1 and name.find("o3")==-1):return True;
   elif(type==4):
       return False;
#      if(name.find("DY")==-1 or name.find("o3")>=0):return True;
   elif(type==5):
      if(name.find("DY")==-1 or (name.find("1o3")==-1 and name.find("2o3")==-1 and name.find("Q1")==-1)):return True;
   
   return False

def initProxy():
   if(not os.path.isfile(os.path.expanduser('~/x509_user_proxy/x509_proxy')) or ((time.time() - os.path.getmtime(os.path.expanduser('~/x509_user_proxy/x509_proxy')))>600)):
      print "You are going to run on a sample over grid using either CRAB or the AAA protocol, it is therefore needed to initialize your grid certificate"
      os.system('mkdir -p ~/x509_user_proxy; voms-proxy-init --voms cms -valid 192:00 --out ~/x509_user_proxy/x509_proxy')#all must be done in the same command to avoid environement problems.  Note that the first sourcing is only needed in Louvain


if sys.argv[1]=='1':	
        if UseRemoteSamples:
           initProxy()
        print 'ANALYSIS'
        FarmDirectory = "FARM"
        JobName = "HscpAnalysis"
        LaunchOnCondor.subTool = 'condor'
        LaunchOnCondor.Jobs_RunHere = 0
        LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
        f= open('Analysis_Samples.txt','r')
        index = -1
        for line in f :
           index+=1           
           if(line.startswith('#')):continue
           vals=line.split(',')
           if((vals[0].replace('"','')) in CMSSW_VERSION):
              for Type in AnalysesToRun:
                 #LaunchOnCondor.Jobs_FinalCmds = ['mv *.root %s/src/SUSYBSMAnalysis/HSCP/test/AnalysisCode/Results/Type%i/' % (os.environ['CMSSW_BASE'], Type)]
                 LaunchOnCondor.Jobs_FinalCmds = ['cp -r Results %s/src/SUSYBSMAnalysis/HSCP/test/%s/ && rm -rf Results' % (os.environ['CMSSW_BASE'], os.path.basename(os.getcwd()) if os.getcwd().find('AnalysisCode') > 0 else 'AnalysisCode')]
                 if(UseRemoteSamples):
                    LaunchOnCondor.Jobs_InitCmds = ['ulimit -c 0', 'export HOME=%s' % os.environ['HOME'], 'export X509_USER_PROXY=$HOME/x509_user_proxy/x509_proxy; voms-proxy-init --noregen;', 'export REMOTESTORAGESERVER='+RemoteServer, 'export REMOTESTORAGEPATH='+RemoteStorageDir.replace('/storage/data/cms/store/', '//store/')]
                 else: LaunchOnCondor.Jobs_InitCmds = ['ulimit -c 0']
                 if(int(vals[1])>=2 and skipSamples(Type, vals[2])==True):continue
                 if vals[2].find("13TeV16") < 0: continue # skip everything that is not 2016 -- namely the 2015 samples
#                 if(int(vals[1])!=1):continue #Skip everything except for background MC
                 if(int(vals[1])>0):continue #Skip all data
                 LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step1_EventLoop.C", '"ANALYSE_'+str(index)+'_to_'+str(index)+'"'  , Type, vals[2].rstrip() ])
        f.close()
        LaunchOnCondor.SendCluster_Submit()

elif sys.argv[1]=='2':
        print 'MERGING FILE AND PREDICTING BACKGROUNDS'  
        FarmDirectory = "FARM"
        JobName = "HscpPred"
        LaunchOnCondor.Jobs_RunHere = 1
        LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
        for Type in AnalysesToRun:
           Path = "Results/Type"+str(Type)+"/"
           os.system('rm -f ' + Path + 'Histos.root')
           #os.system('hadd -f ' + Path + 'Histos.root ' + Path + '*.root')           
           smallFiles = commands.getstatusoutput('find ' + Path + 'Histos_*.root  -type f -size -1024c -exec ls -lSh {} +')[1]
           if(len(smallFiles)>1):
              print("Small files have been found, these are generally due to either crashed jobs, or to still running jobs.\nThe following files will NOT be hadd:\n" + smallFiles + "\n\n")           
           os.system('find ' + Path + 'Histos_*.root  -type f -size +1024c | xargs hadd -n 50 -f ' + Path + 'Histos.root ')
           LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step2_BackgroundPrediction.C", '"'+Path+'"'])
        LaunchOnCondor.SendCluster_Submit()

elif sys.argv[1]=='3':
        print 'PLOTTING'
	os.system('root Analysis_Step3_MakePlots.C++ -l -b -q')

elif sys.argv[1]=='4o':
        print 'OPTIMIZATION & LIMIT COMPUTATION'
        FarmDirectory = "FARM"
        JobName = "HscpLimits"
        LaunchOnCondor.Jobs_RunHere = 1
        LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)

        f= open('Analysis_Samples.txt','r')
        for line in f :
           if(line.startswith('#')):continue
           vals=line.split(',')
           if(int(vals[1])!=2):continue
           for Type in AnalysesToRun:
              if(int(vals[1])>=2 and skipSamples(Type, vals[2])==True):continue
#              if(vals[2].find("8TeV")<0):continue
              Path = "Results/Type"+str(Type)+"/"
              LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step4_LimitComputation.C", '"OPTIMIZE"', '"'+Path+'"', vals[2] ])
        f.close()
        LaunchOnCondor.SendCluster_Submit()

elif sys.argv[1]=='4':
        print 'LIMIT COMPUTATION (ONLY)'
        FarmDirectory = "FARM"
        JobName = "HscpLimits"
        LaunchOnCondor.Jobs_Queue   = '8nh'
        LaunchOnCondor.Jobs_RunHere = 1
        LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)

        if not os.path.islink('%s/src/HiggsAnalysis' % base80X):
           if not os.path.isdir('%s/src/HiggsAnalysis' % base74X):
              print 'Error: %s/src/HiggsAnalysis does not exist' % base74X
              print 'Make sure you set up the correct CMSSW_7_4_X version and path in base74X variable.'
              print 'Exiting.'
              sys.exit()
           if not os.path.islink('%s/src/HiggsAnalysis' % base80X):
              print 'Creating the symlink to HiggsAnalysis/CombinedLimit tool from 74X'
              print 'cd %s/src && ln -s %s/src/HiggsAnalysis' % (base80X, base74X)
              os.system('cd %s/src && ln -s %s/src/HiggsAnalysis' % (base80X, base74X))

        LaunchOnCondor.Jobs_InitCmds += ['cd %s/src' % base74X, 'eval `scramv1 runtime -sh`', 'cd -']

        f = open('Analysis_Samples.txt','r')
        for line in f :
           if(line.startswith('#')):continue
           vals=line.split(',')
           if(int(vals[1])<2):continue
           if vals[2].find("13TeV16") == -1: continue # remove if you want to run on 2015 samples also
           if vals[2].find("13TeV16G") > 0: continue # only run on the 13TeV16 pre-G for combine
           for Type in AnalysesToRun:            
              if(int(vals[1])>=2 and skipSamples(Type, vals[2])==True):continue     
              skip = False

#              #skip 8TeV samples that have already been processed together with the  7TeV (since for each sample we do 7TeV+8TeV+Comb)
#              if(vals[2].find("8TeV")>=0):
#                  key = vals[2]
#                  key = key.replace("8TeV","7TeV")
#                  f2= open('Analysis_Samples.txt','r')
#                  for line2 in f2 :
#                     vals2=line2.split(',')
#                     if(vals2[1]==vals[1] and vals2[2] == key): skip = True;
#                  if(skip==True): continue;
#                  f2.close()
#              #print vals[2] + "   " + str(skip)

              Path = "Results/Type"+str(Type)+"/"
#              LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step4_LimitComputation.C", '"COMPUTELIMIT2016"' if vals[2].find("13TeV16G") == -1 else '"COMPUTELIMIT2016G"', '"'+Path+'"', vals[2] ]) # compute separately 2016 preG and postG

              LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step4_LimitComputation.C", '"COMBINE2016"', '"'+Path+'"', vals[2] ]) #combine 2016 preG, postG, and combined 2016 in the same job
        f.close()
        LaunchOnCondor.SendCluster_Submit()

elif sys.argv[1]=='5':
        print 'EXCLUSION'
        os.system('sh Analysis_Step4_LimitComputation.sh')
else:
        print 'Unknown case: use an other argument or no argument to get help'




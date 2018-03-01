#!/usr/bin/env python

# Original Author:  Loic Quertenmont


import urllib
import string
import os,sys,time
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor
import glob
import commands
import json
import collections # kind of map

#script parameters #feel free to edit those
#JSON = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt'
#JSON = 'FinalBatch.json'
JSON = '23Sept2016ReReco_Collisions16.json'
#JSON = os.getcwd() + '/Json.txt'
#/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt'   
LOCALTIER   = 'T2_CH_CERN'
DATASETMASK = '/DoubleMuon/Run2016*-PromptReco-v*/AOD'
ISLOCAL     = False #automatically assigned
STORAGEDIR = 'jzobec@ingrid-ui1.cism.ucl.ac.be:/nfs/scratch/fynu/jzobec/Run2Analysis/CMSSW_8_0_8_patch1/src/SUSYBSMAnalysis/HSCP/test/UsefulScripts/MuonTimingStudy/.'  #scp final dir


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

def GetRunFromFile(F):
   return int(F.split('/')[8])*1000+int(F.split('/')[9])   

   
def getChunksFromList(MyList, n):
  return [MyList[x:x+n] for x in range(0, len(MyList), n)]


def initProxy():
#   if(not os.path.isfile(os.path.expanduser('~/x509_user_proxy/x509_proxy')) or ((time.time() - os.path.getmtime(os.path.expanduser('~/x509_user_proxy/x509_proxy')))>600 and  int(commands.getstatusoutput('(export X509_USER_PROXY=~/x509_user_proxy/x509_proxy;voms-proxy-init --noregen;voms-proxy-info -all) | grep timeleft | tail -n 1')[1].split(':')[2])<8 )):
      print "You are going to run on a sample over grid using either CRAB or the AAA protocol, it is therefore needed to initialize your grid certificate"
      os.system('mkdir -p ~/x509_user_proxy; voms-proxy-init --voms cms -valid 192:00 --out ~/x509_user_proxy/x509_proxy')#all must be done in the same command to avoid environement problems.  Note that the first sourcing is only needed in Louvain


def filesFromDataset(dataset):
   ISLOCAL=False
#   command_out = commands.getstatusoutput('das_client --limit=0 --query "site dataset='+dataset+' | grep site.name,site.dataset_fraction"')
#   for site in command_out[1].split('\n'):
#      if(LOCALTIER in site and '100.00%' in site): 
#         ISLOCAL=True
#         break

   Files = {}
   command_out = commands.getstatusoutput('das_client --limit=0 --query "file dataset='+dataset+'"')
   for f in command_out[1].split():
      run = GetRunFromFile(f)
      if(not IsGoodRun(run)):continue
      if(run not in Files):  Files[run] = [] #make sure that we have a collection for this run
      
#      if(ISLOCAL and LOCALTIER=='T2_CH_CERN'): Files[run] += ["root://eoscms//eos/cms"+f]
#      elif(ISLOCAL):                           Files[run] += [f]
#      else       :                             Files[run] += ['root://cms-xrd-global.cern.ch/' + f]
      Files[run] += ['root://cms-xrd-global.cern.ch//' + f]
   return Files
  

if len(sys.argv)==1:
        print "Please pass in argument a number between 0 and 3"
        print "  0  - cleanup the workspace (all results are erased)        --> interactive processing" 
        print "  1  - Run MuonTimingStudy on RECO files from a DAS dataset  --> submitting 1job per file"
        print "  2  - Hadd root files containing the histograms (1file/run) --> interactive processing" 
        print "  3  - run the plotter on the hadded root files              --> interactive processing" 
        sys.exit()



LoadJson(JSON)
if sys.argv[1]=='0':
   os.system("rm -f *.d *.so *.pcm *.root *.pyc") 
   os.system("rm -rdf FARM_TEMP out pictures") 



if sys.argv[1]=='1':

   #get the list of sample to process from das and datasetmask query
   print("Initialize your grid proxy in case you need to access remote samples\n")
   initProxy()


   print("compile the MuonTimingStudy code")
   os.system("sh " + os.getcwd() + "/MuonTimingStudy.sh ") #just compile


   command_out = commands.getstatusoutput('das_client --limit=0 --query "dataset='+DATASETMASK+'"')
   print 'das_client --limit=0 --query "dataset='+DATASETMASK+'"'
   print command_out
   datasetList = command_out[1].split()

   #get the list of samples to process from a local file
   #datasetList= open('DatasetList','r')
   JobName = "CSCTimeStudy"
   FarmDirectory = "FARM_TEMP"
   LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
   LaunchOnCondor.Jobs_Queue = '8nh'

   os.system("mkdir -p out");
   for DATASET in datasetList :
      DATASET = DATASET.replace('\n','')
      FILELIST = filesFromDataset(DATASET)
      LaunchOnCondor.Jobs_InitCmds = []
      if(not ISLOCAL):LaunchOnCondor.Jobs_InitCmds = ['export X509_USER_PROXY=~/x509_user_proxy/x509_proxy; voms-proxy-init --noregen;']

      print DATASET + " : " 
      for RUN in FILELIST:
         print str(RUN) + " --> %i files to process"  % len(FILELIST[RUN])
         INDEX=0
         for inFileList in getChunksFromList(FILELIST[RUN],5):
            InputListCSV = ''
            for inFile in inFileList:
               InputListCSV+= inFile + ','
            InputListCSV = InputListCSV[:-1] #remove the last duplicated comma
            LaunchOnCondor.SendCluster_Push  (["BASH", "sh " + os.getcwd() + "/MuonTimingStudy.sh " + InputListCSV + " out.root " + str(RUN)+"; mv out.root " + os.getcwd() + "/out/Histos_%i_%i.root" % (RUN,INDEX)])
            INDEX+=1

   LaunchOnCondor.SendCluster_Submit()



if sys.argv[1]=='2':
   os.system('find out/Histos_*.root  -type f -size +1024c | xargs hadd -f  Histos.root')

if sys.argv[1]=='3':
   os.system('sh MakePlot.sh')

if sys.argv[1]=='4':
    counter=0
    indexfile=open("pictures/index.html",'w')
    PictureList = glob.glob("pictures/*.png")
    for Plot in PictureList:
        Plot = Plot.replace('pictures/', '')
	if(counter%3==0): indexfile.write('<TR> <TD align=center> <a href="'+Plot+'"><img src="'+Plot+'"hspace=5 vspace=5 border=0 style="width: 80%" ALT="'+Plot+'"></a><br> '+Plot+' </TD>\n')
        if(counter%3==1): indexfile.write('     <TD align=center> <a href="'+Plot+'"><img src="'+Plot+'"hspace=5 vspace=5 border=0 style="width: 80%" ALT="'+Plot+'"></a><br> '+Plot+' </TD>\n')
        if(counter%3==2): indexfile.write('     <TD align=center> <a href="'+Plot+'"><img src="'+Plot+'"hspace=5 vspace=5 border=0 style="width: 80%" ALT="'+Plot+'"></a><br> '+Plot+' </TD> </TR>\n')
        counter+=1
    indexfile.close()

if sys.argv[1]=='5':
    os.system('scp pictures/* ' + STORAGEDIR)


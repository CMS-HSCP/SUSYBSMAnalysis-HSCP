#!/usr/bin/env python

import urllib
import string
import os,sys,time
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor  
import glob
import commands
import json
import collections # kind of map
from pdb import set_trace



goodLumis = {}
def LoadJson(jsonPath):
   jsonFile = open(jsonPath,'r')
   runList=json.load(jsonFile,encoding='utf-8').items()
   runList.sort()
   for run in runList :
      goodLumis[int(run[0])] = []
      for lumi in run[1] : goodLumis[int(run[0])].append(lumi)

def getChunksFromList(MyList, n):
  return [MyList[x:x+n] for x in range(0, len(MyList), n)]

def filesFromDataset(dataset,nof=-1,start=0):
   ISLOCAL=False
   command_out = commands.getstatusoutput('das_client --limit=0 --query "site dataset='+dataset+' | grep site.name,site.dataset_fraction"')
   for site in command_out[1].split('\n'):
      if(LOCALTIER in site and '100.00%' in site): 
         ISLOCAL=True
         break

   NumberOfFiles = 0
   counter=0
   maxNofFiles = nof
   Files = []
   command_out = commands.getstatusoutput('das_client --limit=0 --query "file dataset='+dataset+'"')
   for f in command_out[1].split():
#      if(not IsFileWithGoodLumi(f)):continue
      counter += 1
      if counter<start: continue
      if(ISLOCAL): Files += [f]
      else       : Files += ['root://cms-xrd-global.cern.ch/' + f]
      NumberOfFiles += 1
      if maxNofFiles>0 and NumberOfFiles > maxNofFiles: break
   return Files

def sublistFilesFromList(listFiles,segm=5):
    lR = []
    subl = []
    counter=0
    for inList in listFiles:
        if counter<segm:
            subl += [inList]
        if counter>=segm:
            counter=0
            lR += [subl]
            subl = []
        counter += 1
    return lR


def initProxy():
   os.system('rm -rf ~/x509_user_proxy/x509_proxy')
   #os.system('rm -rf ~/private/x509_proxy')
   if(not os.path.isfile(os.path.expanduser('~/private/x509_proxy')) or ((time.time() - os.path.getmtime(os.path.expanduser('~/private/x509_proxy')))>600)):
      print "You are going to run on a sample over grid using either CRAB or the AAA protocol, it is therefore needed to initialize your grid certificate"
      os.system('mkdir -p ~/x509_user_proxy; voms-proxy-init --voms cms -rfc -valid 192:00 --out ~/private/x509_proxy')#all must be done in the same command to avoid environement problems.  Note that the first sourcing is only needed in Louvain


if len(sys.argv)==1:
        print "Please pass in argument a number between 1 and 3"
        print "  1  - Run dEdxStudy on RECO, AOD, or dEdxSKIM files         --> submitting 1job per file"
        print "  2  - Hadd root files containing the histograms             --> interactive processing" 
        print "  3  - run the plotter on the hadded root files              --> interactive processing" 
        sys.exit()

JSON = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
LOCALTIER = 'T2_CH_CERN'



datasetList = [
##  ["Run251252", "/storage/data/cms/store/user/querten/AAA/HSCP/dedxSkim/Run251252/"],
##  ["Run257395", "/storage/data/cms/store/user/jozobec/out/257395/"],
##  ["Run257396", "/storage/data/cms/store/user/jozobec/out/257396/"],
#  ["Run257490", "/storage/data/cms/store/user/jozobec/out/257490/"],
##  ["Run257599", "/storage/data/cms/store/user/jozobec/out/257599/"],
#  ["Run257805", "/storage/data/cms/store/user/jozobec/out/257805/"],
##  ["Run257816", "/storage/data/cms/store/user/jozobec/out/257816/"],
##  ["Run257819", "/storage/data/cms/store/user/jozobec/out/257819/"],
#  ["Run257823", "/storage/data/cms/store/user/jozobec/out/257823/"],
##  ["Run258443", "/storage/data/cms/store/user/jozobec/out/258443/"],
##  ["Run258446", "/storage/data/cms/store/user/jozobec/out/258446/"],
##  ["Run258655", "/storage/data/cms/store/user/jozobec/out/258655/"],
##  ["Run258694", "/storage/data/cms/store/user/jozobec/out/258694/"],
##  ["Run258705", "/storage/data/cms/store/user/jozobec/out/258705/"],
##  ["Run258741", "/storage/data/cms/store/user/jozobec/out/258741/"],
##  ["Run258750", "/storage/data/cms/store/user/jozobec/out/258750/"],
#  ["Run278018", "/storage/data/cms/store/user/jozobec/dEdxCalib/278018/","listFiles_278018.txt"],
#  ["Run278308", "/storage/data/cms/store/user/jozobec/dEdxCalib/278308/","listFiles_278308.txt"],
#  ["Run279931", "/storage/data/cms/store/user/jozobec/dEdxCalib/279931/","listFiles_279931.txt"],
#  ["RunPostG", "/eos/user/j/jpriscia/out/postG/"]
#  ["Run305040", "/eos/user/d/dapparu/thesis/hscp/DeDxCalib/out/"]
#  ["2017B", "/eos/user/d/dapparu/thesis/hscp/DeDxCalib/out2017B/"]
#  ["2017D", "/eos/user/d/dapparu/thesis/hscp/DeDxCalib/out2017D/"]
#  ["2017F", "/eos/user/d/dapparu/thesis/hscp/DeDxCalib/out2017F/"]
#  ["RunPreG", "/eos/user/j/jpriscia/out/preG/"],
#  ["Run278018", "/eos/user/j/jpriscia/out/278018/"],
#  ["Run278308", "/eos/user/j/jpriscia/out/278308/"],
#  ["Run279931", "/eos/user/j/jpriscia/out/279931/"],
##  ["Run280385", "/storage/data/cms/store/user/jozobec/dEdxCalib/280385/"],

#  ["MCGluino_M1000_f10", "Gluino_13TeV_M1000_f10"],
#  ["MCGluino_M1400_f10", "Gluino_13TeV_M1400_f10"],
#  ["MCGluino_M1800_f10", "Gluino_13TeV_M1800_f10"],
#  ["MCGluino_M1000_f50", "Gluino_13TeV_M1000_f50"],
#  ["MCGluino_M1400_f50", "Gluino_13TeV_M1400_f50"],
#  ["MCGluino_M1800_f50", "Gluino_13TeV_M1800_f50"],
#  ["MCGMStau_M494",      "GMStau_13TeV_M494"],
#  ["MCStop_M1000",       "Stop_13TeV_M1000"],
#  ["MCDYM2600Q2",        "DY_13TeV_M2600_Q2"],

##["2017A","/SingleMuon/Run2017A-09Aug2019_UL2017-v1/AOD"], #does not exist ? Commissionning
# ["2017B","/SingleMuon/Run2017B-09Aug2019_UL2017-v1/AOD"],
# ["2017B_all_janv21","/SingleMuon/Run2017B-09Aug2019_UL2017-v1/AOD"],
# ["2017C_5_janv21","/SingleMuon/Run2017E-09Aug2019_UL2017-v1/AOD"],
# ["2017C","/SingleMuon/Run2017C-09Aug2019_UL2017-v1/AOD"],
# ["2017D","/SingleMuon/Run2017D-09Aug2019_UL2017-v1/AOD"],
# ["2017E","/SingleMuon/Run2017E-09Aug2019_UL2017-v1/AOD"],
# ["2017F","/SingleMuon/Run2017F-09Aug2019_UL2017-v1/AOD"],
# ["2017F_100files_0_100","/SingleMuon/Run2017F-09Aug2019_UL2017-v1/AOD"],
# ["2017F_100files_100_200","/SingleMuon/Run2017F-09Aug2019_UL2017-v1/AOD"],
# ["2017F_100files_200_300","/SingleMuon/Run2017F-09Aug2019_UL2017-v1/AOD"],
# ["2017F_100files_300_400","/SingleMuon/Run2017F-09Aug2019_UL2017-v1/AOD"],
# ["2017F_400files_0_400","/SingleMuon/Run2017F-09Aug2019_UL2017-v1/AOD"],
# ["2017F_40files_0_40","/SingleMuon/Run2017F-09Aug2019_UL2017-v1/AOD"],
# ["2017Fsecond","/SingleMuon/Run2017F-09Aug2019_UL2017-v1/AOD"],
# ["2017Fthird","/SingleMuon/Run2017F-09Aug2019_UL2017-v1/AOD"],
# ["2017Ffourth","/SingleMuon/Run2017F-09Aug2019_UL2017-v1/AOD"],
# ["2017G","/SingleMuon/Run2017G-09Aug2019_UL2017-v1/AOD"],


 ["2017B_all_march21_500jobs_0320","/SingleMuon/Run2017B-09Aug2019_UL2017-v1/AOD"],
 ["2017C_all_march21_500jobs_0320","/SingleMuon/Run2017C-09Aug2019_UL2017-v1/AOD"],
 ["2017D_all_march21_500jobs_0320","/SingleMuon/Run2017D-09Aug2019_UL2017-v1/AOD"],
 ["2017E_all_march21_500jobs_0320","/SingleMuon/Run2017E-09Aug2019_UL2017-v1/AOD"],
 ["2017F_all_march21_500jobs_0320","/SingleMuon/Run2017F-09Aug2019_UL2017-v1/AOD"],

# ["2017B_all_march21_50jobs","/SingleMuon/Run2017B-09Aug2019_UL2017-v1/AOD"],
# ["2017C_all_march21_50jobs","/SingleMuon/Run2017C-09Aug2019_UL2017-v1/AOD"],
# ["2017D_all_march21_50jobs","/SingleMuon/Run2017D-09Aug2019_UL2017-v1/AOD"],
# ["2017E_all_march21_50jobs","/SingleMuon/Run2017E-09Aug2019_UL2017-v1/AOD"],
# ["2017F_all_march21_50jobs","/SingleMuon/Run2017F-09Aug2019_UL2017-v1/AOD"],

# ["2017B_all_march21_500jobs_0320","/SingleMuon/Run2017B-09Aug2019_UL2017-v1/AOD"],
# ["2017C_all_march21_500jobs_0320","/SingleMuon/Run2017C-09Aug2019_UL2017-v1/AOD"],
# ["2017D_all_march21_500jobs_0320","/SingleMuon/Run2017D-09Aug2019_UL2017-v1/AOD"],
# ["2017E_all_march21_500jobs_0320","/SingleMuon/Run2017E-09Aug2019_UL2017-v1/AOD"],
# ["2017F_all_march21_500jobs_0320","/SingleMuon/Run2017F-09Aug2019_UL2017-v1/AOD"],


# ["2017B_500files_march21_50jobs","/SingleMuon/Run2017B-09Aug2019_UL2017-v1/AOD"],
 
# ["2017B_all_march21_1000jobs","/SingleMuon/Run2017B-09Aug2019_UL2017-v1/AOD"],
# ["2017B_all_march21_200jobs_binsLumi","/SingleMuon/Run2017B-09Aug2019_UL2017-v1/AOD"],
]

isLocal = True  #allow to access data in Louvain from remote sites
#if(commands.getstatusoutput("hostname -f")[1].find("ucl.ac.be")!=-1): isLocal = True
os.system('rm -rf ~/x509_user_proxy/x509_proxy')

LoadJson(JSON)

if sys.argv[1]=='1':
        os.system("sh " + os.getcwd() + "/DeDxStudy.sh ") #just compile

	for DATASET in datasetList :
       #outdir = "/afs/cern.ch/work/d/dapparu/private/thesis/hscp/dEdxCalib/Results" + "/Histos2017/"+DATASET[0]+"/"
	   outdir =  os.getcwd() + "/Histos2017/"+DATASET[0]+"/"
	   os.system('mkdir -p ' + outdir)

	   JobName = "DEDXHISTO_"+DATASET[0]
	   FarmDirectory = "FARM_DEDXHISTO_"+DATASET[0]
           LaunchOnCondor.subTool = 'condor'
	   LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)

 	   FILELIST = []        
           if(DATASET[1][-1]=='/'): #file path is a directory, consider all files from the directory
              if(isLocal):
      	         FILELIST = LaunchOnCondor.GetListOfFiles('', DATASET[1]+'/*.root', '')
              else:
                 initProxy()
                 initCommand = 'export X509_USER_PROXY=~/private/x509_proxy; voms-proxy-init --noregen;'
                 LaunchOnCondor.Jobs_InitCmds = [initCommand]
                 #print initCommand+'lcg-ls -b -D srmv2 "srm://ingrid-se02.cism.ucl.ac.be:8444/srm/managerv2?SFN='+DATASET[1]+'" | xargs -I {} basename {}'
                 #print commands.getstatusoutput(initCommand+'lcg-ls -b -D srmv2 "srm://ingrid-se02.cism.ucl.ac.be:8444/srm/managerv2?SFN='+DATASET[1]+'" | xargs -I {} basename {}')
                 #LocalFileList = commands.getstatusoutput(initCommand+'lcg-ls -b -D srmv2 "srm://ingrid-se02.cism.ucl.ac.be:8444/srm/managerv2?SFN='+DATASET[1]+'" | xargs -I {} basename {}')[1].split('\n')

                 with open(DATASET[2], 'r') as fileList:
                   for f in fileList:
                     #if(f[-5:].find('.root')==-1):continue #only .root file considered
                     FILELIST += ["root://cms-xrd-global.cern.ch/"+DATASET[1].replace('/storage/data/cms/store/','/store/')+f.rstrip()]
    
           else: #file path is an HSCP sample name, use the name to run the job
             FILELIST += [DATASET[1]]
           
           print FILELIST
           initProxy()
           initCommand = 'export X509_USER_PROXY=~/private/x509_proxy;'
           LaunchOnCondor.Jobs_InitCmds = [initCommand]
           FILELIST = filesFromDataset(DATASET[1],-1) # run overall
           #FILELIST = filesFromDataset(DATASET[1],500) 
           print FILELIST
           for inFileList in getChunksFromList(FILELIST,max(1,len(FILELIST)/1000)): #50 jobs, this is a trade off between hadding time and processing time
              InputListCSV = ''
  	      for inFile in inFileList:
                 InputListCSV+= inFile + ','
              InputListCSV = InputListCSV[:-1] #remove the last duplicated comma
              LaunchOnCondor.SendCluster_Push  (["BASH", "sh " + os.getcwd() + "/DeDxStudy.sh " + InputListCSV + " out.root; mv out.root " + outdir+"dEdxHistos_%i.root" %  LaunchOnCondor.Jobs_Count ])
	   LaunchOnCondor.SendCluster_Submit()

elif sys.argv[1]=='2':
        for DATASET in datasetList :#+signalList :
           indir =  os.getcwd() + "/Histos2017/"+DATASET[0]+'/'
           #os.system('rm -f Histos_'+DATASET[0]+'.root')
           os.system('find ' + indir + '*.root  -type f -size +1024c | xargs hadd -f Histos_'+DATASET[0]+'.root')
	# finally merge all the runs into the histogram with data
	#os.system('rm -f Histos_Data.root')
	#os.system('hadd -f Histos_Data.root Histos_Run*.root')

elif sys.argv[1]=='3':
        os.system('sh MakePlot.sh')

else:
   print "Invalid argument"
   sys.exit()


#!/usr/bin/env python

import urllib
import string
import os
import sys
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor  
import glob
import time
import commands
def initProxy():
   if(not os.path.isfile(os.path.expanduser('~/x509_user_proxy/x509_proxy')) or ((time.time() - os.path.getmtime(os.path.expanduser('~/x509_user_proxy/x509_proxy')))>600 and  int(commands.getstatusoutput('(export X509_USER_PROXY=~/x509_user_proxy/x509_proxy;voms-proxy-init --noregen;voms-proxy-info -all) | grep timeleft | tail -n 1')[1].split(':')[2])<8 )):
      print "You are going to run on a sample over grid using either CRAB or the AAA protocol, it is therefore needed to initialize your grid certificate"
      os.system('mkdir -p ~/x509_user_proxy; voms-proxy-init --voms cms -valid 192:00 --out ~/x509_user_proxy/x509_proxy')#all must be done in the same command to avoid environement problems.  Note that the first sourcing is only needed in Louvain



"""
check that a file exist and is not corrupted
"""
def checkInputFile(url):
    if(url.startswith('/store')==True):
       url= 'root://eoscms//eos/cms'+url
    command_out = commands.getstatusoutput("root -l -b -q " + url)
    if(command_out[1].find("Error")>=0 or command_out[1].find("probably not closed")>=0 or command_out[1].find("Corrupted")>=0):return False
    return True


def getChunksFromList(MyList, n):
  return [MyList[x:x+n] for x in range(0, len(MyList), n)]

samples = [
#SampleName, generator config file, NJobs, NEvents/Job, NumberOfFinalEDMMergedfiles
  ['MC_13TeV_DYToMuMu' , 'ZMM_13TeV_TuneCUETP8M1_cfi', 4000, 2500, 15],  #carefull this is a huge number!
  ['MC_13TeV_WToLNu' , 'WToLNu_13TeV_pythia8_cff', 4000, 2500, 10],  #carefull this is a huge number
#  ['MC_13TeV_MinBias'  , 'MinBIas_13TeV_pythia8_TuneCUETP8M1_cfi', 1000, 1000, 1],
]


print("Initialize your grid proxy in case you need to access remote samples (for PU for instance)\n")
initProxy()

if sys.argv[1]=='1':  #GEN-SIM-DIGI-RECO-AOD in one shot
   #this is common to all samples, so we need to do it only once
   os.system('cmsDriver.py --filein file:step1.root --fileout step2.root --mc --eventcontent AODSIM --datatier GEN-SIM-DIGI-AOD --conditions MCRUN2_74_V9 --step HLT:@frozen25ns,RAW2DIGI,L1Reco,RECO --python_filename RECO_Template_cfg.py --magField 38T_PostLS1 --pileup_input "dbs:/MinBias_TuneCUETP8M1_13TeV-pythia8/RunIIWinter15GS-MCRUN2_71_V1-v1/GEN-SIM" --pileup 2015_25ns_Startup_PoissonOOTPU --geometry Extended2015 --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --no_exec -n -1')

   for S in samples:
      if('HSCP' in S[1]):
         os.system('cmsDriver.py ' + S[1] + ' --fileout file:step1.root --mc --eventcontent RAWSIM --datatier GEN-SIM-RAW --conditions MCRUN2_74_V9 --step GEN,SIM,DIGI,L1,DIGI2RAW --python_filename GEN_SIM_Template_cfg.py --magField 38T_PostLS1 --pileup_input "dbs:/MinBias_TuneCUETP8M1_13TeV-pythia8/RunIIWinter15GS-MCRUN2_71_V1-v1/GEN-SIM" --pileup 2015_25ns_Startup_PoissonOOTPU  --geometry Extended2015 --customise SimG4Core/CustomPhysics/Exotica_HSCP_SIM_cfi.customise,SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --no_exec -n 10')
      else:
         os.system('cmsDriver.py ' + S[1] + ' --fileout file:step1.root --mc --eventcontent RAWSIM --datatier GEN-SIM-RAW --conditions MCRUN2_74_V9 --step GEN,SIM,DIGI,L1,DIGI2RAW --python_filename GEN_SIM_Template_cfg.py --magField 38T_PostLS1 --pileup_input "dbs:/MinBias_TuneCUETP8M1_13TeV-pythia8/RunIIWinter15GS-MCRUN2_71_V1-v1/GEN-SIM" --pileup 2015_25ns_Startup_PoissonOOTPU --geometry Extended2015 --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --no_exec -n 10')

      with open("GEN_SIM_Template_cfg.py", "a") as f:
         f.write('process.RandomNumberGeneratorService.generator.initialSeed = cms.untracked.uint32(XXX_SEED_XXX)\n')
         f.write('process.RandomNumberGeneratorService.mix.initialSeed = cms.untracked.uint32(XXX_SEED_XXX)\n')
         f.write('process.maxEvents.input = cms.untracked.int32(XXX_NEVENTS_XXX)\n')
#         f.write('process.RAWSIMoutput.fileName = cms.untracked.string("file:XXX_OUTPUT_XXX_XXX_I_XXX.root")\n')
         f.write('process.source.firstLuminosityBlock =  cms.untracked.uint32(1XXX_I_XXX)\n')
         f.close()


      with open("HSCPEDM_Template_cfg.py", "w") as f:
         f.write("import sys, os\n")
         f.write("import FWCore.ParameterSet.Config as cms\n")
         f.write("\n")
         if('HSCP' in S[1]):
            f.write("isSignal = True\n")
            f.write("isBckg = False\n")
         else:
            f.write("isSignal = False\n")
            f.write("isBckg = True\n")
         f.write("isData = False\n")
         f.write("isSkimmedSample = False\n")
         f.write("GTAG = 'MCRUN2_74_V8'\n")
         f.write("OUTPUTFILE = 'XXX_OUTPUT_XXX_XXX_I_XXX.root'\n")
         f.write("InputFileList = cms.untracked.vstring()\n")
         f.write("\n")
         f.write("InputFileList.extend(['file:"+"step2.root"+"'])\n")

         f.write("\n")
         f.write("#main EDM tuple cfg that depends on the above parameters\n")
         f.write("execfile( os.path.expandvars('${CMSSW_BASE}/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/HSCParticleProducer_cfg.py') )\n")
         f.close()


      JobName = S[0]+"_SIMEDM"
      FarmDirectory = "FARM_"+JobName
      LaunchOnCondor.Jobs_NEvent = S[3]
      LaunchOnCondor.Jobs_Skip = 0
      LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
      LaunchOnCondor.Jobs_InitCmds = ['export X509_USER_PROXY=~/x509_user_proxy/x509_proxy; voms-proxy-init --noregen;']
      for i in range(0,S[2]):
           LaunchOnCondor.Jobs_Count = i
           LaunchOnCondor.Jobs_Skip+=LaunchOnCondor.Jobs_NEvent
           LaunchOnCondor.SendCluster_Push  (["CMSSW", ["GEN_SIM_Template_cfg.py", "RECO_Template_cfg.py", "HSCPEDM_Template_cfg.py"] ])
           LaunchOnCondor.Jobs_FinalCmds = ['rm step1.root; rm step2.root']
      LaunchOnCondor.SendCluster_Submit()


elif sys.argv[1]=='2':  #MergeAll EDM files into one
   for S in samples:
      InputFiles = LaunchOnCondor.GetListOfFiles('"file:',os.getcwd()+"/FARM_"+S[0]+"_SIMEDM/outputs/*.root",'"')
      chunks = getChunksFromList(InputFiles, (len(InputFiles)/S[4])+1)
      for I in range(0,S[4]):
         if(I>=len(chunks)):continue
         JobName = S[0]
         if(I>0):JobName+="_" + str(I)
         FarmDirectory = "FARM_"+JobName
         LaunchOnCondor.SendCMSMergeJob(FarmDirectory, JobName, chunks[I],  '"'+JobName+'.root"', '"keep *"')



elif sys.argv[1]=='3':  #Transfert final EDM files from your place to CERN CMST3
   for S in samples:
      for I in range(0,S[4]):
         fileName = S[0]
         if(I>0):fileName+="_" + str(I)
#         os.system('lcg-cp --verbose -b -D srmv2 "srm://ingrid-se02.cism.ucl.ac.be:8444/srm/managerv2?SFN=/storage/data/cms/store/user/quertenmont/15_05_20_HSCP_AODSIM/'+fileName+'.root" "srm://srm-eoscms.cern.ch:8443/srm/v2/server?SFN=/eos/cms//store/cmst3/user/querten/15_03_25_HSCP_Run2EDMFiles/'+fileName+'.root" &> LOG_'+fileName+'.log &')
         os.system('lcg-cp --verbose -b -D srmv2 "srm://ingrid-se02.cism.ucl.ac.be:8444/srm/managerv2?SFN=/storage/data/cms/users/quertenmont/HSCP/2015/'+fileName+'.root" "srm://srm-eoscms.cern.ch:8443/srm/v2/server?SFN=/eos/cms//store/cmst3/user/querten/15_03_25_HSCP_Run2EDMFiles/'+fileName+'.root" &> LOG_'+fileName+'.log &')




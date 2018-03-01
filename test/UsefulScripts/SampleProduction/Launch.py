#!/usr/bin/env python

import urllib
import string
import os
import sys
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor  
import glob
import time
"""
check that a file exist and is not corrupted
"""
def checkInputFile(url):
    if(url.startswith('/store')==True):
       url= 'root://eoscms//eos/cms'+url
    command_out = commands.getstatusoutput("root -l -b -q " + url)
    if(command_out[1].find("Error")>=0 or command_out[1].find("probably not closed")>=0 or command_out[1].find("Corrupted")>=0):return False
    return True


samples = [
#SampleName, generator config file, NJobs, NEvents/Job
#  ['Gluino_13TeV_M200' , 'Configuration/GenProduction/python/ThirteenTeV/HSCPgluino_M_200_TuneCUETP8M1_13TeV_pythia8_cff.py', 100, 100],
#  ['Gluino_13TeV_M600' , 'Configuration/GenProduction/python/ThirteenTeV/HSCPgluino_M_600_TuneCUETP8M1_13TeV_pythia8_cff.py', 100, 100],
#  ['Gluino_13TeV_M1200', 'Configuration/GenProduction/python/ThirteenTeV/HSCPgluino_M_1200_TuneCUETP8M1_13TeV_pythia8_cff.py', 100, 100],
#  ['Gluino_13TeV_M1800', 'Configuration/GenProduction/python/ThirteenTeV/HSCPgluino_M_1800_TuneCUETP8M1_13TeV_pythia8_cff.py', 100, 100],
#  ['Gluino_13TeV_M2400', 'Configuration/GenProduction/python/ThirteenTeV/HSCPgluino_M_2400_TuneCUETP8M1_13TeV_pythia8_cff.py', 100, 100],

  ['Stop_13TeV_M200' , 'Configuration/GenProduction/python/ThirteenTeV/HSCPstop_M_200_TuneCUETP8M1_13TeV_pythia8_cff.py', 1, 100], #FIXME
#  ['Stop_13TeV_M600' , 'Configuration/GenProduction/python/ThirteenTeV/HSCPstop_M_600_TuneCUETP8M1_13TeV_pythia8_cff.py', 100, 100],
#  ['Stop_13TeV_M1200', 'Configuration/GenProduction/python/ThirteenTeV/HSCPstop_M_1200_TuneCUETP8M1_13TeV_pythia8_cff.py', 100, 100],
  ['Stop_13TeV_M1800', 'Configuration/GenProduction/python/ThirteenTeV/HSCPstop_M_1800_TuneCUETP8M1_13TeV_pythia8_cff.py', 1, 100],  #FIXME
#  ['Stop_13TeV_M2400', 'Configuration/GenProduction/python/ThirteenTeV/HSCPstop_M_2400_TuneCUETP8M1_13TeV_pythia8_cff.py', 100, 100],

#  ['GMStau_13TeV_M156', 'Configuration/GenProduction/python/ThirteenTeV/HSCPgmstau_M_156_TuneZ2star_13TeV_pythia6_cff.py', 100, 100],
#  ['GMStau_13TeV_M308', 'Configuration/GenProduction/python/ThirteenTeV/HSCPgmstau_M_308_TuneZ2star_13TeV_pythia6_cff.py', 100, 100],
#  ['GMStau_13TeV_M651', 'Configuration/GenProduction/python/ThirteenTeV/HSCPgmstau_M_651_TuneZ2star_13TeV_pythia6_cff.py', 100, 100],
#  ['GMStau_13TeV_M1029', 'Configuration/GenProduction/python/ThirteenTeV/HSCPgmstau_M_1029_TuneZ2star_13TeV_pythia6_cff.py', 100, 100],
#  ['GMStau_13TeV_M1599', 'Configuration/GenProduction/python/ThirteenTeV/HSCPgmstau_M_1599_TuneZ2star_13TeV_pythia6_cff.py', 100, 100],

#  ['PPStau_13TeV_M156', 'Configuration/GenProduction/python/ThirteenTeV/HSCPppstau_M_156_TuneZ2star_13TeV_pythia6_cff.py', 100, 100],
#  ['PPStau_13TeV_M308', 'Configuration/GenProduction/python/ThirteenTeV/HSCPppstau_M_308_TuneZ2star_13TeV_pythia6_cff.py', 100, 100],
#  ['PPStau_13TeV_M651', 'Configuration/GenProduction/python/ThirteenTeV/HSCPppstau_M_651_TuneZ2star_13TeV_pythia6_cff.py', 100, 100],
  ['PPStau_13TeV_M1029', 'Configuration/GenProduction/python/ThirteenTeV/HSCPppstau_M_1029_TuneZ2star_13TeV_pythia6_cff.py', 100, 100],
#  ['PPStau_13TeV_M1599', 'Configuration/GenProduction/python/ThirteenTeV/HSCPppstau_M_1599_TuneZ2star_13TeV_pythia6_cff.py', 100, 100],
]

doTimeShift = False
if sys.argv[1]=='1':  #GEN-SIM-DIGI-RECO-AOD in one shot
   #this is common to all samples, so we need to do it only once
#   os.system('cmsDriver.py --filein file:step1.root --fileout XXX_OUTPUT_XXX_XXX_I_XXX.root --mc --eventcontent AODSIM --datatier GEN-SIM-DIGI-AOD --conditions auto:run2_mc --step HLT:GRun,RAW2DIGI,L1Reco,RECO --python_filename RECO_Template_cfg.py --magField 38T_PostLS1 --geometry Extended2015 --era Run2_25ns --no_exec -n -1')
   os.system('cmsDriver.py --filein file:step1.root --fileout XXX_OUTPUT_XXX_XXX_I_XXX.root --mc --eventcontent AODSIM --datatier GEN-SIM-DIGI-AOD --conditions auto:run2_mc --step HLT:GRun,RAW2DIGI,L1Reco,RECO --python_filename RECO_Template_cfg.py  --era Run2_2016 --no_exec -n -1')

   for S in samples:
      if('HSCP' in S[1]):
         os.system('cmsDriver.py ' + S[1] + ' --fileout file:step1.root --mc --eventcontent RAWSIM --datatier GEN-SIM-RAW --conditions auto:run2_mc --step GEN,SIM,DIGI,L1,DIGI2RAW --python_filename GEN_SIM_Template_cfg.py  --era Run2_2016 --no_exec -n 10')
      else:
         os.system('cmsDriver.py ' + S[1] + ' --fileout file:step1.root --mc --eventcontent RAWSIM --datatier GEN-SIM-RAW --conditions auto:run2_mc --step GEN,SIM,DIGI,L1,DIGI2RAW --python_filename GEN_SIM_Template_cfg.py  --era Run2_2016 --no_exec -n 10')

      lines = ""
      with open("GEN_SIM_Template_cfg.py", "r") as f:
         linesList = f.readlines()
         for line in linesList:
            lines+=line;

         f.close()

      lines += ('process.RandomNumberGeneratorService.generator.initialSeed = cms.untracked.uint32(XXX_SEED_XXX)\n')
      lines += ('process.RandomNumberGeneratorService.mix.initialSeed = cms.untracked.uint32(XXX_SEED_XXX)\n')
      lines += ('process.maxEvents.input = cms.untracked.int32(XXX_NEVENTS_XXX)\n')
      #lines += ('process.RAWSIMoutput.fileName = cms.untracked.string("file:XXX_OUTPUT_XXX_XXX_I_XXX.root")\n')
      lines += ('process.source.firstLuminosityBlock =  cms.untracked.uint32(1XXX_I_XXX)\n')


      if(doTimeShift==True):
         repla  = 'process.shifter = cms.EDProducer("SimHitShifterRun2",\n'
         repla += '   ShiftFileName = cms.untracked.string("/nfs/home/fynu/quertenmont/scratch/15_05_13_HSCP_SampleProd/CMSSW_7_4_14/src/SUSYBSMAnalysis/HSCP/data/MuonTimeShift_2015.txt"),\n'
         repla += ')\n'
         repla += '\n'
         repla += 'process.mix.mixObjects.mixSH.input = cms.VInputTag(  # note that this list needs to be in the same order as the subdets\n'
         repla += '        #cms.InputTag("g4SimHits","BSCHits"), cms.InputTag("g4SimHits","BCM1FHits"), cms.InputTag("g4SimHits","PLTHits"), cms.InputTag("g4SimHits","FP420SI"),\n'
         repla += '        cms.InputTag("shifter","MuonCSCHits"), cms.InputTag("shifter","MuonDTHits"), cms.InputTag("shifter","MuonRPCHits"), \n'
         repla += '        #cms.InputTag("g4SimHits","TotemHitsRP"), cms.InputTag("g4SimHits","TotemHitsT1"), cms.InputTag("g4SimHits","TotemHitsT2Gem"),\n'
         repla += '        cms.InputTag("g4SimHits","TrackerHitsPixelBarrelHighTof"), cms.InputTag("g4SimHits","TrackerHitsPixelBarrelLowTof"), \n'
         repla += '        cms.InputTag("g4SimHits","TrackerHitsPixelEndcapHighTof"), cms.InputTag("g4SimHits","TrackerHitsPixelEndcapLowTof"), cms.InputTag("g4SimHits","TrackerHitsTECHighTof"), cms.InputTag("g4SimHits","TrackerHitsTECLowTof"), cms.InputTag("g4SimHits","TrackerHitsTIBHighTof"), \n'
         repla += '        cms.InputTag("g4SimHits","TrackerHitsTIBLowTof"), cms.InputTag("g4SimHits","TrackerHitsTIDHighTof"), cms.InputTag("g4SimHits","TrackerHitsTIDLowTof"), cms.InputTag("g4SimHits","TrackerHitsTOBHighTof"), cms.InputTag("g4SimHits","TrackerHitsTOBLowTof"))\n'
         repla += '\n'
         repla += '\n'
         repla += 'process.simMuonCSCDigis.InputCollection = cms.string(\'shifterMuonCSCHits\')\n'
         repla += 'process.simMuonDTDigis.InputCollection = cms.string(\'shifterMuonDTHits\')\n'
         repla += 'process.simMuonRPCDigis.InputCollection = cms.string(\'shifterMuonRPCHits\')\n'
         repla += '\n'
         repla += '# Path and EndPath definitions\n'

      
         lines = lines.replace('# Path and EndPath definitions', repla)
         lines = lines.replace('cms.Path(process.psim)', 'cms.Path(process.psim * process.shifter)')


      with open("GEN_SIM_Template_cfg.py", "w") as f:
         f.write(lines)
         f.close()

      JobName = S[0]+"_SIMAOD"
      FarmDirectory = "FARM_"+JobName
      LaunchOnCondor.Jobs_NEvent = S[3]
      LaunchOnCondor.Jobs_Skip = 0
      LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
      for i in range(0,S[2]):
           LaunchOnCondor.Jobs_Count = i
           LaunchOnCondor.Jobs_Skip+=LaunchOnCondor.Jobs_NEvent
           LaunchOnCondor.SendCluster_Push  (["CMSSW", ["GEN_SIM_Template_cfg.py", "RECO_Template_cfg.py"] ])
           #LaunchOnCondor.Jobs_FinalCmds = ['mv step1.root ' + os.getcwd()+"/FARM_"+S[0]+"_SIMAOD/outputs/"+S[0]+"_SIM_%04i.root" % i]
      LaunchOnCondor.SendCluster_Submit()




elif sys.argv[1]=='2':  #AOD --> EDM files
   for S in samples:
      JobName = S[0]
      FarmDirectory = "FARM_EDM"
      LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)

      f = open("HSCPEDM_cfg.py", "w")
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
      f.write("OUTPUTFILE = 'XXX_OUTPUT_XXX.root'\n")
      f.write("InputFileList = cms.untracked.vstring()\n")
      f.write("\n")
      for i in range(0,S[2]):
         inFile = os.getcwd()+"/FARM_"+S[0]+"_SIMAOD/outputs/"+S[0]+"_SIMAOD_%04i.root" % i
         f.write("InputFileList.extend(['file:"+inFile+"'])\n")
#         if(checkInputFile(inFile) ):
#            f.write("InputFileList.extend(['file:"+inFile+"'])\n")
#         else:
#            print "missing "+ inFile
      f.write("\n")
      f.write("#main EDM tuple cfg that depends on the above parameters\n")
      f.write("execfile( os.path.expandvars('${CMSSW_BASE}/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/HSCParticleProducer_cfg.py') )\n")
      f.close()

      LaunchOnCondor.SendCluster_Push  (["CMSSW", "HSCPEDM_cfg.py" ])
      LaunchOnCondor.SendCluster_Submit()

elif sys.argv[1]=='3':  #Transfert final EDM files from your place to CERN CMST3
   for S in samples:
      os.system('lcg-cp --verbose -b -D srmv2 "srm://ingrid-se02.cism.ucl.ac.be:8444/srm/managerv2?SFN=/storage/data/cms/store/user/quertenmont/15_05_20_HSCP_AODSIM/'+S[0]+'.root" "srm://srm-eoscms.cern.ch:8443/srm/v2/server?SFN=/eos/cms//store/cmst3/user/querten/15_03_25_HSCP_Run2EDMFiles/'+S[0]+'.root" &> LOG_'+S[0]+'.log &')

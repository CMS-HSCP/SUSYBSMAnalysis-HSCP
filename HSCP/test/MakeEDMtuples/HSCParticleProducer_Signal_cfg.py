import sys, os
import FWCore.ParameterSet.Config as cms

isSignal = True
isBckg = False
isData = False
isSkimmedSample = False
GTAG = 'MCRUN2_74_V7::All'
OUTPUTFILE = 'HSCP.root'

InputFileList = cms.untracked.vstring()

#debug input files 
#this list is overwritten by CRAB
InputFileList = cms.untracked.vstring(
#   'root://cmseos.fnal.gov//eos/uscms/store/user/aackert/HSCP/AODFiles/HSCPgluinoM1000_AOD-SIM.root',
#   'root://cmseos.fnal.gov//eos/uscms/store/user/aackert/HSCP/AODFiles/HSCPgluinoM1500_AOD-SIM.root',
#   'root://cmseos.fnal.gov//eos/uscms/store/user/aackert/HSCP/AODFiles/HSCPstopM500_AOD-SIM.root', 
#   'root://cmseos.fnal.gov//eos/uscms/store/user/aackert/HSCP/AODFiles/HSCPstopM900_AOD-SIM.root',
#    'root://cmseos.fnal.gov//eos/uscms/store/user/aackert/HSCP/AODFiles/HSCPgmstauM308_AOD-SIM.root',
#   'root://cmseos.fnal.gov//eos/uscms/store/user/aackert/HSCP/AODFiles/HSCPgmstauM494_AOD-SIM.root',

#'root://cmseos.fnal.gov//eos/uscms/store/user/aackert/HSCP/AODFiles/HSCPgluinoM1500_AOD-SIM_740_MCRUN2_74_v7.root',

#'root://cmseos.fnal.gov///eos/uscms/store/user/aackert/HSCP/AODFiles/HSCPmchamp3_M2600_AOD-SIM_740_MCRUN2_74_V7.root',
'root://cmseos.fnal.gov//eos/uscms/store/user/aackert/HSCP/AODgen/condorjdls/step2_condortest.root',
)


#main EDM tuple cfg that depends on the above parameters
execfile( os.path.expandvars('${CMSSW_BASE}/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/HSCParticleProducer_cfg.py') )

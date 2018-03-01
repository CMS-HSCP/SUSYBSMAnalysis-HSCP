import FWCore.ParameterSet.Config as cms

process = cms.Process("HSCPAnalysis")

#The following parameters need to be provided
#isSignal, isBckg, isData, isSkimmedSample, GTAG, InputFileList
#isSignal = True
#isBckg = False
#isData = False
#isSkimmedSample = False
#GTAG = 'START72_V1::All'

isSignal = False
isBckg = False
isData = True
isSkimmedSample = False
GTAG = 'GR_R_72_V2::All'

#debug input files 
#this list is overwritten by CRAB
if('XXX_SAMPLE_XXX'=='MU'):
   InputFileList = cms.untracked.vstring(
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/00346527-C91C-E411-AB5E-02163E00ECEF.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/00ADAB1A-BC1C-E411-8EF1-002590494C40.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/0200D6AF-D41C-E411-8420-0025904B0FC0.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/0268317F-AB1C-E411-8022-02163E00ECFB.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/02F71E4D-B61C-E411-888D-02163E00CFB4.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/02FCA65D-AF1C-E411-BED1-18A90555637A.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/0405F2E2-B01C-E411-BC20-02163E00E5B2.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/041A1E18-C21C-E411-A77A-00259029EF3E.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/044F878E-E91C-E411-B9D8-02163E00CAA2.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/04C9C6B8-D31C-E411-88E9-003048F0E7BE.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/0634987C-A91C-E411-A9B5-02163E009C1E.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/0637E986-BE1C-E411-ACFF-02163E00EF94.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/0644EC20-C91C-E411-A376-02163E008EEA.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/068E6015-CB1C-E411-96C8-02163E00E95C.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/06D99999-BD1C-E411-896E-0025B3203748.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/08361980-BC1C-E411-97CB-003048C9C1D4.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/0836D993-B61C-E411-B245-02163E009BA7.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/08A1A887-C41C-E411-9BCB-02163E00FEC3.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/08B88B4E-D11C-E411-9EDB-02163E00B7A3.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/08CE58CE-EB1C-E411-BABA-02163E010110.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/08D1A881-9F1C-E411-8E09-02163E00ECE6.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/08F06953-BB1C-E411-8F27-003048C9C1D0.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleMu/RECO/GR_R_72_V2_frozenHLT_RelVal_mu2012D-v1/00000/0A44FE81-A71C-E411-862A-02163E0104C0.root',
   )
elif('XXX_SAMPLE_XXX'=='MET'):
   InputFileList = cms.untracked.vstring(
   '/store/relval/CMSSW_7_2_0_pre3/MET/RECO/GR_R_72_V2_frozenHLT_RelVal_met2012A-v1/00000/00F06DCF-A21C-E411-8188-02163E00BA3A.root',
   '/store/relval/CMSSW_7_2_0_pre3/MET/RECO/GR_R_72_V2_frozenHLT_RelVal_met2012A-v1/00000/10EBA9AC-9C1C-E411-9F46-02163E00EB28.root',
   '/store/relval/CMSSW_7_2_0_pre3/MET/RECO/GR_R_72_V2_frozenHLT_RelVal_met2012A-v1/00000/1ABAB2CE-A21C-E411-A5A7-0025904B26B0.root',
   '/store/relval/CMSSW_7_2_0_pre3/MET/RECO/GR_R_72_V2_frozenHLT_RelVal_met2012A-v1/00000/1C7BCC9D-9C1C-E411-AFFF-02163E008CD6.root',
   '/store/relval/CMSSW_7_2_0_pre3/MET/RECO/GR_R_72_V2_frozenHLT_RelVal_met2012A-v1/00000/20F5465B-A01C-E411-BA2F-02163E00F444.root',
   '/store/relval/CMSSW_7_2_0_pre3/MET/RECO/GR_R_72_V2_frozenHLT_RelVal_met2012A-v1/00000/24A35A0B-9E1C-E411-8744-02163E00C765.root',
   '/store/relval/CMSSW_7_2_0_pre3/MET/RECO/GR_R_72_V2_frozenHLT_RelVal_met2012A-v1/00000/2A75D387-A31C-E411-86EC-002481E0DC6C.root',
   '/store/relval/CMSSW_7_2_0_pre3/MET/RECO/GR_R_72_V2_frozenHLT_RelVal_met2012A-v1/00000/30146E2D-9B1C-E411-8B1B-02163E00EAA9.root',
   '/store/relval/CMSSW_7_2_0_pre3/MET/RECO/GR_R_72_V2_frozenHLT_RelVal_met2012A-v1/00000/3A3491B5-9D1C-E411-9B26-02163E00CE93.root',
   '/store/relval/CMSSW_7_2_0_pre3/MET/RECO/GR_R_72_V2_frozenHLT_RelVal_met2012A-v1/00000/3AE32D1D-9A1C-E411-ADDA-002481E0E5E8.root',
   )
elif('XXX_SAMPLE_XXX'=='ELEC'):
   InputFileList = cms.untracked.vstring(
   '/store/relval/CMSSW_7_2_0_pre3/SingleElectron/RECO/GR_R_72_V2_bis_RelVal_electron2012D-v1/00000/0085343A-6B29-E411-851C-02163E00E6BB.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleElectron/RECO/GR_R_72_V2_bis_RelVal_electron2012D-v1/00000/00A09E7F-7829-E411-8A30-001E67ABF674.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleElectron/RECO/GR_R_72_V2_bis_RelVal_electron2012D-v1/00000/00B51569-8A29-E411-BCCE-00237DDC5C24.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleElectron/RECO/GR_R_72_V2_bis_RelVal_electron2012D-v1/00000/00BEDA60-6F29-E411-901F-02163E00E911.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleElectron/RECO/GR_R_72_V2_bis_RelVal_electron2012D-v1/00000/00C9021B-5F29-E411-9AB1-02163E00A40D.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleElectron/RECO/GR_R_72_V2_bis_RelVal_electron2012D-v1/00000/00E087C5-6D29-E411-80C3-02163E00EA47.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleElectron/RECO/GR_R_72_V2_bis_RelVal_electron2012D-v1/00000/02177604-8F29-E411-AF45-02163E00F48B.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleElectron/RECO/GR_R_72_V2_bis_RelVal_electron2012D-v1/00000/02226B31-7B29-E411-BBA1-02163E00E683.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleElectron/RECO/GR_R_72_V2_bis_RelVal_electron2012D-v1/00000/0275D623-4829-E411-987D-FA163E6DB06F.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleElectron/RECO/GR_R_72_V2_bis_RelVal_electron2012D-v1/00000/02BAEBC6-6329-E411-86D5-02163E00A0E5.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleElectron/RECO/GR_R_72_V2_bis_RelVal_electron2012D-v1/00000/04F2CBCE-6929-E411-B67C-002481E0D084.root',
   '/store/relval/CMSSW_7_2_0_pre3/SingleElectron/RECO/GR_R_72_V2_bis_RelVal_electron2012D-v1/00000/04FDD18E-8429-E411-A1BE-02163E00F51D.root',
   )
else:
   GTAG = 'START72_V1::All'
   InputFileList = cms.untracked.vstring(
   'file:/afs/cern.ch/user/q/querten/workspace/public/14_08_12_Run2HSCP/CMSSW_7_2_X_2014-08-18-0200/src/SUSYBSMAnalysis/HSCP/test/UsefulScripts/SampleProduction/FARM_XXX_SAMPLE_XXX_RECO/outputs/XXX_SAMPLE_XXX_RECO_0000.root',
   'file:/afs/cern.ch/user/q/querten/workspace/public/14_08_12_Run2HSCP/CMSSW_7_2_X_2014-08-18-0200/src/SUSYBSMAnalysis/HSCP/test/UsefulScripts/SampleProduction/FARM_XXX_SAMPLE_XXX_RECO/outputs/XXX_SAMPLE_XXX_RECO_0001.root',
   'file:/afs/cern.ch/user/q/querten/workspace/public/14_08_12_Run2HSCP/CMSSW_7_2_X_2014-08-18-0200/src/SUSYBSMAnalysis/HSCP/test/UsefulScripts/SampleProduction/FARM_XXX_SAMPLE_XXX_RECO/outputs/XXX_SAMPLE_XXX_RECO_0002.root',
   'file:/afs/cern.ch/user/q/querten/workspace/public/14_08_12_Run2HSCP/CMSSW_7_2_X_2014-08-18-0200/src/SUSYBSMAnalysis/HSCP/test/UsefulScripts/SampleProduction/FARM_XXX_SAMPLE_XXX_RECO/outputs/XXX_SAMPLE_XXX_RECO_0003.root',
   'file:/afs/cern.ch/user/q/querten/workspace/public/14_08_12_Run2HSCP/CMSSW_7_2_X_2014-08-18-0200/src/SUSYBSMAnalysis/HSCP/test/UsefulScripts/SampleProduction/FARM_XXX_SAMPLE_XXX_RECO/outputs/XXX_SAMPLE_XXX_RECO_0004.root',
   'file:/afs/cern.ch/user/q/querten/workspace/public/14_08_12_Run2HSCP/CMSSW_7_2_X_2014-08-18-0200/src/SUSYBSMAnalysis/HSCP/test/UsefulScripts/SampleProduction/FARM_XXX_SAMPLE_XXX_RECO/outputs/XXX_SAMPLE_XXX_RECO_0005.root',
   'file:/afs/cern.ch/user/q/querten/workspace/public/14_08_12_Run2HSCP/CMSSW_7_2_X_2014-08-18-0200/src/SUSYBSMAnalysis/HSCP/test/UsefulScripts/SampleProduction/FARM_XXX_SAMPLE_XXX_RECO/outputs/XXX_SAMPLE_XXX_RECO_0006.root',
   'file:/afs/cern.ch/user/q/querten/workspace/public/14_08_12_Run2HSCP/CMSSW_7_2_X_2014-08-18-0200/src/SUSYBSMAnalysis/HSCP/test/UsefulScripts/SampleProduction/FARM_XXX_SAMPLE_XXX_RECO/outputs/XXX_SAMPLE_XXX_RECO_0007.root',
   'file:/afs/cern.ch/user/q/querten/workspace/public/14_08_12_Run2HSCP/CMSSW_7_2_X_2014-08-18-0200/src/SUSYBSMAnalysis/HSCP/test/UsefulScripts/SampleProduction/FARM_XXX_SAMPLE_XXX_RECO/outputs/XXX_SAMPLE_XXX_RECO_0008.root',
   'file:/afs/cern.ch/user/q/querten/workspace/public/14_08_12_Run2HSCP/CMSSW_7_2_X_2014-08-18-0200/src/SUSYBSMAnalysis/HSCP/test/UsefulScripts/SampleProduction/FARM_XXX_SAMPLE_XXX_RECO/outputs/XXX_SAMPLE_XXX_RECO_0009.root',
   )


process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Services_cff')

process.options   = cms.untracked.PSet(
      wantSummary = cms.untracked.bool(True),
      SkipEvent = cms.untracked.vstring('ProductNotFound'),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.source = cms.Source("PoolSource",
   fileNames = InputFileList,
   inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")
)
if(isSignal): process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")


#for i in range(0,25):
#   process.source.fileNames.extend(["file:/afs/cern.ch/user/q/querten/workspace/public/14_08_12_Run2HSCP/CMSSW_7_2_X_2014-08-18-0200/src/SUSYBSMAnalysis/HSCP/test/BuildHSCParticles/Signals/../../../../../SampleProd/FARM_RECO/outputs/gluino1TeV_RECO_%04i.root" % i])

process.GlobalTag.globaltag = GTAG

process.HSCPTuplePath = cms.Path()

########################################################################
#Run the Skim sequence if necessary
if(not isSkimmedSample):
   process.nEventsBefSkim  = cms.EDProducer("EventCountProducer")

   process.load('Configuration.Skimming.PDWG_EXOHSCP_cff')
   process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
   process.HSCPTrigger = process.hltHighLevel.clone()
   process.HSCPTrigger.TriggerResultsTag = cms.InputTag( "TriggerResults", "", "HLT" )
   process.HSCPTrigger.andOr = cms.bool( True ) #OR
   process.HSCPTrigger.throw = cms.bool( False )
   if(isData):
      process.HSCPTrigger.HLTPaths = [
         "HLT_*_dEdx*",
         "HLT_Mu40_eta2p1*",
         "HLT_Mu50_eta2p1*",
         "HLT_HT650_*",
         "HLT_MET80_*",
         "HLT_L2Mu*MET*",
         "HLT_L2Mu*NoBPTX*",
         "HLT_PFMET150_*",
      ]
      process.HSCPTrigger.HLTPaths = ["*"]
   elif(isBckg):
      #to be updated to Run2 Triggers, in the meanwhile keep all of them to study trigger efficiency
      process.HSCPTrigger.HLTPaths = ["*"]
   else:
      #do not apply trigger filter on signal
      process.HSCPTrigger.HLTPaths = ["*"]  
   
process.HSCPTuplePath += process.nEventsBefSkim + process.HSCPTrigger + process.exoticaHSCPSeq

########################################################################

#Run the HSCP EDM-tuple Sequence on skimmed sample
process.nEventsBefEDM   = cms.EDProducer("EventCountProducer")
process.HSCPTuplePath += process.nEventsBefEDM

process.generalTracksSkim.ptMin = cms.double(XXX_PT_XXX)
process.generalTracksSkim.nHitMin = XXX_NH_XXX
process.generalTracksSkim.chi2nMax = 10.0
process.DedxFilter.trkPtMin = cms.double(XXX_PT_XXX)
process.DedxFilter.dedxMin =cms.double(XXX_DEDX_XXX)
process.DedxFilter.dedxMaxLeft =cms.double(2.8)
#process.DedxFilter.etaMin= cms.double(-2.4)
#process.DedxFilter.etaMax= cms.double(-2.4)
process.DedxFilter.SAMuPtMin =  cms.double(9999)


########################################################################

#make the pool output
process.Out = cms.OutputModule("PoolOutputModule",
     outputCommands = cms.untracked.vstring(
         "drop *",
         "keep EventAux_*_*_*",
         "keep LumiSummary_*_*_*",
         "keep edmMergeableCounter_*_*_*",
         "keep PileupSummaryInfos_addPileupInfo_*_*"
    ),
    fileName = cms.untracked.string('HSCP.root'),
    SelectEvents = cms.untracked.PSet(
       SelectEvents = cms.vstring('')
    ),
)

process.Out.SelectEvents.SelectEvents =  cms.vstring('HSCPTuplePath')


########################################################################
#process.Out.outputCommands = process.EXOHSCPSkim_EventContent.outputCommands

#schedule the sequence
process.endPath1 = cms.EndPath(process.Out)
process.schedule = cms.Schedule(process.HSCPTuplePath, process.endPath1)



process.maxEvents.input = cms.untracked.int32(XXX_NEVENTS_XXX)
process.Out.fileName = cms.untracked.string("file:XXX_OUTPUT_XXX_XXX_SAMPLE_XXX_pTXXX_PT_XXX_nHXXX_NH_XXX_dEdxXXX_DEDX_XXX.root")


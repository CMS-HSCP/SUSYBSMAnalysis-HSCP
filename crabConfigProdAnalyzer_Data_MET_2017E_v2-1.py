from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Analysis_MET_UL2017E'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HSCParticleProducerAnalyzer_cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.JobType.inputFiles = ['HSCP/data/Data13TeVGains_v2.root','HSCP/data/dEdxTemplate_harm2_SO_in_noC_CCC_MG_2017E.root','HSCP/data/CMS_GeomTree.root','HSCP/data/MuonTimeOffset.txt']

test=False
if test:
    config.Data.userInputFiles = ['root://cms-xrd-global.cern.ch//store/data/Run2017E/MET/AOD/09Aug2019_UL2017_rsb-v1/060000/13E6F967-FC6F-BD4A-81C5-C4D1A610B8B8.root']
    config.Data.outputPrimaryDataset = 'store'
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 1
else:
    config.Data.inputDataset = '/MET/Run2017E-09Aug2019_UL2017_rsb-v1/AOD'    
    config.Data.inputDBS = 'global'
    #config.Data.splitting = 'Automatic'
    #config.Data.splitting = 'LumiBased'
    #config.Data.unitsPerJob = 1 #20
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 1

    #config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
    config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
    #config.Data.lumiMask = '/afs/cern.ch/work/d/dapparu/private/thesis/hscp/CMSSW_10_6_20/src/SUSYBSMAnalysis/crab_projects/crab_EDMProd_MET_UL2017B_297047_297484/results/processedLumis.json'

    config.Data.runRange = '297047-306462'

config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
config.Data.outLFNDirBase = '/store/user/dapparu/HSCP/Analysis/v2-1/prodJuly2021_CMSSW_10_6_20/'

config.Site.storageSite = 'T2_FR_IPHC'

if __name__=='__main__':

    from CRABAPI.RawCommand import crabCommand
    
    if test:
        crabCommand('submit', config = config)
        exit(0)

    datasetlist = ['/MET/Run2017E-09Aug2019_UL2017_rsb-v1/AOD']

    for dataset in datasetlist:
        config.Data.inputDataset = dataset
        Name = dataset.split('/')[1]+'_'+dataset.split('/')[2]
        config.General.requestName = Name
        crabCommand('submit', config = config)

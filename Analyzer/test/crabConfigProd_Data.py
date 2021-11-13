from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'EDMProd_SingleMuon_UL2017C_299368_300284_reRunJSON'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HSCParticleProducer_cfg.py'
config.JobType.allowUndistributedCMSSW = True

test=False
if test:
    config.Data.userInputFiles = ['root://cms-xrd-global.cern.ch//store/data/Run2017B/MET/AOD/09Aug2019_UL2017_rsb-v1/00000/AA1FC1E6-1E88-204D-B867-4637AEAC4BEA.root']
    config.Data.outputPrimaryDataset = 'store'
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 1
else:
    config.Data.inputDataset = '/SingleMuon/Run2017C-09Aug2019_UL2017-v1/AOD'
    config.Data.inputDBS = 'global'
    #config.Data.splitting = 'Automatic'
    #config.Data.splitting = 'LumiBased'
    #config.Data.unitsPerJob = 1 #20
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 1

    #config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
    #config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
    config.Data.lumiMask = '/afs/cern.ch/work/d/dapparu/private/thesis/hscp/CMSSW_10_6_20/src/SUSYBSMAnalysis/crab_projects/crab_Analysis_SingleMuon_UL2017C_299368_300284/results/processedLumis.json'
    config.Data.runRange = '299368-300284'
#    config.Data.runRange = '297047-306462'

config.Data.publication = True
config.Data.outputDatasetTag = config.General.requestName

config.Site.storageSite = 'T2_FR_IPHC'

if __name__=='__main__':

    from CRABAPI.RawCommand import crabCommand
    
    if test:
        crabCommand('submit', config = config)
        exit(0)

    datasetlist = ['/MET/Run2017B-09Aug2019_UL2017_rsb-v1/AOD']

    for dataset in datasetlist:
        config.Data.inputDataset = dataset
        Name = dataset.split('/')[1]+'_'+dataset.split('/')[2]
        config.General.requestName = Name
        crabCommand('submit', config = config)

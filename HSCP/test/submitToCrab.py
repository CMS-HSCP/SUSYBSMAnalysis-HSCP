#!/usr/bin/env python
import sys, argparse, os
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parser = argparse.ArgumentParser(prog='python '+sys.argv[0])
parser.add_argument('--dataset', dest='dataset', metavar='<Dataset>')
parser.add_argument('--inputFiles', dest='inputFiles', metavar='<file>', nargs='+')
parser.add_argument('--name', dest='requestName', metavar='<request-name>')
parser.add_argument('--sample', dest='sample', metavar='<isSignal or isBckg or isData>',choices=['isSignal', 'isBckg', 'isData'])
parser.add_argument('--lumiToProcess', dest='lumiToProcess', metavar='<JSON file>',default='')#certified_lumi_file
parser.add_argument('--dryrun', dest='dryrun', action='store_true')
parser.add_argument('--skip-estimates', dest='skipEstimates', action='store_true')
parser.add_argument('--njobs', dest='njobs', metavar='<unitsPerJob>')
args = parser.parse_args()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_copy(cfg_file,sample):
    import random,string
    rdm_str=''.join(random.choice(string.ascii_letters) for i in range(6))
    temp_cfg_file='tmp{}_{}'.format(rdm_str,cfg_file)
    with open(cfg_file) as fi: lines=fi.readlines()
    print('Creating temporally file "{}"'.format(temp_cfg_file))
    with open(temp_cfg_file,'w') as fo:
        for line in lines:
            if 'options.register' in line and 'SAMPLE' in line:
                line="options.register('SAMPLE', '{}',\n".format(sample)
            if sample=='isData' and 'options.register' in line and 'LUMITOPROCESS' in line:
                line="options.register('LUMITOPROCESS', '{}',\n".format(args.lumiToProcess)
            fo.write(line)
    return temp_cfg_file
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def delete_tmp_file(f):
    print('Deleting temporally file "{}"'.format(f))
    os.system('rm '+f)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def print_help(tmp_file):
    parser.print_help()
    delete_tmp_file(tmp_file)
    exit(0)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

try:
    from CRABClient.UserUtilities import config
except ImportError:
    exit('!!!You need to run: source /cvmfs/cms.cern.ch/common/crab-setup.sh')
config = config()

config.General.requestName = 'test'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
psetName=get_copy('HSCParticleProducer_cfg.py', args.sample)
config.JobType.psetName = psetName
config.JobType.allowUndistributedCMSSW = True

if args.inputFiles:
    config.Data.userInputFiles = args.inputFiles
    config.Data.outputPrimaryDataset = 'store'
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 1
elif args.dataset:
    config.Data.inputDataset = args.dataset
    config.Data.inputDBS = 'global'
    if args.sample=='isData':
        config.Data.splitting = 'FileBased'
        config.Data.unitsPerJob = 1
        config.Data.lumiMask = args.lumiToProcess #'HSCP/test/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
        #config.Data.runRange = '297047-306462'
    else:
        config.Data.splitting = 'FileBased'
        config.Data.unitsPerJob = 1
else:
    print_help(psetName)

config.Data.publication = True

config.Site.storageSite = 'T2_FR_IPHC'

if __name__=='__main__':

    from CRABAPI.RawCommand import crabCommand
    
    if not args.requestName:
        print_help(psetName)
    else:
        config.General.requestName = args.requestName
        config.Data.outputDatasetTag = 'CRAB3_'+args.requestName
        if args.njobs:
            config.Data.unitsPerJob = args.njobs
        ##
        if args.skipEstimates :
            crabCommand('submit', dryrun=True, config=config,**{"skip-estimates": True})
        elif args.dryrun and not args.skipEstimates:
            crabCommand('submit', dryrun=True, config=config)
        else:
            crabCommand('submit', config=config)
    ############################################
    delete_tmp_file(psetName)

# coding: utf-8
import ROOT as r
import subprocess

class hscpChain():
    def __init__(self, fileslist, treenamePath):
        self.fileslist    = fileslist
        self.treenamePath = treenamePath
        self.chain = self.__createTChain()

    def __createTChain(self):
        chain = r.TChain(self.treenamePath)
        for num, rootfile in enumerate(self.fileslist):
            chain.AddFile(rootfile)
        return chain

    def getRunAndEventID (self, masscut):
        print('+'*20)
        print(self.treenamePath)
        for event in self.chain:
            if event.Pt > 60 and event.I > 0.1:
                if event.Mass > masscut:
                    print('{0} - {1}'.format(event.Run, event.Event))
    
#listOfFiles = get_ipython().getoutput('ls /afs/cern.ch/user/j/jpriscia/work/CMSSW_8_0_30/src/SUSYBSMAnalysis/HSCP/test/AnalysisCode/Results/Type0/Histo*Data13TeV16*_Run2016_*.root')
result = subprocess.run(['ls /afs/cern.ch/user/j/jpriscia/work/CMSSW_8_0_30/src/SUSYBSMAnalysis/HSCP/test/AnalysisCode/Results/Type0/Histo*Data13TeV16_Run2016_*.root'], shell=True, stdout=subprocess.PIPE) 
result_runG = subprocess.run(['ls /afs/cern.ch/user/j/jpriscia/work/CMSSW_8_0_30/src/SUSYBSMAnalysis/HSCP/test/AnalysisCode/Results/Type0/Histo*Data13TeV16G_Run2016_*.root'], shell=True, stdout=subprocess.PIPE) 
listOfFiles   = result.stdout.decode().split()
listOfFiles_G = result_runG.stdout.decode().split()

chain      = hscpChain(listOfFiles, 'Data13TeV16/HscpCandidates')
chain_runG = hscpChain(listOfFiles_G, 'Data13TeV16G/HscpCandidates')

chain.getRunAndEventID(900)
chain_runG.getRunAndEventID(900)
            

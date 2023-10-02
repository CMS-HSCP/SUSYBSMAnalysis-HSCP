import ROOT
import numpy as np

def makeRegion(name, iHist, syst):
    '''
    takes in a 20x20 HSCP ProbQvsIas histogram and returns the 20x1 Pass and Fail histos
    '''
    # first do Pass
    hPassTempName = 'hpass{}_temp'.format(syst)
    hPassTemp = ROOT.TH2F(hPassTempName,'ProbQvsIas;ProbQ;Ias',1,0.0,1.0,20,0.0,1.0) # need to change x [0,0.1]?
    # loop over 
    for yBin in range(1,21):
        ySum = 0.0
        for xBin in range(18,21):
            print('Getting bin content: ({},{})'.format(xBin,yBin))
            ySum += iHist.GetBinContent(xBin,yBin)
        hPassTemp.SetBinContent(1,yBin,ySum)
#    hPassTemp.SetBinContent(1,21,10)
#    hPassTemp.SetBinContent(1,22,10)

    # now do Fail
    hFailTempName = 'hfail{}_temp'.format(syst)
    hFailTemp = ROOT.TH2F(hFailTempName,'ProbQvsIas;ProbQ;Ias',1,0.0,1.0,20,0.0,1.0) # need to change x [0.1,0.7]?
    for yBin in range(1,21):
        ySum = 0.0
        for xBin in range(3,18):
            print('Getting bin content: ({},{})'.format(xBin,yBin))
            ySum += iHist.GetBinContent(xBin,yBin)
        hFailTemp.SetBinContent(1,yBin,ySum)
#    hFailTemp.SetBinContent(1,21,10)
#    hFailTemp.SetBinContent(1,22,10)

    f = ROOT.TFile.Open('HSCP_{}.root'.format(name),'UPDATE')
    hPassName = 'hpass{}'.format(syst)
    hFailName = 'hfail{}'.format(syst)
    # 2DAlphabet expects signal/blinding to be on X axis only.. So need to invert
    hPass = ROOT.TH2F(hPassName,'IasVsProbQ;Ias;ProbQ',20,0.0,1.0,1,0.0,1.0) # need to change y [0.0, 0.1]?
    hFail = ROOT.TH2F(hFailName,'IasVsProbQ;Ias;ProbQ',20,0.0,1.0,1,0.0,1.0) # need to change y [0.1, 0.7]?
    for yBin in range(1,21):
      hPNew = hPassTemp.GetBinContent(1,yBin)
      hPass.SetBinContent(yBin,1,hPNew)
      hFNew = hFailTemp.GetBinContent(1,yBin)
      hFail.SetBinContent(yBin,1,hFNew)

    hFail.Write()
    hPass.Write()
    f.Close()
    
#f = ROOT.TFile.Open('crab_Analysis_2018_HSCPgluino_M-1800_CodeV40p4_v1.root')
#h = f.Get('HSCParticleAnalyzer/BaseName/PostPreS_ProbQNoL1VsIas')

d = {'Signal': 'crab_Analysis_2018_HSCPgluino_M-1800_CodeV40p9_v1.root',
 'Data': 'crab_Analysis_SingleMuon_Run2018C_CodeV40p9_v1.root'}

for name, fName in d.items():
    fTemp = ROOT.TFile.Open(fName)
    makeRegion(name, fTemp.Get('HSCParticleAnalyzer/BaseName/PostPreS_ProbQNoL1VsIas'),'')
    if("Signal" in name) :
      makeRegion(name, fTemp.Get('HSCParticleAnalyzer/BaseName/PostPreS_ProbQNoL1VsIas__Pileup_down'),'_PUsyst_down')
      makeRegion(name, fTemp.Get('HSCParticleAnalyzer/BaseName/PostPreS_ProbQNoL1VsIas__Pileup_up'),'_PUsyst_up')
    fTemp.Close()

print("scp HSCP_*root vami@ui3.kfki.hu:/data/vami/CMSSW_10_6_14/src/DataStudies/.")

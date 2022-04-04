#!/usr/bin/env python

import urllib
import string
import os
import sys
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor 
import glob
import ROOT
import math

#SStauMassPoints=[100, 126, 156, 200, 247, 308, 370, 432, 494, 557, 651, 745, 871, 1029, 1200, 1400, 1600, 1800, 2000, 2500]
SStauMassPoints=[100, 126, 156, 200, 247, 308, 370, 432, 494, 557, 651, 745, 871, 1029, 1200, 1600, 2500]
StauMassPoints=[100, 126, 156, 200, 247, 308, 370, 432, 494, 557, 651, 745, 871, 1029]
StauMassCut   =[10, 20, 50, 90, 130, 180, 230, 280, 320, 370, 440, 500, 570, 660]
DYMassPoints=[100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
DYMassCut   =[40, 120, 190, 270, 340, 400, 470, 530, 590, 650]
TeresaM=[100,200,300,400,500,600,700,800,900]
TeresaW=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
#TeresaDecayWidth = [9.0926593E-14, 4.7183150E-15, 1.0390462E-15, 4.7807750E-16, 6.7618337E-20, 9.0151758E-14, 4.6781076E-15, 1.0301919E-15, 4.7400354E-16, 6.7042124E-20, 8.8845352E-14, 4.6103163E-15, 1.0152632E-15, 4.6713466E-16, 6.6070604E-20, 8.6983425E-14, 4.5136982E-15, 9.9398637E-16, 4.5734495E-16, 6.4685966E-20, 8.4529278E-14, 4.3863489E-15, 9.6594207E-16, 4.4444144E-16, 6.2860919E-20, 8.1429382E-14, 4.2254907E-15, 9.3051860E-16, 4.2814268E-16, 6.0555655E-20, 7.7606394E-14, 4.0271102E-15, 8.8683215E-16, 4.0804202E-16, 5.7712657E-20]
#TeresaDecayWidth = [4.7807800E-16, 3.4004500E-16, 2.7204700E-16, 2.0405000E-16, 1.3605200E-16, 6.8055100E-17, 6.7618300E-20, 4.7400400E-16, 3.4004500E-16, 2.7204700E-16, 2.0405000E-16, 1.3605200E-16, 6.8055100E-17, 6.7042100E-20, 4.6713500E-16, 3.4004500E-16, 2.7204700E-16, 2.0405000E-16, 1.3605200E-16, 6.8055100E-17, 6.6070600E-20, 4.5734500E-16, 3.4004500E-16, 2.7204700E-16, 2.0405000E-16, 1.3605200E-16, 6.8055100E-17, 6.4686000E-20, 4.4444100E-16, 3.4004500E-16, 2.7204700E-16, 2.0405000E-16, 1.3605200E-16, 6.8055100E-17, 6.2860900E-20, 4.2814300E-16, 3.4004500E-16, 2.7204700E-16, 2.0405000E-16, 1.3605200E-16, 6.8055100E-17, 6.0555700E-20, 4.0804200E-16, 3.4004500E-16, 2.7204700E-16, 2.0405000E-16, 1.3605200E-16, 6.8055100E-17, 5.7712700E-20]
TeresaCTau = [0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 6.0, 8.0, 10., 15., 20., 30.]

pMSSMINC = ["pMSSM12_MCMC1_30_549144_sgsgToChipmsq_m200", "pMSSM12_MCMC1_30_549144_sgsgToChipmsq_m700", "pMSSM12_MCMC1_30_549144_sgsgToChipmtb_m200", "pMSSM12_MCMC1_30_549144_sgsgToChipmtb_m700", "pMSSM12_MCMC1_30_549144_sqsq_m200", "pMSSM12_MCMC1_30_549144_sqsq_m700"]


Energy = ["8TeV"]#, "7TeV"]

CMSSW_VERSION = os.getenv('CMSSW_VERSION','CMSSW_VERSION')
if CMSSW_VERSION == 'CMSSW_VERSION':
  print 'please setup your CMSSW environement'
  sys.exit(0)


if len(sys.argv)==1:
	print "Please pass in argument a number between 0 and 2"
	sys.exit()

elif sys.argv[1]=='0':	
        print 'Build Efficiency maps'
        FarmDirectory = "FARM"
        JobName = "HscpBuildEffMaps"
	LaunchOnCondor.Jobs_RunHere = 1
	LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)	
        for m in SStauMassPoints :
           LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/GetEfficiencyMaps.C", '"SStau'+str(m)+'"', "1", '"root://eoscms//eos/cms/store/user/querten/ModelIndepSept/ModelIndep_SingleStau'+str(m)+'.root"'])
	LaunchOnCondor.SendCluster_Submit()

elif sys.argv[1]=='1':
        print 'Merge efficiencies'
        fileList = ''
        for m in SStauMassPoints :
           fileList+=' pictures/Histos_SStau'+str(m)+'.root'        
        os.system('hadd -f pictures/Histos.root' + fileList)
        os.system('root MakePlot.C++ -l -b -q');


elif sys.argv[1]=='2':
        print 'Compute model independent acceptance'
        FarmDirectory = "FARM"
        JobName = "HscpModelIndepAccep"
        LaunchOnCondor.Jobs_RunHere = 1
        LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
        for E in Energy:           
           Suffix=''
           if(E=="7TeV"):Suffix="_BX1";
           for m in StauMassPoints :
              LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/ModelIndependent_Acceptance.C", '"MI_PPStau_'+E+'_M'+str(m)+'"', '"root://eoscms//eos/cms//store/cmst3/user/querten/12_08_30_HSCP_EDMFiles/PPStau_'+E+'_M'+str(m)+Suffix+'.root"'])
              LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/ModelIndependent_Acceptance.C", '"MI_GMStau_'+E+'_M'+str(m)+'"', '"root://eoscms//eos/cms//store/cmst3/user/querten/12_08_30_HSCP_EDMFiles/GMStau_'+E+'_M'+str(m)+Suffix+'.root"'])
           for m in DYMassPoints :
              LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/ModelIndependent_Acceptance.C", '"MI_DY_'+E+'_M'+str(m)+'"', '"root://eoscms//eos/cms//store/cmst3/user/querten/12_08_30_HSCP_EDMFiles/HSCP_'+E+'_M'+str(m)+Suffix+'_Q1.root"'])
              LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/ModelIndependent_Acceptance.C", '"MI_pMSSM_'+E+'_M'+str(m)+'"', '"root://eoscms//eos/cms//store/cmst3/user/querten/12_08_30_HSCP_EDMFiles/pMSSM_'+E+'_ChipmChipm_mChi'+str(m)+Suffix+'.root"'])
              LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/ModelIndependent_Acceptance.C", '"MI_pMSSM_MassOnly_'+E+'_M'+str(m)+'"', '"root://eoscms//eos/cms//store/cmst3/user/querten/12_08_30_HSCP_EDMFiles/pMSSM_MassOnly_'+E+'_ChipmChipm_mChi'+str(m)+Suffix+'.root"'])
           I=0;
           for m in TeresaM:
              for w in TeresaW:
                LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/ModelIndependent_Acceptance.C", '"MI_pMSSM12_MCMC1_30_549144_m'+str(m)+'_width'+str(w)+'"', '"root://eoscms//eos/cms//store/cmst3/user/querten/12_08_30_HSCP_EDMFiles/pMSSM12_MCMC1_30_549144_new_m'+str(m)+'_width'+str(w)+Suffix+'.root"', TeresaCTau[w]])
                I+=1;

        for n in pMSSMINC:
           LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/ModelIndependent_Acceptance.C", '"MI_'+str(n)+'"', '"root://eoscms//eos/cms//store/cmst3/user/querten/12_08_30_HSCP_EDMFiles/'+str(n)+'.root"'])

        LaunchOnCondor.SendCluster_Submit()
elif sys.argv[1]=='3':
        print 'Compute acceptance for standard analysis'
        FarmDirectory = "FARM"
        JobName = "HscpStandardAccep"
        LaunchOnCondor.Jobs_RunHere = 1
        LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
        for E in Energy:
           Suffix=''
           if(E=="7TeV"):Suffix="_BX1";
           I=0;
           for m in StauMassPoints :
#              LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/StandardAnalysis_Acceptance.C", '"PPStau_'+E+'_M'+str(m)+Suffix+'"', '1', str(StauMassCut[I])])
              I+=1;
           I=0;
           for m in StauMassPoints :
#              LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/StandardAnalysis_Acceptance.C", '"GMStau_'+E+'_M'+str(m)+Suffix+'"', '1', str(StauMassCut[I])])
              I+=1;
           I=0;
           for m in DYMassPoints :
#              LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/StandardAnalysis_Acceptance.C", '"DY_'+E+'_M'+str(m)+'_Q1'+Suffix+'"', '1', str(DYMassCut[I])])
              I+=1;
           I=0;
           for m in DYMassPoints :
#              LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/StandardAnalysis_Acceptance.C", '"pMSSM_'+E+'_M'+str(m)+Suffix+'"', '1', str(DYMassCut[I])])
#              LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/StandardAnalysis_Acceptance.C", '"pMSSM_MassOnly_'+E+'_M'+str(m)+Suffix+'"', '1', str(DYMassCut[I])])
              I+=1;
           for m in TeresaM:
              for w in TeresaW:
#                 LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/StandardAnalysis_Acceptance.C", '"TeresapMSSM12_MCMC1_30_549144_m'+str(m)+'_width'+str(w)+Suffix+'"', '1', "100"])
                 I+=1;

        for n in pMSSMINC:
           LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/StandardAnalysis_Acceptance.C", '"'+str(n)+'"', '1', "100"])

        LaunchOnCondor.SendCluster_Submit()



elif sys.argv[1]=='4':
        print 'Compute acceptance for all pMSSM points'
        FarmDirectory = "FARM"
        JobName = "HscpModelIndepAccepPMSSM"
        LaunchOnCondor.Jobs_RunHere = 1
        LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)

        Rfile = ROOT.TFile("pMSSM/pMSSM12_56K_full.root")
        Rtree = Rfile.Get('MCMC') 
        i = 0
        for Revent in Rtree:     
#           if(Revent.ctau>1000):continue #only keep >1m ctau
           if(Revent.ctau<500):continue #only keep >1m ctau
           m = Revent.mw1
           i = i+1;   
#           if(i>100):continue
           LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/ModelIndependent_Acceptance.C", '"pMSSM_'+str(Revent.itr)+'"', '"file:/afs/cern.ch/user/q/querten/workspace/public/13_06_01_2012_HSCP_Analysis_Archive/SingleHSCP_Gen/GEN-SIM/CMSSW_5_2_9/src/pMSSM_Scan/InFiles/%i/out.root"' % Revent.itr, str(Revent.ctau/1000.0)])
        LaunchOnCondor.SendCluster_Submit()


elif sys.argv[1]=='5':
        print 'Compute acceptance for all pMSSM points'
        FarmDirectory = "FARM"
        JobName = "HscpModelIndepAccepPMSSM"
        LaunchOnCondor.Jobs_RunHere = 1
        LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)

        Rfile = ROOT.TFile("pMSSM/pMSSM12_56K_full.root")
        Rtree = Rfile.Get('MCMC') 
        i = 0
        for Revent in Rtree:     
           if(Revent.ctau<500):continue #only keep >1m ctau
           m = Revent.mw1
           i = i+1;   
           LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/ModelIndependent_Acceptance.C", '"pMSSMInc_'+str(Revent.itr)+'"', '"file:/afs/cern.ch/user/q/querten/workspace/public/13_06_01_2012_HSCP_Analysis_Archive/SingleHSCP_Gen/GEN-SIM/CMSSW_5_2_9/src/pMSSM_ScanInc/InFiles/%i/out.root"' % Revent.itr, str(Revent.ctau/1000.0)])
        LaunchOnCondor.SendCluster_Submit()


elif sys.argv[1]=='9':
        Rfile = ROOT.TFile("pMSSM/pMSSM12_56K_full.root")
        Rtree = Rfile.Get('MCMC')
        i = 0
        j = 0
        k = 0
        for Revent in Rtree:
           i = i+1
           if(math.fabs(Revent.mh-125)<2.5): j= j+1
           if(Revent.ctau>10 and Revent.ctau<500): k = k+1
        print i
        print j
        print k

else:
	print 'Unknown case: use an other argument or no argument to get help'




#!/usr/bin/env python

import string
import os
import sys
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor  

SAMPLES = ['MU', 'MET', 'ELEC', 'Gluino300', 'Gluino1000', 'Gluino1500', 'GMStau126', 'GMStau494']
PT   = [20.0, 40.0, 60.0, 100]
NH   = [3, 5, 7]
DEDX = [2.7, 3.0, 3.3]

if sys.argv[1]=='1':
   JobName = "SkimEff"
   FarmDirectory = "FARM_SkimEff"
   LaunchOnCondor.Jobs_NEvent = 100000
   LaunchOnCondor.Jobs_Skip = 0
   LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
   for s in SAMPLES:
      for p in PT:
         for n in NH:
            for d in DEDX:
              LaunchOnCondor.SendCluster_Push  (["CMSSW", "HSCPSkim_cfg.py", "XXX_SAMPLE_XXX", str(s), "XXX_PT_XXX", str(p), "XXX_NH_XXX", str(n), "XXX_DEDX_XXX", str(d)])
   LaunchOnCondor.SendCluster_Submit()

if sys.argv[1]=='2':
   #LaunchOnCondor.runInteractively = True
   LaunchOnCondor.LSFlog = False
   JobName = "Analyze"
   FarmDirectory = "FARM_Analyze"
   LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
   for s in SAMPLES:
      for p in PT:
         for n in NH:
            for d in DEDX:
              LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/analyzeEDM.C", '"SkimEff_%s_pT%s_nH%s_dEdx%s"' % (str(s),  str(p), str(n),  str(d)), '"%s"' % os.getcwd()])
   LaunchOnCondor.SendCluster_Submit()

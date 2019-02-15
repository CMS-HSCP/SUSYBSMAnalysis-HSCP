#! /usr/bin/env python

import os, time, string

indir = "Results/"
jobfilesdir = "FARM/inputs/" 
samplesFiles = "Analysis_Samples.txt"
writeCMDFile = True
types = set()

#for root, dirs, files in os.walk(indir):
#    if "Type" not in root: continue
#    anaType = int(root.split("/")[-1].replace("Type",""))
#    types.add(anaType)

def initProxy():
   if(not os.path.isfile(os.path.expanduser('~/x509_user_proxy/x509_proxy')) or ((time.time() - os.path.getmtime(os.path.expanduser('~/x509_user_proxy/x509_proxy')))>600)):
      print "You are going to run on a sample over grid using either CRAB or the AAA protocol, it is therefore needed to initialize your grid certificate"
      os.system('mkdir -p ~/x509_user_proxy; voms-proxy-init --voms cms -valid 192:00 --out ~/x509_user_proxy/x509_proxy')#all must be done in the same command to avoid environement problems.  Note that the first sourcing is only needed in Louvain

types.add(0)
types.add(2)

todo = []
if writeCMDFile:
   with open(samplesFiles) as ifile:
       for iline, l in enumerate(ifile):
           line = l.strip()
           if len(line)==0 or line[0]=='#' : continue
           spl = [l.strip().strip('"') for l in line.split(",")]
           if int(spl[1])>0: continue
           expectedFileName= "Histos_{}_{}.root".format(spl[2], spl[3])
           sampleString = "ANALYSE_{}_to_{}".format(iline, iline)
           for t in types:
               fp = indir+"Type{}".format(t)+"/"+expectedFileName
###            if not os.path.isfile(fp) or os.path.getsize(fp)<1024:
               if (os.path.isfile(fp) and os.path.getsize(fp)<1024) or not os.path.isfile(fp):
                   typeString = ", {},".format(t)
                   for root, _, files in os.walk(jobfilesdir):
                       for f in files:
                           if "_HscpAnalysis.sh" not in f: continue
                           contents = open(os.path.join(root, f)).read()
                           if sampleString not in contents or typeString not in contents:
                               continue
                           todo.append(f)
                           print f

   newcmd = open("newcmd.cmd", "w")
#   newcmd.write('#!/bin/bash\n')
   newcmd.write('Universe                = vanilla\nEnvironment             = CONDORJOBID=$(Process)\nnotification            = Error\nrequirements            = (Memory > 200)\nshould_transfer_files   = YES\nwhen_to_transfer_output = ON_EXIT\n\n')


   for i, executable in enumerate(todo):
   #   if os.environ['HOST'].find('lxplus')>=0:
   #      toWrite = 'bsub -q 8nh -J HSCPResub_%i \'%s%s 0 ele\'' % (i, os.environ['CMSSW_BASE'], executable)
   #   elif os.environ['HOST'].find('ingrid'):
   #      toWrite = 'sbatch --partition=Def --qos=normal --wckey=cms %s\n' % executable
#      toWrite = 'sbatch --partition=Def --qos=normal --wckey=cms FARM/inputs/%s' % executable
#      toWrite = 'bsub -q 8nh -R "type==SLC6_64&&pool>30000" -J HscpAnalysis'+string.replace(executable,"HscpAnalysis.sh","")+' -oo /afs/cern.ch/work/j/jozobec/public/Latest/CMSSW_8_0_30/src/SUSYBSMAnalysis/HSCP/test/AnalysisCode/FARM/outputs/'+string.replace(executable,"sh","cout")+' \'/afs/cern.ch/work/j/jozobec/public/Latest/CMSSW_8_0_30/src/SUSYBSMAnalysis/HSCP/test/AnalysisCode/FARM/inputs/%s\'' % (executable)
      toWrite = 'Executable              = FARM/inputs/%s\n' % executable
      toWrite+= 'output                  = FARM/logs/%s\n' % string.replace(executable, "sh", "cout")
      toWrite+= 'error                   = FARM/logs/%s\n' % string.replace(executable, "sh", "cerr")

      toWrite+= 'log                     = FARM/logs/%s\n' % string.replace(executable, "sh", "log")

      toWrite+='Queue 1\n\n'


      print toWrite
      newcmd.write(toWrite+'\n')
   newcmd.close()

initProxy()
os.system('condor_submit newcmd.cmd')


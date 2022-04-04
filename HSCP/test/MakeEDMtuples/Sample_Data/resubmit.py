#/bin/env/python

import os, sys, string, json
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor

#
# option 1 just lists runs that have to be rerun
# option 2 creates a new .cmd file that can be used to resubmit on condor
# for extra functionality consult options
#
# modes for 1st option
#################################
printOnlySuccess      = True    # list only runs ready for merging
printNumberOfFiles    = True    # list number of good MET/SingleMuon/DoubleMuon files vs number of expected files
mergeAvailableRuns    = True    # if we print only success, merge the available runs
mergeOnTheFly         = True    # make sure not to resubmit runs that are already being merged
overwriteMerged       = True    # proceed with merging, even if the file exists
# modes for 2nd option
#################################
resubmit              = True    # resubmit on the cluster at the end
resubmitOnTheFly      = True    # creates .cmd file in such a way, that the already running jobs are not resubmitted
checkFileSize         = True    # checks not only if the target file exists, but also its size in determining whether or not the job failed
checkLumiContent      = False   # checks not only the target file's size, but also its lumi content in determining whether or not the job failed
ensureNecessaryRuns   = True    # skips jobs in FARM/inputs that are obsolete
EmptyFileSize         = 2321123 # empty nTuple file size (no lumi content whatsoever)
resubmitScriptName    = "resubmit.cmd" if not resubmitOnTheFly else "resubmitOnTheFly.cmd"

def initProxy():
   print "You are going to run on a sample over grid using either CRAB or the AAA protocol, it is therefore needed to initialize your grid certificate"
   os.system('mkdir -p ~/x509_user_proxy; voms-proxy-init --voms cms -valid 192:00 --out ~/x509_user_proxy/x509_proxy')#all must be done in the same command to avoid environement problems.  Note that the first sourcing is only needed in Louvain

if len(sys.argv)==1 or sys.argv[1]=="1":
   runs=os.listdir("out")
   runsToPrint = []
   finalJudgementString = "pending." if not printOnlySuccess else "ready for merging!"
   for run in runs:
      if not os.path.isdir("out/"+run):
         continue
      doubleMuon=0
      singleMuon=0
      MET       =0
      correctMET       =int(os.popen("grep %s FARM/inputs/*.sh | grep MET        | wc -l" % run).read())
      correctSingleMuon=int(os.popen("grep %s FARM/inputs/*.sh | grep SingleMuon | wc -l" % run).read())
      correctDoubleMuon=int(os.popen("grep %s FARM/inputs/*.sh | grep DoubleMuon | wc -l" % run).read())
      listOfFiles = os.listdir("out/%s" % run)
      for file in listOfFiles:
         fileStats=os.stat("out/%s/%s" % (run, file))
         if file.find("SingleMuon")!=-1 and fileStats.st_size>EmptyFileSize:
            singleMuon += 1
            continue
         elif file.find("DoubleMuon")!=-1 and fileStats.st_size>EmptyFileSize:
            doubleMuon += 1
            continue
         elif file.find("MET")!=-1 and fileStats.st_size>EmptyFileSize:
            MET += 1
            continue
         elif file.find("MET")==-1 and file.find("SingleMuon")==-1 and file.find("DoubleMuon")==-1:
            print "Some weird file was found ..."
            continue

      if ((correctMET != MET) or (correctSingleMuon != singleMuon) or (correctDoubleMuon != doubleMuon)) and not printOnlySuccess:
         runsToPrint.append([run, MET, correctMET, singleMuon, correctSingleMuon, doubleMuon, correctDoubleMuon])
      elif correctMET==MET and correctSingleMuon==singleMuon and correctDoubleMuon==doubleMuon and printOnlySuccess:
         runsToPrint.append([run, MET, correctMET, singleMuon, correctSingleMuon, doubleMuon, correctDoubleMuon])
   if not printNumberOfFiles:
      for entry in runsToPrint:
         print entry[0], finalJudgementString
   else:
      for entry in runsToPrint:
         print "Run %s %s: MET (%i/%i), SingleMuon (%i/%i), DoubleMuon (%i/%i)" % (entry[0], finalJudgementString, entry[1], entry[2], entry[3], entry[4], entry[5], entry[6])

   if printOnlySuccess and mergeAvailableRuns:

      FarmDirectory = 'MERGEAvailable'
      LaunchOnCondor.SendCluster_Create(FarmDirectory, "HSCPEdmMergeAvail")
      LaunchOnCondor.Jobs_Queue = '8nh'
      AlreadyMergingRuns = []

      if mergeOnTheFly:
         runningScripts=os.popen('condor_q %s -long | grep Cmd | grep HSCPEdmMergeAvail' % (os.environ['USER'])).read()
         runningScripts=runningScripts.split('\n')
         runningScripts.pop()
         for script in runningScripts:
            if script.find(os.getcwd())==-1: continue
            script=script.split('=')[1]
            script=script.replace('"', '')
            script=script.split('/')[len(script.split('/'))-1]
            script=script.split('_')[2]
            script=script.split('.')[0]
            AlreadyMergingRuns.append(int(script))

      for entry in runsToPrint:
         RUN = int(entry[0])
         if mergeOnTheFly and RUN in AlreadyMergingRuns: continue
         if not overwriteMerged and os.path.isfile('out/Run2016_%i.root' % RUN): continue
         LaunchOnCondor.Jobs_InitCmds   = ['export HOME=%s' % os.environ['HOME']]
         LaunchOnCondor.Jobs_FinalCmds  = ["edmLumisInFiles.py Run2016_%i.root --output=%s/out/Run2016_%i.json" % (RUN, os.getcwd(), RUN)]
         LaunchOnCondor.Jobs_FinalCmds += ["mv Run2016_%i.root %s/out/Run2016_%i.root" % (RUN, os.getcwd(), RUN)]
         LaunchOnCondor.ListToFile(LaunchOnCondor.GetListOfFiles('"file:','%s/out/%i/*_HSCP_*.root' % (os.getcwd(), RUN),'",'), FarmDirectory + "InputFile.txt")
         LaunchOnCondor.SendCMSJobs(FarmDirectory, "HSCPEdmMergeAvail_%i"%RUN, "Merge_cfg.py", FarmDirectory + "InputFile.txt", 1, ['XXX_SAVEPATH_XXX','Run2016_%i.root' % RUN])
         os.system("rm " +  FarmDirectory + "InputFile.txt")

if sys.argv[1]=="2":
   shellJobs = os.listdir("FARM/inputs")
   shellJobs.sort()
   jobsToRerun=[]
   runsQuotedInJSON=[]
   runsQuotedInFARM=[]
   runsToRemove    =[]
   for file in shellJobs:
      fileIsGood=False
      if not file.find(".sh")!=-1: continue
      endPath=os.popen('grep HSCP_ FARM/inputs/%s' % file).read().split()[2]

      if ensureNecessaryRuns:
         run = endPath.split('/')[len(endPath.split('/'))-2]
         if not run in runsQuotedInFARM:
            runsQuotedInFARM.append(run)

      if os.path.isfile(endPath):
         if not checkFileSize: 
            fileIsGood = True
            continue

         fileStats=os.stat(endPath)
         if int(fileStats.st_size) > EmptyFileSize:
            if not checkLumiContent:
               fileIsGood = True
               continue
       
            runs = []
            os.system('edmLumisInFiles.py %s --out tmp.json' % endPath)
            runList = json.load('tmp.json').items()
            for entry in runList:
               runs.append("%s" % str(entry[0]))
            os.system('rm tmp.json')
            if not len(runs)<1: fileIsGood = True
	    
      if not fileIsGood:
         jobsToRerun.append(file)

   if ensureNecessaryRuns:
      print "\n"
      files = os.listdir("out/")
      files.sort()
      for file in files:
         if os.path.isdir("out/%s" % file):
            runsQuotedInJSON.append(file)

      for run in runsQuotedInFARM:
         if run in runsQuotedInJSON: continue
         else: runsToRemove.append(run)
      
      for run in runsToRemove: print run
      print "\n"
      for scriptFile in jobsToRerun:
         for run in runsToRemove:
            if int(os.popen('grep %s FARM/inputs/%s | wc -l' % (run, scriptFile)).read())==1:
               jobsToRerun.remove(scriptFile)

   if resubmitOnTheFly:
      runningJobs = os.popen('condor_q %s -long | grep Cmd' % os.environ['USER']).read()
      runningJobs = runningJobs.split('\n')
      runningJobs.pop()
      for job in runningJobs:
         if job.find(os.getcwd())==-1: continue
         runningScript = job.split('/')[len(job.split('/'))-1]
         runningScript = runningScript.replace('"', '')
         if runningScript in jobsToRerun:
            jobsToRerun.remove(runningScript)

   print "%i jobs to rerun" % len(jobsToRerun)
   if len(jobsToRerun) > 0:
      print "creating resubmit command script:"
      f = open(resubmitScriptName, "w")
      f.write('Universe                = vanilla\n')
      f.write('Environment             = CONDORJOBID=$(Process)\n')
      f.write('notification            = Error\n')
      f.write('requirements            = (CMSFARM=?=True)&&(Memory > 200)\n')
      f.write('should_transfer_files   = YES\n')
      f.write('when_to_transfer_output = ON_EXIT\n')
      for job in jobsToRerun:
         job = job.split('.')[0]
         f.write('\n')
         f.write('Executable              = FARM/inputs/%s.sh\n' % job)
         f.write('output                  = FARM/logs/%s.out\n'  % job)
         f.write('error                   = FARM/logs/%s.err\n'  % job)
         f.write('log                     = FARM/logs/%s.log\n'  % job)
         f.write('Queue 1\n')
      f.close()
      print "%s written" % resubmitScriptName

      if resubmit:
         initProxy()
         os.system('condor_submit %s' % resubmitScriptName)


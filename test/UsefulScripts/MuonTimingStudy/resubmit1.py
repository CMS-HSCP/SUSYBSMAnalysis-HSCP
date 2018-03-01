#!/bin/env python

import os, string, sys, json, os.path
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor

def initProxy():
   print "You are going to run on a sample over grid using either CRAB or the AAA protocol, it is therefore needed to initialize your grid certificate"
   os.system('mkdir -p ~/x509_user_proxy; voms-proxy-init --voms cms -valid 192:00 --out ~/x509_user_proxy/x509_proxy')#all must be done in the same command to avoid environement problems.  Note that the first sourcing is only needed in Louvain

def ListRootFiles():
	return os.popen('find out/*.root | xargs -I% basename %').read().split()

def ListShellFiles():
	return os.popen('find FARM_TEMP/inputs/*.sh ').read().split()

def ListFailedJobs():
	ShellFiles  = ListShellFiles()
	files = ListRootFiles()
	toReprocess = []
	for f in files:
		print '%i/%i/%i' % (len(toReprocess), len(ShellFiles), len(files))
		for script in ShellFiles:
			success = True if os.stat('out/%s' % f).st_size > 323 else False
			if f in open(script, 'r').read():
				ShellFiles.remove(script)
				if not success: toReprocess.append(script)
	return toReprocess + ShellFiles

def RemoveSuccessfulJobs():
	files = os.popen('find out | xargs -I% basename %').read().split()

if sys.argv[1] == '1':
	f = open('reprocess.cmd', 'w')
	f.write('#!/bin/bash\n\n')
	f.flush()
	jobs = ListFailedJobs()
	for job in jobs:
		print job + '\n'
		f.write('sbatch --partition=Def --qos=normal --wckey=cms %s\n' % job)
		f.flush()
	f.close()

if sys.argv[1] == '2':
	if not os.path.exists('reprocess.cmd'):
		os.system('python resubmit1.py 1')
	
	initProxy()
	os.system('sh reprocess.cmd')

#!/bin/env python

import os, string, sys, json, os.path
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor

def initProxy():
   print "You are going to run on a sample over grid using either CRAB or the AAA protocol, it is therefore needed to initialize your grid certificate"
   os.system('mkdir -p ~/x509_user_proxy; voms-proxy-init --voms cms -valid 192:00 --out ~/x509_user_proxy/x509_proxy')#all must be done in the same command to avoid environement problems.  Note that the first sourcing is only needed in Louvain

def ListShellFiles():
	return os.popen('find FARM/inputs/*.sh ').read().split()

def JobsToFilesMap():
	scripts = ListShellFiles()
	theMap  = []
	for script in scripts:
		files_index = os.popen('grep \'sh /home\' %s | awk \'{print $4}\'' % script).read().split()
		theMap.append ([script, 'pictures/Histos_%s.root' % files_index[0]])
	return theMap


def ListFailedJobs():
	theMap = JobsToFilesMap()
	toReprocess = []
	for i in range(0, len(theMap)):
		print '%i/%i' % (i+1, len(theMap))
		if not os.path.isfile(theMap[i][1]):
			toReprocess.append(theMap[i][0])
		elif os.path.isfile(theMap[i][1]) and os.stat(theMap[i][1]).st_size < 1024:
			toReprocess.append(theMap[i][0])
	return toReprocess

if sys.argv[1] == '1':
	f = open('reprocess.cmd', 'w')
	f.write('#!/bin/bash\n\n')
	f.flush()
	jobs = ListFailedJobs()
	for job in jobs:
		f.write('sbatch --partition=cp3 --qos=cp3 --wckey=cms %s\n' % job)
		f.flush()
	f.close()

if sys.argv[1] == '2':
	if not os.path.exists('reprocess.cmd'):
		os.system('python resubmit1.py 1')
	
	initProxy()
	os.system('sh reprocess.cmd')

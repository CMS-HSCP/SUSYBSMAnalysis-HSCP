#!/bin/env python

import os, string, sys, json
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor


def ListFailedFiles ():
	return os.popen('find out/*.root -size -1024c | xargs -I% basename %').read().split()

def ListShellFiles ():
	return os.popen('find FARM_TEMP/inputs/*.sh ').read().split()

def ListFailedJobs ():
	ShellFiles  = ListShellFiles()
	FailedFiles = ListFailedFiles()
	toReprocess = []
	for element in FailedFiles:
		print '%i/%i/%i' % (len(toReprocess), len(ShellFiles), len(FailedFiles))
		for script in ShellFiles:
			if element in open(script, 'r').read():
				toReprocess.append(script)
				ShellFiles.remove(script)
	return toReprocess + ShellFiles # scripts with no matches also have to be included

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


# original author: Joze Zobec
# joze.zobec@cern.ch
#############################

#!/bin/env python

import os
import string
import hashlib

NEW_VERSION = '8_0_30'
NEW_BRANCH  = 'Run2HSCP16_v4'
CMSSW_Old   = '/home/ucl/cp3/jzobec/Run2HSCP16/CMSSW_8_0_24_patch1/src'
TMP_DIR     = '/home/ucl/cp3/jzobec/Run2HSCP16/ToMigrate'
CMSSW_New   = TMP_DIR + '/CMSSW_' + NEW_VERSION + '/src'
OLDLIST     = 'V3.txt'
NEWLIST     = 'VOld.txt'
DIFF        = 'DIFF.patch'
CHECKSUMS   = 'sha256_sums.txt'

cdAndCmsenv = 'cd %s && eval `scramv1 runtime -sh`' % CMSSW_New

def execAndPrint (command):
   print command
   os.system(command)

def sha256(fname):
   method = hashlib.sha256()
   with open(fname, 'rb') as f:
      for chunk in iter(lambda: f.read(4096), b''):
         method.update(chunk)
   return method.hexdigest()

if not os.path.isdir(CMSSW_New+'/SUSYBSMAnalysis'):
   execAndPrint( 'mkdir -p %s && cd %s && scramv1 project CMSSW_%s && cd CMSSW_%s/src && eval `scramv1 runtime -sh` && git cms-init && git cms-addpkg SUSYBSMAnalysis && git checkout -b %s' % (TMP_DIR, TMP_DIR, NEW_VERSION, NEW_VERSION, NEW_BRANCH))
execAndPrint( 'find %s/SUSYBSMAnalysis > %s'       % (CMSSW_Old, OLDLIST)       )
execAndPrint( 'find %s/SUSYBSMAnalysis > %s'       % (CMSSW_New, NEWLIST)       )

# remove unnecessary lines
toRemove = ['PreHIP', 'PostHIP', 'Results', 'pictures', '\.root', 'FARM', '\.pyc', '_C\.d', '_C\.so','\.cmd', '__init__', 'ACLiC_dict_rdict\.pcm', '\.log', '\.list', '\.err', 'Harmonic', 'Failed', 'Type', '_bak', 'backup', 'Update', 'List']
for mask in toRemove:
   execAndPrint( 'sed -i \'/%s/d\' %s' % (mask, OLDLIST))
   execAndPrint( 'sed -i \'/%s/d\' %s' % (mask, NEWLIST))

# prepare file paths to be from the same root
execAndPrint( 'cat %s | cut -c %i- > DIFF1.txt && mv DIFF1.txt %s' % (NEWLIST, os.popen('find %s/SUSYBSMAnalysis -maxdepth 0' % CMSSW_New).read().find('SUSYBSMAnalysis'), NEWLIST))
execAndPrint( 'cat %s | cut -c %i- > DIFF2.txt && mv DIFF2.txt %s' % (OLDLIST, os.popen('find %s/SUSYBSMAnalysis -maxdepth 0' % CMSSW_Old).read().find('SUSYBSMAnalysis'), OLDLIST))
execAndPrint( 'diff -u %s %s > %s'% (NEWLIST, OLDLIST, DIFF))

# you just have to jack off sometimes
with open(DIFF, 'r') as f, open(CHECKSUMS, 'w') as r:
   dirsToRemove = []
   for i, line in enumerate(f):
      name = line[2:].strip()
      if len(line) == 0 or line.find('+++')>=0 or line.find('---')>=0 or line[0]==' ': continue
      if line[0] == '-':
         if os.path.isdir('%s/%s' % (CMSSW_New, name)):
            dirsToRemove.append('%s/%s' % (CMSSW_New, name))
         elif os.path.isfile('%s/%s' % (CMSSW_New, name)):
            execAndPrint('%s && git rm -f %s/%s' % (cdAndCmsenv, CMSSW_New, name))
      if line[0] == '+':
         if os.path.isdir('%s/%s' % (CMSSW_Old, name)):
            execAndPrint('mkdir -p %s/%s' % (CMSSW_New, name))
	 elif os.path.isfile('%s/%s' % (CMSSW_Old, name)) and (
              not name.find('.txt')>=0 or name.find('Analysis_Cuts.txt')>=0 or name.find('Analysis_Samples.txt')>=0
              or name.find('/data/')>=0 or name.find('README')>=0 or name.find('MakeEDMtuples')>=0
            ):
            r.write(name+'\n')
            execAndPrint('cp %s/%s %s/%s && %s && git add %s/%s' % (CMSSW_Old, name, CMSSW_New, name, cdAndCmsenv, CMSSW_New, name))
         else:
            print 'Can\'t find %s/%s' % (CMSSW_Old, name)

for entry in dirsToRemove:
   execAndPrint('rm -rf %s' % entry) # cleanup empty directories

# loop across all the files and compare their checksums
checkSums = []
print '#########################\n# Testing the checksums #\n#########################'
execAndPrint( 'find %s/SUSYBSMAnalysis > %s' % (CMSSW_New, NEWLIST))
execAndPrint( 'cat %s | cut -c %i- > DIFF1.txt && mv DIFF1.txt %s' %
  (NEWLIST, os.popen('find %s/SUSYBSMAnalysis -maxdepth 0' % CMSSW_New).read().find('SUSYBSMAnalysis'), NEWLIST))
with open(NEWLIST, 'r') as f, open(CHECKSUMS, 'a') as r:
   for line in f:
      name = line[1:].strip()
      if os.path.isfile('%s/%s' % (CMSSW_New, name)) and os.path.isfile('%s/%s' % (CMSSW_Old, name)):
         sum_old = sha256('%s/%s' % (CMSSW_New, name)) 
         sum_new = sha256('%s/%s' % (CMSSW_Old, name)) 

         checkSums.append ([name, sum_old, sum_new])

         # the new file should be updated to match the old one
         if sum_old != sum_new:
	    print 'Mismatch in %s:\t%s --> %s'     % (os.path.basename(name), sum_old, sum_new)
            execAndPrint('cp %s/%s %s/%s && %s && git add %s/%s' % (CMSSW_Old, name, CMSSW_New, name, cdAndCmsenv, CMSSW_New, name))
            r.write(name+'\n')

#execAndPrint('rm -rf %s %s %s' % (NEWLIST, OLDLIST, DIFF))


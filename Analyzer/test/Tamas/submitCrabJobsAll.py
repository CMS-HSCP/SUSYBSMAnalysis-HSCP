import sys, os, time, re
import numpy as np
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog codeVersion")
(opt,args) = parser.parse_args()


codeVersion = sys.argv[1]
#just the number, like 18p2

didVoms = input("Push enter if you alread did voms-proxy-init -rfc -voms cms -valid 192:00 otherwise say no and do it\n")
if(didVoms):
 sys.exit()

os.system("python3 submitCrabJobsData.py "+str(codeVersion))
os.system("python3 submitCrabJobsBackgroundWJetsOnly.py "+str(codeVersion))
os.system("python3 submitCrabJobsBackgroundQCDOnly.py "+str(codeVersion))
os.system("python3 submitCrabJobsBackgroundTTbarOnly.py "+str(codeVersion))
os.system("python3 submitCrabJobsSignalGluinoOnly.py "+str(codeVersion))





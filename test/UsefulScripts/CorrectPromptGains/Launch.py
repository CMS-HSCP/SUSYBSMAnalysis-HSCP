#!/bin/env/python

import string
import os,sys,time
import collections # kind of map

gains        = "/afs/cern.ch/cms/tracker/sistrvalidation/WWW/CalibrationValidation/ParticleGain"
#whitelist    = [272760, 273592, 274094, 274198, 274387, 274421, 275282, 275370, 275376, 276097, 276244, 276542, 276585, 276659]
#blacklist    = [273402, 273592]
whitelist    = [278406, 278509, 278769, 279715, 279966, 280242, 281975, 283043, 285090, 285368]
blacklist    = [283877, 285090, 285216]
useWhitelist = True
useBlacklist = False

def GetAppliedIntervals ():
    A = os.popen('conddb list SiStripApvGain_FromParticles_GR10_v1_express | grep SiStripApvGain').read().split()
    B = []
    previousRun = 0
    for i in range(5,len(A),5):
        B.append([int(A[i-5]), int(A[i])-1, A[i-2]])
    B.append([int(A[len(A)-5]), 999999, A[len(A)-2]])
    return B

def GetSQLiteCkSum (path):
    A = os.popen("sqlite3 -line %s/%s/sqlite/Gains_Sqlite.db 'select PAYLOAD_HASH from IOV;'" % (gains, path)).read().split()
    return "invalid" if (len(A) < 3) else A[2]

def AssociateChecksums (ListOfAppliedIntervals):
    UpdatedList = []
    for entry in ListOfAppliedIntervals:
        UpdatedList.append(entry)

        for directory in os.listdir(gains):
            if os.path.isfile("%s/%s/sqlite/Gains_Sqlite.db" % (gains, directory)):
                cksum = GetSQLiteCkSum(directory)
                if cksum == "invalid": continue
                if cksum == entry[2]:
                    UpdatedList[len(UpdatedList)-1].append(directory)
                    break

    return UpdatedList

def CleanAssociatedList (UpdatedList):
    toReturn = []
    for entry in UpdatedList:
        if len(entry) == 4: toReturn.append(entry)
    for i in range (1, len(toReturn)):
        if toReturn[i-1][1] != toReturn[i][0]-1: toReturn[i-1][1] = toReturn[i][0]-1
    toReturn[len(toReturn)-1][1] = 999999
    return toReturn

def GetAppliedGains ():
    toReturn = GetAppliedIntervals()
    toReturn = AssociateChecksums(toReturn)
    toReturn = CleanAssociatedList(toReturn)
    return toReturn

def GetAvailableGains ():
    toReturn = []
    for directory in os.listdir(gains):
        if not os.path.isdir ("%s/%s" % (gains, directory)): continue
        if directory.find("run") == -1 and directory.find("Run") == -1: continue
        if directory.find("Aag") >= 0: continue
        if directory.find("PCL") == -1 and directory.find("CalibTree") == -1: continue
        if len(os.listdir("%s/%s/plots_gain" % (gains, directory))) == 1: continue
        if not os.path.isfile("%s/%s/sqlite/Gains.root" % (gains, directory)): continue
        toManipulate = directory.split("_")
        if int(toManipulate[1]) < 207885: continue # only 2015 onwards
        toReturn.append([int(toManipulate[1]), int(toManipulate[3]), directory])
    for i in range(1, len(toReturn)):
        toReturn[i-1][1] = toReturn[i][0]-1
    toReturn[len(toReturn)-1][1] = 999999
    return toReturn

def ForceWhitelist (inputList):
    toReturn = []
    for entry in inputList:
       for j in range(0, len(whitelist)):
           if (entry[0] == whitelist[j]):
               toReturn.append(entry)
               break
    return toReturn

def ForceBlacklist (inputList):
    toReturn = []
    for entry in inputList:
       Match = False
       for j in range(0, len(blacklist)):
           if (entry[0] == blacklist[j]):
               Match = True
               break
       if not Match:
          toReturn.append(entry)
    return toReturn

def CombineGainsList (correctGains, appliedGains):
    toReturn = []
    for corGain in correctGains:
       for j in range(0, len(appliedGains)):
          if corGain[0] >= appliedGains[j][0] and corGain[1] <= appliedGains[j][1]:
             toReturn.append([corGain[0], corGain[1], corGain[2], appliedGains[j][3]])
          if corGain[0] >= appliedGains[j][0] and corGain[0] < appliedGains[j][1] and corGain[1] > appliedGains[j][1]:
             toReturn.append([corGain[0], appliedGains[j][1], corGain[2], appliedGains[j][3]])
             toReturn.append([appliedGains[j][1]+1,corGain[1], corGain[2], appliedGains[j+1][3]])
    for i in range(1, len(toReturn)):
        toReturn[i-1][1] = toReturn[i][0]-1
    toReturn[len(toReturn)-1][1] = 999999
    return toReturn

def WriteTXT (TXTName, correctGains):
    f = open (TXTName, 'w')
    for entry in correctGains:
        f.write ("Gains_%i_to_%i %s %s\n" % (entry[0], entry[1], entry[2], entry[3]))
    f.close()

if sys.argv[1] == '1':
    appliedGains = GetAppliedGains()
    correctGains = GetAvailableGains()
    if useWhitelist: correctGains = ForceWhitelist(correctGains)
    if useBlacklist: correctGains = ForceBlacklist(correctGains)

    promptGains = CombineGainsList (correctGains, appliedGains)
    WriteTXT ("gains.txt", promptGains)

if sys.argv[1] == '2':
    os.system('rm -rf Gains')
    os.system('sh CombineGains.sh')
    os.system('root -l -b -q ReorganizeGains.C+')
    os.system('hadd -f Stuff.root Gains/*.root')

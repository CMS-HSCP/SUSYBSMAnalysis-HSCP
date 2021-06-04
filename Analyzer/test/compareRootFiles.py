import ROOT
from ROOT import TH1F, TH2F, TFile, TMath, TH1D


#####################################################################################################################
#
#   This Python script takes 2 files in entries (here: HistosEDAnalyzer.root & HistosStep1.root) and 
#   return the differences between those files in differences.txt
#   The tests which are done are: check of the number of entries, of the mean, of the standard deviation & 
#   Kolmogorov test which is equal to 1 if both histograms are the same.
#   There is a flag at -1 in order to know if the Kolmogorov test is bad because of different number of bins.
#   This script does not provide any test on TTrees in the root files. 
#
#####################################################################################################################


# function which returns name and object from a root file
def getall(d, basepath="/"):
    "Generator function to recurse into a ROOT file/dir and yield (path, name, obj_name, obj) tuples"
    if d.ClassName() == 'TTree':
        return 
    for key in d.GetListOfKeys():
        kname = key.GetName()
        obj = key.ReadObj()
        if key.IsFolder():
            for i in getall(d.Get(kname), basepath+kname+"/"):
                yield i
        else:
                yield basepath+kname, kname, d.Get(kname), obj

ROOT.gROOT.SetBatch(True)

ifile1 = TFile("HistosEDAnalyzer.root") # /analyzer/Data/histograms... format to use -> hardcoded @ line 54
ifile2 = TFile("HistosStep1.root") # don't care about folders -> it is the file of reference  

ofile = open("differences.txt","w") # outputfile where are saved the differences observed  

counter=0
for i,j,k,l in getall(ifile2):
    counter+=1

ofile.write("Number of object other than TTree: "+str(counter)+"\n"+"\n") # count & save the number of objects in the root file

for i,j,k,l in getall(ifile2): # loop under all the entries of the root file
            y = ifile1.Get("analyzer/Data/"+str(j))
            test = -1
            passbins = False
            entries1 = y.GetEntries()
            entries2 = l.GetEntries()
            mean1 = y.GetMean()
            mean2 = l.GetMean()
            dev1 = y.GetStdDev()
            dev2 = l.GetStdDev()
            nbinsx1 = y.GetNbinsX()
            nbinsy1 = y.GetNbinsY()
            nbinsx2 = l.GetNbinsX()
            nbinsy2 = l.GetNbinsY()
            if nbinsx1==nbinsx2 and nbinsy1==nbinsy2:
                passbins=True
            if passbins==True:
                test = y.KolmogorovTest(l)
            if test!=1:
                ofile.write(j+'\t entries1: '+str(entries1)+'\t mean1: '+str(mean1)+'\t dev1: '+str(dev1)+'\t entries2: '+str(entries2)+'\t mean2: '+str(mean2)+'\t dev2: '+str(dev2)+'\t'+"kolomogorov_test: "+str(test)+'\n')

ofile.close()
                

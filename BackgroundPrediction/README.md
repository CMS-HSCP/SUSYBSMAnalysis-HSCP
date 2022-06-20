# BackgroundStudies

This repositery is dedicated to background estimate and validation of the method in data

## Run background estimate + validation: read histograms & produce mass plots 

Configure configFile_readHist.txt respecting dedicated columns
```bash
root_file  nofPseudoExperiments     cutIndex     rebin_boolean   rebin_eta   rebin_Ih    rebin_p     rebin_mass  
```
where:
- nofPseudoExperiments is the number of pseudo-experiments done during the background estimate.
- cutIndex is the CutIndex (cut on pT and Ias) on which the prediction is done.
- rebin is a boolean in order to know if we want to rebin the different distributions. Then we give the values of rebinning for each distribution. 

then run,
```bash
root -l -q -b step2_backgroundPrediction.C
```

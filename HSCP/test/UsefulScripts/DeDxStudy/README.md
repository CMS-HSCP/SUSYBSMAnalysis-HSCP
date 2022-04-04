# DeDx Study

This document is about the code used by the HSCP analisys to evaluete the **DeDx** contribution from data.
The original code is `/Users/claup/workspace/Displaced/SUSYBSMAnalysis-HSCP/test/UsefulScripts/DeDxStudy`

- *l2-l23* basic include for ROOT and standard libraries
- *l26-l30* various namespace definition plus different classes enumeration
- *l32-l54* CMSSW dependencies enclosed in `endif`
	- include DataFormats and FWCore
	- include __Analysis Step1__ **!!! <-----**
- *l57-l59* functions declaration: `DistToHSCP`, `isCompatibleWithCosmic`, `GetMass`
- *l62-l70* constant values used somewhere:
```c++
const double P_Min               = 1   ;
const double P_Max               = 16  ; // 1 + 14 + 1; final one is for pixel!
const int    P_NBins             = 15  ; // 15th bin = pixel; 0 is underflow
const double Path_Min            = 0.2 ;
const double Path_Max            = 1.6 ;
const int    Path_NBins          = 42  ;
const double Charge_Min          = 0   ;
const double Charge_Max          = 5000;
const int    Charge_NBins        = 500 ;
```
- *l72-l218* definition fo the `dEdxStudyObj` structure (NOTE: we could convert it in a class)
```c++
void DeDxStudy(string DIRNAME="COMPILE", string INPUT="dEdx.root", string OUTPUT="out.root")
```
	- `DIRNAME`, directory containing the script `DeDxStudy.sh`
		- If the script `DeDxStudy.sh` it is run without argument, the variable `DIRNAME` will assume the default value (it is used for compilation)
	- `INPUT`, could contains the path to a .csv of dEdx.root obtained from AlCaReco
	- `OUTPUT`, it's clear
	- *l224-234* setting plot style
	- *l236-l254* check if running on data or simulation, boolean for _removingCosmics_. Read the input file either picking paths from the .csv or getting them in a differet way **!!! <----** 
	- *l257* `dedxGainCorrector trackerCorrector` find definition 
	- *l260* `HIPemulator`, find definition 
	- *l261-l268* TH3F DeDxTemplate declaretion. Find `loadDeDxTemplate`
```c++
TH3F* dEdxTemplates      = NULL;
TH3F* dEdxTemplatesIn    = NULL;
TH3F* dEdxTemplatesInc   = NULL;
TH3F* dEdxTemplatesCCC   = NULL;
TH3F* dEdxTemplatesCCC16 = NULL;
TH3F* dEdxTemplatesCC    = NULL;
TH3F* dEdxTemplatesCI    = NULL;
```
	- *l276-310* loading DeDxTemplate	
	- *l317* declare vector of _dEdxStudyObj_ named **result**
	- *l318-l408* push back in **result** all the objects, it used among the inputs also _trackerCorrector_
	- *l409* Loic's progess bar
	- *l410-l692* Loop on the input files **!! <-----**
		- *l414-l690* Loop on events
			- *l420-437* check if Run is changed. If data are processed, load the template refered to the specific run and add it the the **results** vector
			- *l439-l451* load dedxCollection, tracks and genParticle (if is signal mc) from the event
```c++
fwlite::Handle<DeDxHitInfoAss> dedxCollH;
dedxCollH.getByLabel(ev, "dedxHitInfo");
fwlite::Handle< std::vector<reco::Track> > trackCollHandle;
trackCollHandle.getByLabel(ev,"RefitterForDeDx");
fwlite::Handle< std::vector<reco::GenParticle> > genCollHandle;
genCollHandle.getByLabel(ev, "genParticlesSkimmed");
```
			- *l471-l688* loop on tracks **<----**
				- *l473-l490* basic track quality cut and loading dedx informations
				- *l493-l570* getting dedx infos ONLY IF for MIPs `track->pt() > 5)`
					- *l494-l569* loop on dedx
						- *l495-l498* get the DetId of the detector module corresponding to the dedx measurement
						- *l500-l568* loop on **results** vector  if related to hit infos **<---**
```c++
if(!results[R]->isHit) continue; //only consider results related to hit infos
```
							- *l503-l505* initialized scale factors and normalization
```c++
double scaleFactor = dEdxSF[0];
if (detid.subdetId()<3) scaleFactor *= dEdxSF[1];
double Norm = (detid.subdetId()<3)?3.61e-06:3.61e-06*265;
```
							- *l508-l514* if condition to force using only results with hit info, plus other cross-checks like skip dedx in pixel if results is not pixel releted. Same cross-check for strips, or if insideTkModule, or clusterCleaning
							- *l517-l531* in case of crossTalk, give the correct cluster charge
						- *l575-l688* loop on **results** if related to estimator/discriminator variables **<---**		 	 
```c++
if(!results[R]->isEstim and !results[R]->isDiscrim) continue; //only consider results related to estimator/discriminator variables here
```


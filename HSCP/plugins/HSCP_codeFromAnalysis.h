// Original Author:  Loic Quertenmont
// Modified by Caroline for the intership
//

#include "TVector3.h"
#include <string>
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <TTree.h>
#include <string.h>
#include "TH3.h"
#include "TFile.h"
#include "TObject.h"

using namespace edm;
using namespace reco;
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// general purpose code

// compute deltaR between two point (eta,phi) (eta,phi)
double deltaR(double eta1, double phi1, double eta2, double phi2);
TObject* GetObjectFromPath(TDirectory* File, std::string Path, bool GetACopy);

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Genertic code related to samples processing in FWLITE --> functions below will be loaded only if FWLITE compiler variable is defined

bool isGoodGenHSCP(const reco::GenParticle& gen, bool onlyCharged);
int HowManyChargedHSCP(const std::vector<reco::GenParticle>& genColl);

TH3F* loadDeDxTemplate(std::string path, bool splitByModuleType = false);
reco::DeDxData computedEdx(const reco::DeDxHitInfo* dedxHits,
                           double* scaleFactors,
                           TH3* templateHisto = NULL,
                           bool usePixel = false,
                           bool useClusterCleaning = true,
                           bool reverseProb = false,
                           bool useTruncated = false,
                           std::unordered_map<unsigned int, double>* TrackerGains = NULL,
                           bool useStrip = true,
                           bool mustBeInside = false,
                           size_t MaxStripNOM = 999,
                           bool correctFEDSat = false,
                           int crossTalkInvAlgo = 0,
                           double dropLowerDeDxValue = 0.0,
                           double* dEdxErr = NULL);

bool clusterCleaning(std::vector<int> ampls, int crosstalkInv = 0, uint8_t* exitCode = NULL);
void printClusterCleaningMessage(uint8_t exitCode);
std::vector<int> convert(const std::vector<unsigned char>& input);
std::vector<int> CrossTalkInv(const std::vector<int>& Q,
                              const float x1 = 0.10,
                              const float x2 = 0.04,
                              bool way = true,
                              float threshold = 20,
                              float thresholdSat = 25);

bool isHitInsideTkModule(const LocalPoint hitPos, const DetId& detid, const SiStripCluster* cluster);

reco::DeDxData computedEdx(const reco::DeDxHitInfo* dedxHits,
                           double* scaleFactors,
                           TH3* templateHisto,
                           bool usePixel,
                           bool useClusterCleaning,
                           bool reverseProb,
                           bool useTruncated,
                           std::unordered_map<unsigned int, double>* TrackerGains,
                           bool useStrip,
                           bool mustBeInside,
                           size_t MaxStripNOM,
                           bool correctFEDSat,
                           int crossTalkInvAlgo,
                           double dropLowerDeDxValue,
                           double* dEdxErr);

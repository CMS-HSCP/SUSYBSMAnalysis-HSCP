// Original Author:  Loic Quertenmont
// Modified by Caroline for the intership
//

#include "HSCP_codeFromAnalysis.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// general purpose code
//

// compute deltaR between two point (eta,phi) (eta,phi)
double deltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = phi1 - phi2;
  while (dphi > M_PI)
    dphi -= 2 * M_PI;
  while (dphi <= -M_PI)
    dphi += 2 * M_PI;
  return sqrt(deta * deta + dphi * dphi);
}

TObject* GetObjectFromPath(TDirectory* File, std::string Path, bool GetACopy = false) {
  size_t pos = Path.find("/");
  if (pos < 256) {
    std::string firstPart = Path.substr(0, pos);
    std::string endPart = Path.substr(pos + 1, Path.length());
    TDirectory* TMP = (TDirectory*)File->Get(firstPart.c_str());
    if (TMP != NULL)
      return GetObjectFromPath(TMP, endPart, GetACopy);

    printf("ObjectNotFound: %s::%s\n", File->GetName(), Path.c_str());
    return NULL;
  } else {
    if (GetACopy) {
      return (File->Get(Path.c_str()))->Clone();
    } else {
      return File->Get(Path.c_str());
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Genertic code related to samples processing in FWLITE --> functions below will be loaded only if FWLITE compiler variable is defined

//selection of gen HSCP
bool isGoodGenHSCP(const reco::GenParticle& gen, bool onlyCharged) {
  if (gen.status() != 1)
    return false;
  if (gen.pt() < 5)
    return false;
  int AbsPdg = abs(gen.pdgId());
  if (AbsPdg < 1000000 && AbsPdg != 17)
    return false;
  if (onlyCharged && (AbsPdg == 1000993 || AbsPdg == 1009313 || AbsPdg == 1009113 || AbsPdg == 1009223 ||
                      AbsPdg == 1009333 || AbsPdg == 1092114 || AbsPdg == 1093214 || AbsPdg == 1093324))
    return false;  //Skip neutral gluino RHadrons
  if (onlyCharged && (AbsPdg == 1000622 || AbsPdg == 1000642 || AbsPdg == 1006113 || AbsPdg == 1006311 ||
                      AbsPdg == 1006313 || AbsPdg == 1006333))
    return false;  //skip neutral stop RHadrons
  if (onlyCharged && AbsPdg == 1000022)
    return false;  //skip neutralino
  return true;
}

// count the number of charged generated HSCP in the event --> this is needed to reweights the events for different gluino ball fraction starting from f=10% samples
int HowManyChargedHSCP(const std::vector<reco::GenParticle>& genColl) {
  int toReturn = 0;
  for (unsigned int g = 0; g < genColl.size(); g++) {
    if (isGoodGenHSCP(genColl[g], true))
      toReturn++;
  }
  return toReturn;
}

TH3F* loadDeDxTemplate(std::string path, bool splitByModuleType) {
  TFile* InputFile = new TFile(path.c_str());
  TH3F* DeDxMap_ = (TH3F*)GetObjectFromPath(InputFile, "Charge_Vs_Path");
  if (!DeDxMap_) {
    printf("dEdx templates in file %s can't be open\n", path.c_str());
    exit(0);
  }

  TH3F* Prob_ChargePath = (TH3F*)(DeDxMap_->Clone("Prob_ChargePath"));
  Prob_ChargePath->Reset();
  Prob_ChargePath->SetDirectory(0);

  if (!splitByModuleType) {
    Prob_ChargePath->RebinX(Prob_ChargePath->GetNbinsX() - 1);  // <-- do not include pixel in the inclusive
  }

  for (int i = 0; i <= Prob_ChargePath->GetXaxis()->GetNbins() + 1; i++) {
    for (int j = 0; j <= Prob_ChargePath->GetYaxis()->GetNbins() + 1; j++) {
      double Ni = 0;
      for (int k = 0; k <= Prob_ChargePath->GetZaxis()->GetNbins() + 1; k++) {
        Ni += DeDxMap_->GetBinContent(i, j, k);
      }

      for (int k = 0; k <= Prob_ChargePath->GetZaxis()->GetNbins() + 1; k++) {
        double tmp = 0;
        for (int l = 0; l <= k; l++) {
          tmp += DeDxMap_->GetBinContent(i, j, l);
        }

        if (Ni > 0) {
          Prob_ChargePath->SetBinContent(i, j, k, tmp / Ni);
        } else {
          Prob_ChargePath->SetBinContent(i, j, k, 0);
        }
      }
    }
  }
  InputFile->Close();
  return Prob_ChargePath;
}

const double TkModGeomThickness[] = {1,
                                     0.029000,
                                     0.029000,
                                     0.047000,
                                     0.047000,
                                     0.029000,
                                     0.029000,
                                     0.029000,
                                     0.029000,
                                     0.029000,
                                     0.029000,
                                     0.029000,
                                     0.047000,
                                     0.047000,
                                     0.047000};
const double TkModGeomLength[] = {1,
                                  5.844250,
                                  5.844250,
                                  9.306700,
                                  9.306700,
                                  5.542900,
                                  4.408000,
                                  5.533000,
                                  4.258000,
                                  4.408000,
                                  5.533000,
                                  5.758000,
                                  7.363125,
                                  9.204400,
                                  10.235775};
const double TkModGeomWidthB[] = {1,
                                  3.072000,
                                  3.072000,
                                  4.684800,
                                  4.684800,
                                  4.579445,
                                  5.502407,
                                  3.158509,
                                  4.286362,
                                  5.502407,
                                  4.049915,
                                  3.561019,
                                  6.002559,
                                  5.235483,
                                  3.574395};
const double TkModGeomWidthT[] = {1,
                                  3.072000,
                                  3.072000,
                                  4.684800,
                                  4.684800,
                                  3.095721,
                                  4.322593,
                                  4.049915,
                                  3.146580,
                                  4.322593,
                                  3.158509,
                                  2.898798,
                                  4.824683,
                                  4.177638,
                                  4.398049};

bool isHitInsideTkModule(const LocalPoint hitPos, const DetId& detid, const SiStripCluster* cluster = NULL) {
  if (detid.subdetId() < 3) {
    return true;
  }  //do nothing for pixel modules
  SiStripDetId SSdetId(detid);
  int moduleGeometry = SSdetId.moduleGeometry();

  //clean along the apv lines
  if (cluster &&
      (cluster->firstStrip() % 128 == 0 || (cluster->firstStrip() + cluster->amplitudes().size() % 128 == 127)))
    return false;

  double nx, ny;
  if (moduleGeometry <= 4) {
    ny = hitPos.y() / TkModGeomLength[moduleGeometry];
    nx = hitPos.x() / TkModGeomWidthT[moduleGeometry];
  } else {
    double offset =
        TkModGeomLength[moduleGeometry] * (TkModGeomWidthT[moduleGeometry] + TkModGeomWidthB[moduleGeometry]) /
        (TkModGeomWidthT[moduleGeometry] -
         TkModGeomWidthB[moduleGeometry]);  // check sign if GeomWidthT[moduleGeometry] < TkModGeomWidthB[moduleGeometry] !!!
    double tan_a = TkModGeomWidthT[moduleGeometry] / std::abs(offset + TkModGeomLength[moduleGeometry]);
    ny = hitPos.y() / TkModGeomLength[moduleGeometry];
    nx = hitPos.x() / (tan_a * std::abs(hitPos.y() + offset));
  }

  // "blacklists" for the gaps and edges
  switch (moduleGeometry) {
    case 0:
      return true;
    case 1:
      if (fabs(ny) > 0.96 || fabs(nx) > 0.98)
        return false;
      break;
    case 2:
      if (fabs(ny) > 0.97 || fabs(nx) > 0.99)
        return false;
      break;
    case 3:
      if (fabs(ny) > 0.98 || fabs(nx) > 0.98 || fabs(ny) < 0.04)
        return false;
      break;
    case 4:
      if (fabs(ny) > 0.98 || fabs(nx) > 0.98 || fabs(ny) < 0.04)
        return false;
      break;
    case 5:
      if (fabs(ny) > 0.98 || fabs(nx) > 0.98)
        return false;
      break;
    case 6:
      if (fabs(ny) > 0.98 || fabs(nx) > 0.99)
        return false;
      break;
    case 7:
      if (fabs(ny) > 0.97 || fabs(nx) > 0.98)
        return false;
      break;
    case 8:
      if (fabs(ny) > 0.98 || fabs(nx) > 0.99)
        return false;
      break;
    case 9:
      if (fabs(ny) > 0.98 || fabs(nx) > 0.99)
        return false;
      break;
    case 10:
      if (fabs(ny) > 0.97 || fabs(nx) > 0.99)
        return false;
      break;
    case 11:
      if (fabs(ny) > 0.97 || fabs(nx) > 0.99)
        return false;
      break;
    case 12:
      if (fabs(ny) > 0.98 || fabs(nx) > 0.99 || (-0.17 < ny && ny < -0.07))
        return false;
      break;
    case 13:
      if (fabs(ny) > 0.97 || fabs(nx) > 0.99 || (-0.10 < ny && ny < -0.01))
        return false;
      break;
    case 14:
      if (fabs(ny) > 0.95 || fabs(nx) > 0.98 || (0.01 < ny && ny < 0.12))
        return false;
      break;
    default:
      std::cerr << "Unknown module geometry! Exiting!" << std::endl;
      exit(EXIT_FAILURE);
  }

  return true;
}

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
                           double* dEdxErr) {
  if (!dedxHits)
    return reco::DeDxData(-1, -1, -1);
  //     if(templateHisto)usePixel=false; //never use pixel for discriminator

  std::vector<double> vect;
  std::vector<double> vectStrip;
  std::vector<double> vectPixel;

  bool debugprint = false;

  unsigned int NSat = 0;
  unsigned int SiStripNOM = 0;
  for (unsigned int h = 0; h < dedxHits->size(); h++) {
    if (debugprint)
      std::cout << "look on dedxHits in computedEdx " << h << std::endl;
    DetId detid(dedxHits->detId(h));
    if (!usePixel && detid.subdetId() < 3)
      continue;  // skip pixels
    if (!useStrip && detid.subdetId() >= 3)
      continue;  // skip strips
    //        if(useClusterCleaning && !clusterCleaning(dedxHits->stripCluster(h), crossTalkInvAlgo))continue;
    bool bool_cleaning = true;
    if (useStrip && detid.subdetId() >= 3) {
      std::vector<int> amps = convert(dedxHits->stripCluster(h)->amplitudes());
      if (crossTalkInvAlgo == 1) {
        amps = CrossTalkInv(amps, 0.10, 0.04, true);
      }
      if (debugprint)
        std::cout << " after amps " << std::endl;
      bool_cleaning = clusterCleaning(amps, crossTalkInvAlgo);
    }
    if (useClusterCleaning && !bool_cleaning)
      continue;
    //if(useClusterCleaning && !clusterCleaning(amps, crossTalkInvAlgo))continue;
    //printStripCluster(stdout, dedxHits->stripCluster(h), dedxHits->detId(h));

    if (debugprint)
      std::cout << " after clusterCleaning " << std::endl;
    if (mustBeInside &&
        !isHitInsideTkModule(dedxHits->pos(h), detid, detid.subdetId() >= 3 ? dedxHits->stripCluster(h) : NULL))
      continue;
    if (debugprint)
      std::cout << " after mustBeInside " << std::endl;
    if (detid.subdetId() >= 3 && ++SiStripNOM > MaxStripNOM)
      continue;  // skip remaining strips, but not pixel
    if (debugprint)
      std::cout << " after MaxStripNOM " << std::endl;

    int ClusterCharge = dedxHits->charge(h);

    if (detid.subdetId() >= 3) {  //for strip only
      const SiStripCluster* cluster = dedxHits->stripCluster(h);
      std::vector<int> amplitudes = convert(cluster->amplitudes());
      if (crossTalkInvAlgo)
        amplitudes = CrossTalkInv(amplitudes, 0.10, 0.04, true);
      if (debugprint)
        std::cout << " after amps2 " << std::endl;
      int firstStrip = cluster->firstStrip();
      int prevAPV = -1;
      double gain = 1.0;

      bool isSatCluster = false;
      ClusterCharge = 0;
      for (unsigned int s = 0; s < amplitudes.size(); s++) {
        if (TrackerGains != NULL) {  //don't reload the gain if unnecessary  since map access are slow
          int APV = (firstStrip + s) / 128;
          if (APV != prevAPV) {
            gain = TrackerGains->at(detid.rawId() << 3 | APV);
            prevAPV = APV;
          }
        }

        int StripCharge = amplitudes[s];
        if (StripCharge < 254) {
          StripCharge = (int)(StripCharge / gain);
          if (StripCharge >= 1024) {
            StripCharge = 255;
          } else if (StripCharge >= 254) {
            StripCharge = 254;
          }
        }

        if (StripCharge >= 254) {
          isSatCluster = true;
        }
        if (StripCharge >= 255 && correctFEDSat) {
          StripCharge = 512;
        }
        ClusterCharge += StripCharge;
      }
      if (isSatCluster)
        NSat++;
    }
    if (debugprint)
      std::cout << " after gain corrected " << std::endl;

    double scaleFactor = scaleFactors[0];
    if (detid.subdetId() < 3)
      scaleFactor *= scaleFactors[1];  // add pixel scaling
    if (debugprint)
      std::cout << " after SF " << std::endl;

    if (templateHisto) {  //save discriminator probability
      double ChargeOverPathlength =
          scaleFactor * ClusterCharge / (dedxHits->pathlength(h) * 10.0 * (detid.subdetId() < 3 ? 265 : 1));
      //           if(fakeHIP && detid.subdetId()>=3 && rand()%1000<35)ChargeOverPathlength = ( 0.5 + ((rand()%15000)/10000.0) ) / (3.61e-06*265*10);
      //           if(fakeHIP && detid.subdetId() <3 && rand()%1000<20)ChargeOverPathlength = ( 0.3 + ((rand()%12000)/10000.0) ) / (3.61e-06*265*10*265);

      int moduleGeometry = 0;  // underflow for debug
      if (detid.subdetId() < 3)
        moduleGeometry = 15;  // 15 == pixel
      else {
        SiStripDetId SSdetId(detid);
        moduleGeometry = SSdetId.moduleGeometry();
      }
      int BinX = templateHisto->GetXaxis()->FindBin(moduleGeometry);
      int BinY = templateHisto->GetYaxis()->FindBin(dedxHits->pathlength(h) * 10.0);  //*10 because of cm-->mm
      int BinZ = templateHisto->GetZaxis()->FindBin(ChargeOverPathlength);
      double Prob = templateHisto->GetBinContent(BinX, BinY, BinZ);
      //printf("%i %i %i  %f\n", BinX, BinY, BinZ, Prob);
      if (reverseProb)
        Prob = 1.0 - Prob;
      vect.push_back(Prob);  //save probability
      if (debugprint)
        std::cout << " after Prob vect.push_back " << std::endl;
    } else {
      double Norm = (detid.subdetId() < 3) ? 3.61e-06 : 3.61e-06 * 265;
      double ChargeOverPathlength = scaleFactor * Norm * ClusterCharge / dedxHits->pathlength(h);

      vect.push_back(ChargeOverPathlength);  //save charge
      if (detid.subdetId() < 3)
        vectPixel.push_back(ChargeOverPathlength);
      if (detid.subdetId() >= 3)
        vectStrip.push_back(ChargeOverPathlength);
      //           printf("%i - %f / %f = %f\n", h, scaleFactor*Norm*dedxHits->charge(h), dedxHits->pathlength(h), ChargeOverPathlength);
      if (debugprint)
        std::cout << " after ChargeOverPathlength vect.push_back " << std::endl;
    }
  }

  if (dropLowerDeDxValue > 0) {
    std::vector<double> tmp(vect.size());
    std::copy(vect.begin(), vect.end(), tmp.begin());
    std::sort(tmp.begin(), tmp.end(), std::greater<double>());
    int nTrunc = tmp.size() * dropLowerDeDxValue;
    vect.clear();
    for (unsigned int t = 0; t + nTrunc < tmp.size(); t++) {
      vect.push_back(tmp[t]);
    }
  }
  if (debugprint)
    std::cout << " after dropLowerDeDxValue " << std::endl;

  double result;
  int size = vect.size();

  if (size > 0) {
    if (templateHisto) {  //dEdx discriminator
      //Prod discriminator
      //result = 1;
      //for(int i=0;i<size;i++){
      //   if(vect[i]<=0.0001){result *= pow(0.0001 , 1.0/size);}
      //   else               {result *= pow(vect[i], 1.0/size);}
      //}

      //Ias discriminator
      result = 1.0 / (12 * size);
      std::sort(vect.begin(), vect.end(), std::less<double>());
      for (int i = 1; i <= size; i++) {
        result += vect[i - 1] * pow(vect[i - 1] - ((2.0 * i - 1.0) / (2.0 * size)), 2);
      }
      result *= (3.0 / size);
      if (debugprint)
        std::cout << " Ias discriminator " << result << std::endl;
    } else {  //dEdx estimator
      if (useTruncated) {
        //truncated40 estimator
        std::sort(vect.begin(), vect.end(), std::less<double>());
        result = 0;
        int nTrunc = size * 0.40;
        for (int i = 0; i + nTrunc < size; i++) {
          result += vect[i];
        }
        result /= (size - nTrunc);
        if (debugprint)
          std::cout << " truncated40 discriminator " << result << std::endl;
      } else {
        //harmonic2 estimator
        result = 0;
        double expo = -2;
        if (dEdxErr)
          *dEdxErr = 0;
        for (int i = 0; i < size; i++) {
          result += pow(vect[i], expo);
          if (dEdxErr)
            *dEdxErr += pow(vect[i], 2 * (expo - 1)) * pow(0.01, 2);
        }
        result = pow(result / size, 1. / expo);
        if (dEdxErr)
          *dEdxErr = result * result * result * sqrt(*dEdxErr) / size;
        if (debugprint)
          std::cout << " harmonic2 discriminator " << result << std::endl;
      }
      //           printf("Ih = %f\n------------------\n",result);
    }
  } else {
    result = -1;
  }
  if (debugprint)
    std::cout << " ok finished computeDeDx " << std::endl;
  return reco::DeDxData(result, NSat, size);
}

std::vector<int> convert(const std::vector<unsigned char>& input) {
  std::vector<int> output;
  for (unsigned int i = 0; i < input.size(); i++) {
    output.push_back((int)input[i]);
  }
  return output;
}

std::vector<int> CrossTalkInv(
    const std::vector<int>& Q, const float x1, const float x2, bool way, float threshold, float thresholdSat) {
  const unsigned N = Q.size();
  std::vector<int> QII;
  std::vector<float> QI(N, 0);
  Double_t a = 1 - 2 * x1 - 2 * x2;
  //  bool debugbool=false;
  TMatrix A(N, N);

  //---  que pour 1 max bien net
  if (Q.size() < 2 || Q.size() > 8) {
    for (unsigned int i = 0; i < Q.size(); i++) {
      QII.push_back((int)Q[i]);
    }
    return QII;
  }
  if (way) {
    std::vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end());
    if (*mQ > 253) {
      if (*mQ == 255 && *(mQ - 1) > 253 && *(mQ + 1) > 253)
        return Q;
      if (*(mQ - 1) > thresholdSat && *(mQ + 1) > thresholdSat && *(mQ - 1) < 254 && *(mQ + 1) < 254 &&
          abs(*(mQ - 1) - *(mQ + 1)) < 40) {
        QII.push_back((10 * (*(mQ - 1)) + 10 * (*(mQ + 1))) / 2);
        return QII;
      }
    }
  }
  //---

  for (unsigned int i = 0; i < N; i++) {
    A(i, i) = a;
    if (i < N - 1) {
      A(i + 1, i) = x1;
      A(i, i + 1) = x1;
    } else
      continue;
    if (i < N - 2) {
      A(i + 2, i) = x2;
      A(i, i + 2) = x2;
    }
  }

  if (N == 1)
    A(0, 0) = 1 / a;
  else
    A.InvertFast();

  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = 0; j < N; j++) {
      QI[i] += A(i, j) * (float)Q[j];
    }
  }

  for (unsigned int i = 0; i < QI.size(); i++) {
    if (QI[i] < threshold)
      QI[i] = 0;
    QII.push_back((int)QI[i]);
  }

  return QII;
}

//bool clusterCleaning(const SiStripCluster*   cluster,  int crosstalkInv, uint8_t * exitCode)
bool clusterCleaning(std::vector<int> ampls, int crosstalkInv, uint8_t* exitCode) {
  /*
   if(!cluster) return true;
   std::vector<int>  ampls = convert(cluster->amplitudes());
   if(crosstalkInv==1)ampls = CrossTalkInv(ampls,0.10,0.04, true);
*/

  // ----------------  COMPTAGE DU NOMBRE DE MAXIMA   --------------------------
  //----------------------------------------------------------------------------
  Int_t NofMax = 0;
  Int_t recur255 = 1;
  Int_t recur254 = 1;
  bool MaxOnStart = false;
  bool MaxInMiddle = false, MaxOnEnd = false;
  Int_t MaxPos = 0;
  // Début avec max
  if (ampls.size() != 1 &&
      ((ampls[0] > ampls[1]) ||
       (ampls.size() > 2 && ampls[0] == ampls[1] && ampls[1] > ampls[2] && ampls[0] != 254 && ampls[0] != 255) ||
       (ampls.size() == 2 && ampls[0] == ampls[1] && ampls[0] != 254 && ampls[0] != 255))) {
    NofMax = NofMax + 1;
    MaxOnStart = true;
  }

  // Maximum entouré
  if (ampls.size() > 2) {
    for (unsigned int i = 1; i < ampls.size() - 1; i++) {
      if ((ampls[i] > ampls[i - 1] && ampls[i] > ampls[i + 1]) ||
          (ampls.size() > 3 && i > 0 && i < ampls.size() - 2 && ampls[i] == ampls[i + 1] && ampls[i] > ampls[i - 1] &&
           ampls[i] > ampls[i + 2] && ampls[i] != 254 && ampls[i] != 255)) {
        NofMax = NofMax + 1;
        MaxInMiddle = true;
        MaxPos = i;
      }
      if (ampls[i] == 255 && ampls[i] == ampls[i - 1]) {
        recur255 = recur255 + 1;
        MaxPos = i - (recur255 / 2);
        if (ampls[i] > ampls[i + 1]) {
          NofMax = NofMax + 1;
          MaxInMiddle = true;
        }
      }
      if (ampls[i] == 254 && ampls[i] == ampls[i - 1]) {
        recur254 = recur254 + 1;
        MaxPos = i - (recur254 / 2);
        if (ampls[i] > ampls[i + 1]) {
          NofMax = NofMax + 1;
          MaxInMiddle = true;
        }
      }
    }
  }
  // Fin avec un max
  if (ampls.size() > 1) {
    if (ampls[ampls.size() - 1] > ampls[ampls.size() - 2] ||
        (ampls.size() > 2 && ampls[ampls.size() - 1] == ampls[ampls.size() - 2] &&
         ampls[ampls.size() - 2] > ampls[ampls.size() - 3]) ||
        ampls[ampls.size() - 1] == 255) {
      NofMax = NofMax + 1;
      MaxOnEnd = true;
    }
  }
  // Si une seule strip touchée
  if (ampls.size() == 1) {
    NofMax = 1;
  }

  // ---  SELECTION EN FONCTION DE LA FORME POUR LES MAXIMA UNIQUES ---------
  //------------------------------------------------------------------------
  //
  //               ____
  //              |    |____
  //          ____|    |    |
  //         |    |    |    |____
  //     ____|    |    |    |    |
  //    |    |    |    |    |    |____
  //  __|____|____|____|____|____|____|__
  //    C_Mnn C_Mn C_M  C_D  C_Dn C_Dnn
  //
  //   bool shapetest=true;
  bool shapecdtn = false;
  if (exitCode)
    *exitCode = 255;

  if (crosstalkInv == 1) {
    if (NofMax == 1) {
      shapecdtn = true;
      if (exitCode)
        *exitCode = 0;
    }
    return shapecdtn;
  }

  //      Float_t C_M;    Float_t C_D;    Float_t C_Mn;   Float_t C_Dn;   Float_t C_Mnn;  Float_t C_Dnn;
  Float_t C_M = 0.0;
  Float_t C_D = 0.0;
  Float_t C_Mn = 10000;
  Float_t C_Dn = 10000;
  Float_t C_Mnn = 10000;
  Float_t C_Dnn = 10000;
  Int_t CDPos;
  Float_t coeff1 = 1.7;
  Float_t coeff2 = 2.0;
  Float_t coeffn = 0.10;
  Float_t coeffnn = 0.02;
  Float_t noise = 4.0;

  if (NofMax == 1) {
    if (MaxOnStart == true) {
      C_M = (Float_t)ampls[0];
      C_D = (Float_t)ampls[1];
      if (ampls.size() < 3)
        shapecdtn = true;
      else if (ampls.size() == 3) {
        C_Dn = (Float_t)ampls[2];
        if (C_Dn <= coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise || C_D == 255)
          shapecdtn = true;
        else if (exitCode)
          *exitCode = 2;
      } else if (ampls.size() > 3) {
        C_Dn = (Float_t)ampls[2];
        C_Dnn = (Float_t)ampls[3];
        if ((C_Dn <= coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise || C_D == 255) &&
            C_Dnn <= coeff1 * coeffn * C_Dn + coeff2 * coeffnn * C_D + 2 * noise) {
          shapecdtn = true;
        } else if (exitCode)
          *exitCode = 3;
      }
    }

    if (MaxOnEnd == true) {
      C_M = (Float_t)ampls[ampls.size() - 1];
      C_D = (Float_t)ampls[ampls.size() - 2];
      if (ampls.size() < 3)
        shapecdtn = true;
      else if (ampls.size() == 3) {
        C_Dn = (Float_t)ampls[0];
        if (C_Dn <= coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise || C_D == 255)
          shapecdtn = true;
        else if (exitCode)
          *exitCode = 4;
      } else if (ampls.size() > 3) {
        C_Dn = (Float_t)ampls[ampls.size() - 3];
        C_Dnn = (Float_t)ampls[ampls.size() - 4];
        if ((C_Dn <= coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise || C_D == 255) &&
            C_Dnn <= coeff1 * coeffn * C_Dn + coeff2 * coeffnn * C_D + 2 * noise) {
          shapecdtn = true;
        } else if (exitCode)
          *exitCode = 5;
      }
    }

    if (MaxInMiddle == true) {
      C_M = (Float_t)ampls[MaxPos];
      int LeftOfMaxPos = MaxPos - 1;
      if (LeftOfMaxPos <= 0)
        LeftOfMaxPos = 0;
      int RightOfMaxPos = MaxPos + 1;
      if (RightOfMaxPos >= (int)ampls.size())
        RightOfMaxPos = ampls.size() - 1;
      //int after = RightOfMaxPos; int before = LeftOfMaxPos; if (after>=(int)ampls.size() ||  before<0)  std::cout<<"invalid read MaxPos:"<<MaxPos <<"size:"<<ampls.size() <<std::endl;
      if (ampls[LeftOfMaxPos] < ampls[RightOfMaxPos]) {
        C_D = (Float_t)ampls[RightOfMaxPos];
        C_Mn = (Float_t)ampls[LeftOfMaxPos];
        CDPos = RightOfMaxPos;
      } else {
        C_D = (Float_t)ampls[LeftOfMaxPos];
        C_Mn = (Float_t)ampls[RightOfMaxPos];
        CDPos = LeftOfMaxPos;
      }
      if (C_Mn < coeff1 * coeffn * C_M + coeff2 * coeffnn * C_D + 2 * noise || C_M == 255) {
        if (ampls.size() == 3)
          shapecdtn = true;
        else if (ampls.size() > 3) {
          if (CDPos > MaxPos) {
            if (ampls.size() - CDPos - 1 == 0) {
              C_Dn = 0;
              C_Dnn = 0;
            }
            if (ampls.size() - CDPos - 1 == 1) {
              C_Dn = (Float_t)ampls[CDPos + 1];
              C_Dnn = 0;
            }
            if (ampls.size() - CDPos - 1 > 1) {
              C_Dn = (Float_t)ampls[CDPos + 1];
              C_Dnn = (Float_t)ampls[CDPos + 2];
            }
            if (MaxPos >= 2) {
              C_Mnn = (Float_t)ampls[MaxPos - 2];
            } else if (MaxPos < 2)
              C_Mnn = 0;
          }
          if (CDPos < MaxPos) {
            if (CDPos == 0) {
              C_Dn = 0;
              C_Dnn = 0;
            }
            if (CDPos == 1) {
              C_Dn = (Float_t)ampls[0];
              C_Dnn = 0;
            }
            if (CDPos > 1) {
              C_Dn = (Float_t)ampls[CDPos - 1];
              C_Dnn = (Float_t)ampls[CDPos - 2];
            }
            if (ampls.size() - LeftOfMaxPos > 1 && MaxPos + 2 < (int)(ampls.size()) - 1) {
              C_Mnn = (Float_t)ampls[MaxPos + 2];
            } else
              C_Mnn = 0;
          }
          if ((C_Dn <= coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise || C_D == 255) &&
              C_Mnn <= coeff1 * coeffn * C_Mn + coeff2 * coeffnn * C_M + 2 * noise &&
              C_Dnn <= coeff1 * coeffn * C_Dn + coeff2 * coeffnn * C_D + 2 * noise) {
            shapecdtn = true;
          }
        }
      } else if (exitCode)
        *exitCode = 6;
    }
  } else if (NofMax > 1 && exitCode)
    *exitCode = 1;  // more than one maximum
  if (ampls.size() == 1) {
    shapecdtn = true;
  }
  if (shapecdtn && exitCode)
    *exitCode = 0;

  return shapecdtn;
}

void printClusterCleaningMessage(uint8_t exitCode) {
  switch (exitCode) {
    case 0:
      std::cout << "All went well" << std::endl;
      break;
    case 1:
      std::cout << "More than one maximum" << std::endl;
      break;
    case 2:
      std::cout << "MFirst; CSize=3; CDn too big" << std::endl;
      break;
    case 3:
      std::cout << "MFirst; CSize>3; CDn||CDnn too big" << std::endl;
      break;
    case 4:
      std::cout << "MEnd; CSize=3; CDn too big" << std::endl;
      break;
    case 5:
      std::cout << "MEnd; CSize>3; CDn||CDnn too big" << std::endl;
      break;
    case 6:
      std::cout << "MMid; Sides are too big" << std::endl;
      break;
    case 255:
      std::cout << "Failed all shape tests" << std::endl;
      break;
    default:
      std::cout << "Unknown exit code!" << std::endl;
  }
}

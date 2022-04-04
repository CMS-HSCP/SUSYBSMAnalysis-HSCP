// Original Author:  Loic Quertenmont

namespace reco {
  class Vertex;
  class Track;
  class GenParticle;
  class DeDxData;
  class MuonTimeExtra;
}  // namespace reco
namespace susybsm {
  class HSCParticle;
}
namespace fwlite {
  class ChainEvent;
}
namespace trigger {
  class TriggerEvent;
}
namespace edm {
  class TriggerResults;
  class TriggerResultsByName;
  class InputTag;
}  // namespace edm

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"

#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"

using namespace fwlite;
using namespace reco;
using namespace susybsm;
using namespace std;
using namespace edm;
using namespace trigger;

#include "../../AnalysisCode/Analysis_Step1_EventLoop.C"

#endif

double CutMass;
double PtCut;
double ICut;
double TOFCut;

bool isData = true;
bool isMC = !isData;

bool NotARegion(double Pt, double PtCut, double I, double ICut, double TOF, double TOFCut) {
  return (Pt > PtCut || (ICut > -1 && I > ICut) || (TOFCut > -1 && TOF <= TOFCut));
}

bool NotBRegion(double Pt, double PtCut, double I, double ICut, double TOF, double TOFCut) {
  return (Pt > PtCut || (ICut > -1 && I <= ICut) || (TOFCut > -1 && TOF <= TOFCut));
}

bool NotCRegion(double Pt, double PtCut, double I, double ICut, double TOF, double TOFCut) {
  return (Pt <= PtCut || (ICut > -1 && I > ICut) || (TOFCut > -1 && TOF <= TOFCut));
}

bool NotDRegion(double Pt, double PtCut, double I, double ICut, double TOF, double TOFCut) {
  return (Pt <= PtCut || (ICut > -1 && I <= ICut) || (TOFCut > -1 && TOF <= TOFCut));
}

bool NotERegion(double Pt, double PtCut, double I, double ICut, double TOF, double TOFCut) {
  return (Pt > PtCut || (ICut > -1 && I > ICut) || (TOFCut > -1 && TOF > TOFCut));
}

bool NotFRegion(double Pt, double PtCut, double I, double ICut, double TOF, double TOFCut) {
  return (Pt > PtCut || (ICut > -1 && I <= ICut) || (TOFCut > -1 && TOF > TOFCut));
}

bool NotGRegion(double Pt, double PtCut, double I, double ICut, double TOF, double TOFCut) {
  return (Pt <= PtCut || (ICut > -1 && I > ICut) || (TOFCut > -1 && TOF > TOFCut));
}

bool NotHRegion(double Pt, double PtCut, double I, double ICut, double TOF, double TOFCut) {
  return (Pt <= PtCut || (ICut > -1 && I <= ICut) || (TOFCut > -1 && TOF > TOFCut));
}

void DumpCandidateInfo(const susybsm::HSCParticle& hscp,
                       const fwlite::ChainEvent& ev,
                       FILE* pFile,
                       double treeMass,
                       muonTimingCalculator tofCalculator) {
  reco::MuonRef muon = hscp.muonRef();
  reco::TrackRef track = hscp.trackRef();
  if (TypeMode != 3 && track.isNull())
    return;

  reco::TrackRef SAtrack;
  if (!muon.isNull())
    SAtrack = muon->standAloneMuon();
  if (TypeMode == 3 && SAtrack.isNull())
    return;

  fwlite::Handle<std::vector<reco::Vertex> > vertexCollHandle;
  vertexCollHandle.getByLabel(ev, "offlinePrimaryVertices");
  if (!vertexCollHandle.isValid()) {
    printf("Vertex Collection NotFound\n");
    return;
  }
  std::vector<reco::Vertex> vertexColl = *vertexCollHandle;
  if (vertexColl.size() < 1) {
    printf("NO VERTEX\n");
    return;
  }
  const reco::Vertex& vertex = vertexColl[0];

  fwlite::Handle<DeDxHitInfoAss> dedxCollH;
  dedxCollH.getByLabel(ev, "dedxHitInfo");
  if (!dedxCollH.isValid()) {
    printf("Invalid dedxCollH\n");
    return;
  }

  fwlite::Handle<MuonTimeExtraMap> TOFCollH;
  TOFCollH.getByLabel(ev, "muons", TOF_Label.c_str());
  if (!TOFCollH.isValid()) {
    printf("Invalid TOF collection\n");
    return;
  }

  fwlite::Handle<MuonTimeExtraMap> TOFDTCollH;
  TOFDTCollH.getByLabel(ev, "muons", TOFdt_Label.c_str());
  if (!TOFDTCollH.isValid()) {
    printf("Invalid DT TOF collection\n");
    return;
  }

  fwlite::Handle<MuonTimeExtraMap> TOFCSCCollH;
  TOFCSCCollH.getByLabel(ev, "muons", TOFcsc_Label.c_str());
  if (!TOFCSCCollH.isValid()) {
    printf("Invalid CSC TOF collection\n");
    return;
  }

  fwlite::Handle<CSCSegmentCollection> CSCSegmentCollHandle;
  fwlite::Handle<DTRecSegment4DCollection> DTSegmentCollHandle;
  if (!isMC) {  //do not reocmpute TOF on MC background
    CSCSegmentCollHandle.getByLabel(ev, "cscSegments");
    if (!CSCSegmentCollHandle.isValid()) {
      printf("CSC Segment Collection not found!\n");
      return;
    }

    DTSegmentCollHandle.getByLabel(ev, "dt4DSegments");
    if (!DTSegmentCollHandle.isValid()) {
      printf("DT Segment Collection not found!\n");
      return;
    }
  }

  //load quantity associated to this track (TOF and dEdx)
  const DeDxHitInfo* dedxHits = NULL;
  if (TypeMode != 3 && !track.isNull()) {
    DeDxHitInfoRef dedxHitsRef = dedxCollH->get(track.key());
    if (!dedxHitsRef.isNull())
      dedxHits = &(*dedxHitsRef);
  }
  const reco::MuonTimeExtra* tof = NULL;
  const reco::MuonTimeExtra* dttof = NULL;
  const reco::MuonTimeExtra* csctof = NULL;
  if (TypeMode > 1 && TypeMode != 5 && !hscp.muonRef().isNull()) {
    if (isMC) {
      tof = &TOFCollH->get(hscp.muonRef().key());
      dttof = &TOFDTCollH->get(hscp.muonRef().key());
      csctof = &TOFCSCCollH->get(hscp.muonRef().key());
    } else {
      const CSCSegmentCollection& CSCSegmentColl = *CSCSegmentCollHandle;
      const DTRecSegment4DCollection& DTSegmentColl = *DTSegmentCollHandle;
      tofCalculator.computeTOF(
          muon, CSCSegmentColl, DTSegmentColl, isData ? 1 : 0);  //apply T0 correction on data but not on signal MC
      tof = &tofCalculator.combinedTOF;
      dttof = &tofCalculator.dtTOF;
      csctof = &tofCalculator.cscTOF;
    }
  }

  //Compute dE/dx on the fly
  //computedEdx(dedxHits, Data/MC scaleFactor, templateHistoForDiscriminator, usePixel, useClusterCleaning, reverseProb)
  DeDxData dedxSObjTmp = computedEdx(dedxHits,
                                     dEdxSF,
                                     dEdxTemplates,
                                     true,
                                     useClusterCleaning,
                                     TypeMode == 5,
                                     false,
                                     trackerCorrector.TrackerGains,
                                     true,
                                     true,
                                     99,
                                     false,
                                     1,
                                     0.0,
                                     NULL);
  DeDxData dedxMObjTmp = computedEdx(dedxHits,
                                     dEdxSF,
                                     NULL,
                                     true,
                                     useClusterCleaning,
                                     false,
                                     false,
                                     trackerCorrector.TrackerGains,
                                     true,
                                     true,
                                     99,
                                     false,
                                     1,
                                     0.0,
                                     NULL);

  DeDxData* dedxSObj = dedxSObjTmp.numberOfMeasurements() > 0 ? &dedxSObjTmp : NULL;
  DeDxData* dedxMObj = dedxMObjTmp.numberOfMeasurements() > 0 ? &dedxMObjTmp : NULL;
  //if(TypeMode==5)OpenAngle = deltaROpositeTrack(hscpColl, hscp); //OpenAngle is a global variable... that's uggly C++, but that's the best I found so far

  if (TypeMode != 3 && (track->pt() <= PtCut))
    return;  // || dedxSObj->dEdx()<=ICut))return;
  if (TypeMode == 3 && SAtrack->pt() < PtCut)
    return;
  //   if(track->pt()<=PtCut || dedxSObj->dEdx()<=ICut)return;
  if (TOFCut > -1 && tof && tof->inverseBeta() <= TOFCut)
    return;

  double Mass = 0;
  if (!track.isNull() && dedxMObj)
    Mass = GetMass(track->p(), dedxMObj->dEdx(), false);
  //   if(CutMass>=0 && Mass<CutMass)return;

  double v3d = 0;
  double dxy = 0;
  double dz = 0;
  int goodVerts = 0;
  if (!track.isNull()) {
    dz = track->dz(vertex.position());
    dxy = track->dxy(vertex.position());
    for (unsigned int i = 1; i < vertexColl.size(); i++) {
      if (fabs(vertexColl[i].z()) < 15 &&
          sqrt(vertexColl[i].x() * vertexColl[i].x() + vertexColl[i].y() * vertexColl[i].y()) < 2 &&
          vertexColl[i].ndof() > 3) {
        goodVerts++;
      } else {
        continue;
      }
      if (fabs(track->dz(vertexColl[i].position())) < fabs(dz)) {
        dz = track->dz(vertexColl[i].position());
        dxy = track->dxy(vertexColl[i].position());
      }
    }
    v3d = sqrt(dz * dz + dxy * dxy);
  }

  fprintf(pFile, "\n");
  fprintf(pFile, "%s\n", ev.getTFile()->GetName());
  fprintf(pFile,
          "---------------------------------------------------------------------------------------------------\n");
  fprintf(pFile, "Candidate Type = %i --> Mass : %7.2f (recompute-->%7.2f)\n", hscp.type(), treeMass, Mass);
  fprintf(pFile,
          "------------------------------------------ EVENT INFO ---------------------------------------------\n");
  fprintf(pFile,
          "Run=%i Lumi=%i Event=%llu BX=%i  Orbit=%i Store=%i\n",
          ev.eventAuxiliary().run(),
          ev.eventAuxiliary().luminosityBlock(),
          ev.eventAuxiliary().event(),
          ev.eventAuxiliary().luminosityBlock(),
          ev.eventAuxiliary().orbitNumber(),
          ev.eventAuxiliary().storeNumber());
  edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
  //   for(unsigned int i=0;i<tr.size();i++){
  //     fprintf(pFile, "Path %3i %50s --> %1i\n",i, tr.triggerName(i).c_str(),tr.accept(i));
  //   }
  fprintf(pFile,
          "Path %50s --> %1i\n",
          "HLT_PFMET170_NoiseCleaned_v*",
          (int)passTriggerPatterns(tr, "HLT_PFMET170_NoiseCleaned_v*"));
  fprintf(pFile, "Path %50s --> %1i\n", "HLT_Mu45_eta2p1_v*", (int)passTriggerPatterns(tr, "HLT_Mu45_eta2p1_v*"));
  fprintf(pFile, "Path %50s --> %1i\n", "HLT_Mu50_v*", (int)passTriggerPatterns(tr, "HLT_Mu50_v*"));

  if (!track.isNull()) {
    double validFractionTillLast =
        track->found() <= 0
            ? -1
            : track->found() /
                  float(track->found() +
                        track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS) +
                        track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS));
    fprintf(pFile,
            "------------------------------------------ INNER TRACKER ------------------------------------------\n");
    fprintf(pFile,
            "Quality = %i Chi2/NDF=%6.2f dz=+%6.2f dxy=%+6.2f V3D=%+6.2f (nGoodVert=%i) charge:%+i\n",
            track->qualityMask(),
            track->chi2() / track->ndof(),
            dz,
            dxy,
            v3d,
            goodVerts,
            track->charge());
    fprintf(pFile,
            "P=%7.2f  Pt=%7.2f+-%6.2f [%6.2f%%] (Cut=%6.2f) Eta=%+6.2f  Phi=%+6.2f  NOH=%2i FOVH=%6.2f NOMMH=%2i "
            "FOVHTL=%6.2f\n",
            track->p(),
            track->pt(),
            track->ptError(),
            track->ptError() / track->pt(),
            PtCut,
            track->eta(),
            track->phi(),
            track->found(),
            track->validFraction(),
            track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS) +
                track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS),
            validFractionTillLast);

    fprintf(pFile,
            "------------------------------------------ DEDX INFO ----------------------------------------------\n");
    fprintf(pFile,
            "dEdx for selection     :%6.2f (Cut=%6.2f) NOM %2i NOS %2i\n",
            dedxSObj->dEdx(),
            ICut,
            dedxSObj->numberOfMeasurements(),
            dedxSObj->numberOfSaturatedMeasurements());
    fprintf(pFile,
            "dEdx for mass reco     :%6.2f             NOM %2i NOS %2i  --> Beta dEdx = %6.2f\n",
            dedxMObj->dEdx(),
            dedxMObj->numberOfMeasurements(),
            dedxMObj->numberOfSaturatedMeasurements(),
            GetIBeta(dedxMObj->dEdx(), false));
    //fprintf(pFile,"dEdx for mass reco (NP):%6.2f             NOM %2i NOS %2i  --> Beta dEdx = %6.2f\n",dedxMNPObj->dEdx(),dedxMNPObj->numberOfMeasurements(),dedxMNPObj->numberOfSaturatedMeasurements(), GetIBeta(dedxMNPObj->dEdx(), false) );

    fprintf(pFile,
            "dEdx mass error     :%6.2f (1Sigma dEdx) or %6.2f (1Sigma P)\n",
            GetMass(track->p(), 0.95 * dedxMObj->dEdx(), false),
            GetMass(track->p() * (1 - track->ptError() / track->pt()), dedxMObj->dEdx(), false));

    for (unsigned int h = 0; h < dedxHits->size(); h++) {
      DetId detid(dedxHits->detId(h));
      if (detid.subdetId() < 3) {
        fprintf(pFile, "DetId = %7i --> Pixel Hit\n", detid.rawId());
      } else {
        printStripCluster(pFile, dedxHits->stripCluster(h), dedxHits->detId(h), true);
      }
    }
  }

  if (!muon.isNull() && tof != NULL) {
    fprintf(pFile,
            "------------------------------------------ MUON INFO ----------------------------------------------\n");
    double TOFMass;
    if (TypeMode != 3)
      TOFMass = GetTOFMass(track->p(), tof->inverseBeta());
    else
      TOFMass = GetTOFMass(SAtrack->p(), tof->inverseBeta());
    fprintf(pFile, "MassTOF = %7.2fGeV\n", TOFMass);
    fprintf(pFile,
            "Quality=%i type=%i P=%7.2f  Pt=%7.2f Eta=%+6.2f Phi=%+6.2f #Chambers=%i\n",
            muon->isQualityValid(),
            muon->type(),
            muon->p(),
            muon->pt(),
            muon->eta(),
            muon->phi(),
            muon->numberOfChambersNoRPC());
    if (!SAtrack.isNull())
      fprintf(pFile,
              "SA track P=%7.2f  Pt=%7.2f Eta=%+6.2f Phi=%+6.2f #Chambers=%i\n",
              SAtrack->p(),
              SAtrack->pt(),
              SAtrack->eta(),
              SAtrack->phi(),
              muon->numberOfMatchedStations());

    /*      for (int i=0; i<SAtrack->hitPattern().numberOfHits(); i++) {
	uint32_t pattern = SAtrack->hitPattern().getHitPattern(i);
	if (pattern == 0) break;
	if (SAtrack->hitPattern().muonHitFilter(pattern) &&
	    (int(SAtrack->hitPattern().getSubStructure(pattern)) == 1 ||
	     int(SAtrack->hitPattern().getSubStructure(pattern)) == 2) &&
	    SAtrack->hitPattern().getHitType(pattern) == 0) {
	}
      }
*/

    fprintf(pFile,
            "muonTimeDT      : NDOF=%2i InvBeta=%6.2f+-%6.2f (Cut=%6.2f) --> beta=%6.2f FreeInvBeta=%6.2f+-%6.2f "
            "VertexTime=%6.2f+-%6.2f\n",
            dttof->nDof(),
            dttof->inverseBeta(),
            dttof->inverseBetaErr(),
            TOFCut,
            (1.0 / dttof->inverseBeta()),
            dttof->freeInverseBeta(),
            dttof->freeInverseBetaErr(),
            dttof->timeAtIpInOut(),
            dttof->timeAtIpInOutErr());
    fprintf(pFile,
            "muonTimeCSC     : NDOF=%2i InvBeta=%6.2f+-%6.2f (Cut=%6.2f) --> beta=%6.2f FreeInvBeta=%6.2f+-%6.2f "
            "VertexTime=%6.2f+-%6.2f\n",
            csctof->nDof(),
            csctof->inverseBeta(),
            csctof->inverseBetaErr(),
            TOFCut,
            (1.0 / csctof->inverseBeta()),
            csctof->freeInverseBeta(),
            csctof->freeInverseBetaErr(),
            csctof->timeAtIpInOut(),
            csctof->timeAtIpInOutErr());
    fprintf(pFile,
            "muonTimeCombined: NDOF=%2i InvBeta=%6.2f+-%6.2f (Cut=%6.2f) --> beta=%6.2f FreeInvBeta=%6.2f+-%6.2f "
            "VertexTime=%6.2f+-%6.2f\n",
            tof->nDof(),
            tof->inverseBeta(),
            tof->inverseBetaErr(),
            TOFCut,
            (1.0 / tof->inverseBeta()),
            tof->freeInverseBeta(),
            tof->freeInverseBetaErr(),
            tof->timeAtIpInOut(),
            tof->timeAtIpInOutErr());
  }
  if (hscp.hasRpcInfo()) {
    fprintf(pFile,
            "------------------------------------------ RPC INFO -----------------------------------------------\n");
    fprintf(pFile, "isCandidate %i Beta=%6.2f\n", hscp.rpc().isCandidate, hscp.rpc().beta);
  }

  if (hscp.hasCaloInfo() && hscp.caloInfoRef()->ecalTime != -9999) {
    fprintf(pFile,
            "------------------------------------------ CALO INFO ----------------------------------------------\n");
    fprintf(pFile,
            "HCAL: E=%6.2f E3x3=%6.2f E5x5=%6.2f HO E=%6.2f\n",
            hscp.caloInfoRef()->hcalCrossedEnergy,
            hscp.caloInfoRef()->hcal3by3dir,
            hscp.caloInfoRef()->hcal5by5dir,
            hscp.caloInfoRef()->hoCrossedEnergy);
    fprintf(pFile,
            "ECAL: E=%6.2f E3x3=%6.2f E5x5=%6.2f\n",
            hscp.caloInfoRef()->ecalCrossedEnergy,
            hscp.caloInfoRef()->ecal3by3dir,
            hscp.caloInfoRef()->ecal5by5dir);
    fprintf(pFile,
            "ECAL: time=%6.2f beta=%6.2f trkisodr=%6.2f\n",
            hscp.caloInfoRef()->ecalTime,
            hscp.caloInfoRef()->ecalBeta,
            hscp.caloInfoRef()->trkIsoDr);
  }

  fprintf(pFile,
          "------------------------------------------ ISOL INFO ----------------------------------------------\n");
  fwlite::Handle<HSCPIsolationValueMap> IsolationH05;
  IsolationH05.getByLabel(ev, "HSCPIsolation", "R05");  //New format used for data since 17-07-2015
  if (!IsolationH05.isValid()) {
    IsolationH05.getByLabel(ev, "HSCPIsolation05");  //Old format used for first 2015B data, Signal and MC Backgrounds
    if (!IsolationH05.isValid()) {
      printf("Invalid IsolationH\n");
      return;
    }
  }
  const ValueMap<HSCPIsolation>& IsolationMap05 = *IsolationH05.product();

  fwlite::Handle<HSCPIsolationValueMap> IsolationH03;
  IsolationH03.getByLabel(ev, "HSCPIsolation", "R03");  //New format used for data since 17-07-2015
  if (!IsolationH03.isValid()) {
    IsolationH03.getByLabel(ev, "HSCPIsolation03");  //Old format used for first 2015B data, Signal and MC Backgrounds
    if (!IsolationH03.isValid()) {
      printf("Invalid IsolationH\n");
      return;
    }
  }
  const ValueMap<HSCPIsolation>& IsolationMap03 = *IsolationH03.product();

  fwlite::Handle<HSCPIsolationValueMap> IsolationH01;
  IsolationH01.getByLabel(ev, "HSCPIsolation", "R01");  //New format used for data since 17-07-2015
  if (!IsolationH01.isValid()) {
    IsolationH01.getByLabel(ev, "HSCPIsolation01");  //Old format used for first 2015B data, Signal and MC Backgrounds
    if (!IsolationH01.isValid()) {
      printf("Invalid IsolationH\n");
      return;
    }
  }
  const ValueMap<HSCPIsolation>& IsolationMap01 = *IsolationH01.product();

  if (!track.isNull()) {
    HSCPIsolation hscpIso05 = IsolationMap05.get((size_t)track.key());
    HSCPIsolation hscpIso03 = IsolationMap03.get((size_t)track.key());
    HSCPIsolation hscpIso01 = IsolationMap01.get((size_t)track.key());
    fprintf(pFile,
            "Isolation05 --> TkCount=%6.2f TkSumEt=%6.2f EcalE/P=%6.2f HcalE/P=%6.2f --> E/P=%6.2f\n",
            hscpIso05.Get_TK_Count(),
            hscpIso05.Get_TK_SumEt(),
            hscpIso05.Get_ECAL_Energy() / track->p(),
            hscpIso05.Get_HCAL_Energy() / track->p(),
            (hscpIso05.Get_ECAL_Energy() + hscpIso05.Get_HCAL_Energy()) / track->p());
    fprintf(pFile,
            "Isolation03 --> TkCount=%6.2f TkSumEt=%6.2f EcalE/P=%6.2f HcalE/P=%6.2f --> E/P=%6.2f\n",
            hscpIso03.Get_TK_Count(),
            hscpIso03.Get_TK_SumEt(),
            hscpIso03.Get_ECAL_Energy() / track->p(),
            hscpIso03.Get_HCAL_Energy() / track->p(),
            (hscpIso03.Get_ECAL_Energy() + hscpIso03.Get_HCAL_Energy()) / track->p());
    fprintf(pFile,
            "Isolation01 --> TkCount=%6.2f TkSumEt=%6.2f EcalE/P=%6.2f HcalE/P=%6.2f --> E/P=%6.2f\n",
            hscpIso01.Get_TK_Count(),
            hscpIso01.Get_TK_SumEt(),
            hscpIso01.Get_ECAL_Energy() / track->p(),
            hscpIso01.Get_HCAL_Energy() / track->p(),
            (hscpIso01.Get_ECAL_Energy() + hscpIso01.Get_HCAL_Energy()) / track->p());
  }
  fprintf(pFile, "\n");
}

void DumpInfo(string DIRNAME, string Pattern, string CutIndexStr = "0", string MassCutStr = "-1", string Region = "D") {
  bool (*IncorrectRegion)(double, double, double, double, double, double);
  if (Region == "A") {
    IncorrectRegion = NotARegion;
  } else if (Region == "B") {
    IncorrectRegion = NotBRegion;
  } else if (Region == "C") {
    IncorrectRegion = NotCRegion;
  } else if (Region == "D") {
    IncorrectRegion = NotDRegion;
  } else if (Region == "E") {
    IncorrectRegion = NotERegion;
  } else if (Region == "F") {
    IncorrectRegion = NotFRegion;
  } else if (Region == "G") {
    IncorrectRegion = NotGRegion;
  } else if (Region == "H") {
    IncorrectRegion = NotHRegion;
  } else {
    Region = "D";
    IncorrectRegion = NotDRegion;
  }

  int CutIndex;
  double MassCut;
  sscanf(CutIndexStr.c_str(), "%d", &CutIndex);
  sscanf(MassCutStr.c_str(), "%lf", &MassCut);

  InitBaseDirectory();
  BaseDirectory = "/storage/data/cms/store/user/jozobec/HSCP2016/";
  GetSampleDefinition(samples, DIRNAME + "/../../AnalysisCode/Analysis_Samples.txt");
  keepOnlySamplesOfTypeX(samples, 0);
  TypeMode = TypeFromPattern(Pattern);
  if (TypeMode == 4)
    useClusterCleaning = false;

  string Data = "Data13TeV16";

  TFile* InputFile = new TFile((Pattern + "/Histos.root").c_str());
  TH1D* HCuts_Pt = (TH1D*)GetObjectFromPath(InputFile, "HCuts_Pt");
  TH1D* HCuts_I = (TH1D*)GetObjectFromPath(InputFile, "HCuts_I");
  TH1D* HCuts_TOF = (TH1D*)GetObjectFromPath(InputFile, "HCuts_TOF");
  TH1D* H_A = (TH1D*)GetObjectFromPath(InputFile, Data + "/H_A");
  TH1D* H_B = (TH1D*)GetObjectFromPath(InputFile, Data + "/H_B");
  TH1D* H_C = (TH1D*)GetObjectFromPath(InputFile, Data + "/H_C");
  TH1D* H_D = (TH1D*)GetObjectFromPath(InputFile, Data + "/H_D");
  TH1D* H_E = (TH1D*)GetObjectFromPath(InputFile, Data + "/H_E");
  TH1D* H_F = (TH1D*)GetObjectFromPath(InputFile, Data + "/H_F");
  TH1D* H_G = (TH1D*)GetObjectFromPath(InputFile, Data + "/H_G");
  TH1D* H_H = (TH1D*)GetObjectFromPath(InputFile, Data + "/H_H");
  TH1D* H_P = (TH1D*)GetObjectFromPath(InputFile, Data + "/H_P");
  PtCut = HCuts_Pt->GetBinContent(CutIndex + 1);
  ICut = HCuts_I->GetBinContent(CutIndex + 1);
  TOFCut = HCuts_TOF->GetBinContent(CutIndex + 1);

  unsigned int CurrentRun = 0;
  moduleGeom::loadGeometry("../../../data/CMS_GeomTree.root");
  muonTimingCalculator tofCalculator;
  tofCalculator.loadTimeOffset("../../../data/MuonTimeOffset.txt");
  bool is2016 = (Data.find("13TeV16") != string::npos);
  if (isData) {
    dEdxSF[0] = 1.00000;
    dEdxSF[1] = (is2016) ? 1.41822 : 1.21836;
    dEdxTemplates = loadDeDxTemplate(is2016 ? "../../../data/Data13TeV16_dEdxTemplate.root"
                                            : "../../../data/Data13TeV_Deco_SiStripDeDxMip_3D_Rcd_v2_CCwCI.root",
                                     true);
  } else {
    dEdxSF[0] = (is2016) ? 1.41822 : 1.09708;
    dEdxSF[1] = (is2016) ? 1.09256 : 1.01875;
    dEdxTemplates = loadDeDxTemplate(is2016 ? "../../../data/MC13TeV16_dEdxTemplate.root"
                                            : "../../../data/MC13TeV_Deco_SiStripDeDxMip_3D_Rcd_v2_CCwCI.root",
                                     true);
  }

  if (isData) {
    trackerCorrector.LoadDeDxCalibration("../../../data/Data13TeVGains_v2.root");
  } else {
    trackerCorrector.TrackerGains = NULL;  //FIXME check gain for MC
  }

  TTree* tree = (TTree*)GetObjectFromPath(InputFile, Data + "/HscpCandidates");
  printf("Tree Entries=%lli\n", tree->GetEntries());
  cout << "Cut Pt " << PtCut << " Cut I " << ICut << " TOFCut " << TOFCut << endl;

  std::vector<string> FileName;
  for (unsigned int s = 0; s < samples.size(); s++) {
    if (samples[s].Name == Data)
      GetInputFiles(samples[s], BaseDirectory, FileName);
  }
  fwlite::ChainEvent ev(FileName);

  unsigned int Run, Event, HscpI;
  float Pt, I, TOF, Mass;

  tree->SetBranchAddress("Run", &Run);
  tree->SetBranchAddress("Event", &Event);
  tree->SetBranchAddress("Hscp", &HscpI);
  tree->SetBranchAddress("Pt", &Pt);
  tree->SetBranchAddress("I", &I);
  tree->SetBranchAddress("TOF", &TOF);
  tree->SetBranchAddress("Mass", &Mass);

  FILE* pFile = fopen("DumpInfo.txt", "w");
  fprintf(pFile, "A = %6.2E\n", H_A->GetBinContent(CutIndex + 1));
  fprintf(pFile, "B = %6.2E\n", H_B->GetBinContent(CutIndex + 1));
  fprintf(pFile, "C = %6.2E\n", H_C->GetBinContent(CutIndex + 1));
  fprintf(pFile, "D = %6.2E\n", H_D->GetBinContent(CutIndex + 1));
  fprintf(pFile, "E = %6.2E\n", H_E->GetBinContent(CutIndex + 1));
  fprintf(pFile, "F = %6.2E\n", H_F->GetBinContent(CutIndex + 1));
  fprintf(pFile, "G = %6.2E\n", H_G->GetBinContent(CutIndex + 1));
  fprintf(pFile, "H = %6.2E\n", H_H->GetBinContent(CutIndex + 1));
  fprintf(pFile, "OBSERVED  EVENTS = %6.2E\n", H_D->GetBinContent(CutIndex + 1));
  fprintf(pFile, "PREDICTED EVENTS = %6.2E+-%6.2E\n", H_P->GetBinContent(CutIndex + 1), H_P->GetBinError(CutIndex + 1));
  FILE* pickEvent = fopen("PickEvent.txt", "w");
  printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Scanning %s                   :", Region.c_str());

  int TreeStep = tree->GetEntries() / 50;
  if (TreeStep == 0)
    TreeStep = 1;
  for (Int_t i = 0; i < tree->GetEntries(); i++) {
    if (i % TreeStep == 0) {
      printf(".");
      fflush(stdout);
    }
    tree->GetEntry(i);
    //      printf("%6i %9i %1i  %6.2f %6.2f %6.2f\n",Run,Event,HscpI,Pt,I,TOF);

    if (IncorrectRegion((double)Pt, PtCut, (double)I, ICut, (double)TOF, TOFCut) || (MassCut > -1 && Mass <= MassCut))
      continue;
    ev.to(Run, Event);
    fprintf(pickEvent, "%i:%i:%i,\n", Run, ev.eventAuxiliary().luminosityBlock(), Event);

    //if run change, update conditions
    if (CurrentRun != ev.eventAuxiliary().run()) {
      CurrentRun = ev.eventAuxiliary().run();
      tofCalculator.setRun(CurrentRun);
      trackerCorrector.setRun(CurrentRun);
    }

    fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
    hscpCollHandle.getByLabel(ev, "HSCParticleProducer");
    if (!hscpCollHandle.isValid()) {
      printf("HSCP Collection NotFound\n");
      continue;
    }
    const susybsm::HSCParticleCollection& hscpColl = *hscpCollHandle;

    susybsm::HSCParticle hscp = hscpColl[HscpI];
    DumpCandidateInfo(hscp, ev, pFile, Mass, tofCalculator);

    for (unsigned int h = 0; h < hscpColl.size(); h++) {
      if (h == HscpI)
        continue;
      reco::MuonRef muon = hscpColl[h].muonRef();
      reco::TrackRef track = hscpColl[h].trackRef();
      if (!track.isNull()) {
        fprintf(pFile,
                "other tracks P=%7.2f  Pt=%7.2f+-%6.2f (Cut=%6.2f) Eta=%+6.2f  Phi=%+6.2f  NOH=%2i\n",
                track->p(),
                track->pt(),
                track->ptError(),
                PtCut,
                track->eta(),
                track->phi(),
                track->found());
      } else {
        fprintf(pFile, "other tracks muontracks\n");
      }
    }
    fflush(pFile);
  }
  printf("\n");
  fclose(pFile);
  fclose(pickEvent);

  /*
   fwlite::ChainEvent treeD(DataFileName);
   SetWeight(-1);
   printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
   printf("Scanning D                   :");
   TreeStep = treeD.size()/50;if(TreeStep==0)TreeStep=1;

   for(Long64_t ientry=0;ientry<treeD.size();ientry++){
      treeD.to(ientry);
      if(MaxEntry>0 && ientry>MaxEntry)break;
      if(ientry%TreeStep==0){printf(".");fflush(stdout);}

      DataPlots.TotalE->Fill(0.0,Event_Weight);  
      if(!PassTrigger(treeD) )continue;
      DataPlots.TotalTE->Fill(0.0,Event_Weight);

      fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
      hscpCollHandle.getByLabel(treeD,"HSCParticleProducer");
      if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
      const susybsm::HSCParticleCollection& hscpColl = *hscpCollHandle;

      fwlite::Handle<DeDxDataValueMap> dEdxSCollH;
      dEdxSCollH.getByLabel(treeD, dEdxS_Label.c_str());
      if(!dEdxSCollH.isValid()){printf("Invalid dEdx Selection collection\n");continue;}

      fwlite::Handle<DeDxDataValueMap> dEdxMCollH;
      dEdxMCollH.getByLabel(treeD, dEdxM_Label.c_str());
      if(!dEdxMCollH.isValid()){printf("Invalid dEdx Mass collection\n");continue;}

      fwlite::Handle<MuonTimeExtraMap> TOFCollH;
      TOFCollH.getByLabel(treeD, "muontiming",TOF_Label.c_str());
      if(!TOFCollH.isValid()){printf("Invalid TOF collection\n");return;}
      
      bool* HSCPTk = new bool[PtCut.size()]; for(unsigned int CutIndex=0;CutIndex<PtCut.size();CutIndex++){  HSCPTk[CutIndex] = false;   }
      for(unsigned int c=0;c<hscpColl.size();c++){
         susybsm::HSCParticle hscp  = hscpColl[c];
         reco::MuonRef  muon  = hscp.muonRef();
         reco::TrackRef track = hscp.trackRef();
         if(track.isNull())continue;

         const DeDxData& dedxSObj  = dEdxSCollH->get(track.key());
         const DeDxData& dedxMObj  = dEdxMCollH->get(track.key());
         const reco::MuonTimeExtra* tof = NULL;
        if(TypeMode==2 && !hscp.muonRef().isNull()){ tof  = &TOFCollH->get(hscp.muonRef().key()); }


         double MuonTOF = GlobalMinTOF;
         if(tof){MuonTOF = tof->inverseBeta(); }
         if(track->pt()>40 && Mass>75)stPlots_FillTree(DataPlots, treeD.eventAuxiliary().run(),treeD.eventAuxiliary().event(), c, track->pt(), dedxSObj.dEdx(), tof ? tof->inverseBeta() : -1);
      } // end of Track Loop
      for(unsigned int CutIndex=0;CutIndex<PtCut.size();CutIndex++){  if(HSCPTk[CutIndex]){DataPlots.HSCPE->Fill(CutIndex,Event_Weight); }  }
   }// end of Event Loop
   //stPlots_CloseTree(DataPlots);
   printf("\n");
*/
}

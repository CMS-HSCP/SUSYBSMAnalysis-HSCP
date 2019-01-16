#include "SUSYBSMAnalysis/HSCP/interface/HSCPHelpers.h"

namespace hscphelpers {
  TH3F* preprocess_dEdx_template(TH3F* DeDxMap_, bool splitByModuleType) {
    TH3F* Prob_ChargePath  = (TH3F*)(DeDxMap_->Clone("Prob_ChargePath")); 
    Prob_ChargePath->Reset();
    Prob_ChargePath->SetDirectory(0); 

    if(!splitByModuleType){
      Prob_ChargePath->RebinX(Prob_ChargePath->GetNbinsX()-1); // <-- do not include pixel in the inclusive
      // ^ this part is shady, FIXME
    }

    for(int i=0;i<=Prob_ChargePath->GetXaxis()->GetNbins()+1;i++){
      for(int j=0;j<=Prob_ChargePath->GetYaxis()->GetNbins()+1;j++){
	double Ni = 0;
	for(int k=0;k<=Prob_ChargePath->GetZaxis()->GetNbins()+1;k++){Ni+=DeDxMap_->GetBinContent(i,j,k);}

	for(int k=0;k<=Prob_ChargePath->GetZaxis()->GetNbins()+1;k++){
	  double tmp = 0;
	  for(int l=0;l<=k;l++){ tmp+=DeDxMap_->GetBinContent(i,j,l);}

	  if(Ni>0){
	    Prob_ChargePath->SetBinContent (i, j, k, tmp/Ni);
	  }else{
	    Prob_ChargePath->SetBinContent (i, j, k, 0);
	  }
	}
      }
    }
    return Prob_ChargePath;
  }

  int  muonStations(const reco::HitPattern& hitPattern) {
    int stations[4] = { 0,0,0,0 };
    for (int i=0; i<hitPattern.numberOfValidTrackerHits(); i++) {
      uint32_t pattern = hitPattern.getHitPattern(reco::HitPattern::HitCategory::TRACK_HITS, i );
      if(pattern == 0) break;
      if(hitPattern.muonHitFilter(pattern) && (int(hitPattern.getSubStructure(pattern)) == 1 || int(hitPattern.getSubStructure(pattern)) == 2) && hitPattern.getHitType(pattern) == 0){
	stations[hitPattern.getMuonStation(pattern)-1] = 1;
      }
    }
    return stations[0]+stations[1]+stations[2]+stations[3];

  }



  double rescaledPt(const double& pt, const double& eta, const double& phi, const int& charge, const int& TypeMode)
  {
    if(TypeMode!=3) {
      double newInvPt = 1/pt+0.000236-0.000135*pow(eta,2)+charge*0.000282*TMath::Sin(phi-1.337);
      return 1/newInvPt;
    }
    else {
      double newInvPt = (1./pt)*1.1;
      return 1/newInvPt;
    }
  }





}

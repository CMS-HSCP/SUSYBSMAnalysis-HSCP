#ifndef HIPTrackLossEmulator_h
#define HIPTrackLossEmulator_h

#include <vector>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <cstdlib>

class HIPTrackLossEmulator{
  private:
     TH1D*  h;
     double lossRate;
  public:
     HIPTrackLossEmulator(){
        double NVertexBins[] = {0, 5, 10, 15, 22, 28, 35, 9999};
        h = new TH1D("TrackLoss_vs_NVertex", "TrackLoss_vs_NVertex", static_cast<unsigned int> (sizeof(NVertexBins)/sizeof(double)) -1, NVertexBins);
        h->SetBinContent(h->FindBin( 2.5), 1.000);
        h->SetBinContent(h->FindBin( 7.5), 0.998);
        h->SetBinContent(h->FindBin(12.5), 0.995);
        h->SetBinContent(h->FindBin(17.5), 0.992);
        h->SetBinContent(h->FindBin(24.5), 0.985);
        h->SetBinContent(h->FindBin(30.0), 0.980);
        h->SetBinContent(h->FindBin(37.0), 0.970);

        lossRate = 0.0;
     }

     ~HIPTrackLossEmulator(){
//        delete h;
     }

     void SetHIPTrackLossRate(const std::vector<reco::Vertex>& vertexColl){
        if(vertexColl.size()<1){lossRate = 1.0; return;}
        else lossRate = 1.0 - h->GetBinContent(h->FindBin(vertexColl.size()));
     }

     bool TrackSurvivesHIPInefficiency(){
        return (((rand()%999999)*1.0/1000000) > lossRate);
     }
};

#endif

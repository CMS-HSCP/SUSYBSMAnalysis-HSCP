B#include "FWCore/Framework/interface/Event.h"
#include "SUSYBSMAnalysis/HSCP/interface/BigNtuple.h"
#include "TTree.h"

void BigNtuple::set_evtinfo(TTree* tree) {
  tree->Branch("run" , &run_, "run/i");
  tree->Branch("lumi", &lumi_, "lumi/i");
  tree->Branch("evt" , &evt_, "evt/i");
}

void BigNtuple::fill_evtinfo(const edm::EventID& id) {
  lumi_ = id.luminosityBlock();
  run_  = id.run();
  evt_  = id.event();
}

void BigNtuple::set_trackinfo(TTree* tree){

  tree->Branch("eta", &eta_, "eta/f");
  tree->Branch("dz", &dz_, "dz/f"); 
  tree->Branch("dxy", &dxy_, "dxy/f"); 
  tree->Branch("track_Is", &track_Is_, "track_Is/f");
  tree->Branch("track_Ih", &track_Ih_, "track_Ih/f"); 
  tree->Branch("track_TOF", &track_TOF_, "track_TOF/f"); 

}



void BigNtuple::fill_trackinfo(const reco::TrackRef& trk, reco::Vertex& bestVtx, ){

  eta_ = trk->eta();
  dz_  = trk->dz(bestVtx.position());
  dxy_ = trk->dxy(bestVtx.position());

  if(dedxSObj){
    track_Is_ = dedxSObj->dEdx();
    track_Ih_  = dedxMObj->dEdx();
    if (tof){
      track_TOF = tof->inverseBeta();
    }
  }
}



void BigNtuple::set_trigInfo(TTree* tree){



}


void BigNtuple::fill_trigInfo(const edm::TriggerResults& triggerResults, const edm::TriggerNames& trigNames){


}

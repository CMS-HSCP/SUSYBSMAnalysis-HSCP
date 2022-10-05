#include "SUSYBSMAnalysis/Analyzer/interface/TrigToolsFuncs.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"

void trigtools::getP4sOfObsPassingFilter(std::vector<math::XYZTLorentzVector>& p4s,const trigger::TriggerEvent& trigEvent,const std::string& filterName,const std::string& hltProcess)
{
  p4s.clear();

  edm::InputTag filterTag(filterName,"",hltProcess); 
  trigger::size_type filterIndex = trigEvent.filterIndex(filterTag); 
  if(filterIndex<trigEvent.sizeFilters()){ //check that filter is in triggerEvent
    const trigger::Keys& trigKeys = trigEvent.filterKeys(filterIndex); 
    const trigger::TriggerObjectCollection & trigObjColl(trigEvent.getObjects());
    for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
      const trigger::TriggerObject& obj = trigObjColl[*keyIt];
      math::XYZTLorentzVector objP4;
      objP4.SetPxPyPzE(obj.px(),obj.py(),obj.pz(),obj.energy());
      p4s.push_back(objP4);
    }//end loop over keys
  }//end check that filter is valid and in trigEvent
}



void trigtools::getP4sOfObsPassingFilter(std::vector<TLorentzVector>& p4s,const trigger::TriggerEvent& trigEvent,const std::string& filterName,const std::string& hltProcess)
{
  p4s.clear();
 
  edm::InputTag filterTag(filterName,"",hltProcess); 
  trigger::size_type filterIndex = trigEvent.filterIndex(filterTag); 
  if(filterIndex<trigEvent.sizeFilters()){ //check that filter is in triggerEvent
    const trigger::Keys& trigKeys = trigEvent.filterKeys(filterIndex); 
    const trigger::TriggerObjectCollection & trigObjColl(trigEvent.getObjects());
    for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
      const trigger::TriggerObject& obj = trigObjColl[*keyIt];
      TLorentzVector objP4;
      objP4.SetPtEtaPhiM(obj.pt(),obj.eta(),obj.phi(),obj.mass());
      p4s.push_back(objP4);
    }//end loop over keys
  }//end check that filter is valid and in trigEvent
}

void trigtools::dumpTriggerEvent(const trigger::TriggerEvent& trigEvt)
{
  std::cout <<"number of filters in event "<<trigEvt.sizeFilters()<<std::endl;
  for(size_t filterNr=0;filterNr<trigEvt.sizeFilters();filterNr++){
    const std::string filterName(trigEvt.filterTag(filterNr).label());  
    const trigger::Keys& trigKeys = trigEvt.filterKeys(filterNr);//trigger::Keys is actually a vector<uint16_t> holding the position of trigger objects in the trigger collection passing the filter
    std::cout <<"filter "<<filterName<<" has "<<trigKeys.size()<<" passing "<<std::endl;
    const trigger::TriggerObjectCollection & trigObjColl(trigEvt.getObjects());
    for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){
      const trigger::TriggerObject& obj = trigObjColl[*keyIt];
      TLorentzVector p4;
      p4.SetPtEtaPhiM(obj.pt(),obj.eta(),obj.phi(),obj.mass());
    }
  } 

  
}


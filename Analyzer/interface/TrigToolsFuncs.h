#ifndef RHABERLE_TRIGTOOLSFUNCS
#define RHABERLE_TRIGTOOLSFUNCS

#include "DataFormats/Math/interface/LorentzVector.h"

#include "TLorentzVector.h"

#include<vector>
#include<string>


namespace trigger{
  class TriggerEvent;
}

namespace trigtools {
  
  //this function takes the trigger event, a filtername + hlt process name (usally HLT but different if 
  //  //HLT was re-run) and returns a vector of the four-momenta of all objects passing the filter
  //
  void getP4sOfObsPassingFilter(std::vector<math::XYZTLorentzVector>& p4s,const trigger::TriggerEvent& trigEvent,const std::string& filterName,const std::string& hltProcess="HLT");
  
  //a TLorentzVector version for my ntuplising needs
  void getP4sOfObsPassingFilter(std::vector<TLorentzVector>& p4s,const trigger::TriggerEvent& trigEvent,const std::string& filterName,const std::string& hltProcess="HLT");
            
  void dumpTriggerEvent(const trigger::TriggerEvent& trigEvt);

  bool passedFilter(const trigger::TriggerEvent& trigEvt,const std::string& givenFilter);
}
  
#endif

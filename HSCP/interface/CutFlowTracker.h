#ifndef CutFlowTracker_h
#define CutFlowTracker_h

#include <string>
#include <map>
#include "TFile.h"

class CutFlowTracker {
 public:
 CutFlowTracker():
  cutflow_(),
    npoints_(),
    active_(true),
    verbose_(false)
      {}
  void track(std::string point, float weight=1, std::string groupname="");
  void writeTo(TDirectory &file);
  void writeTo(TDirectory &file, std::string group);
  void activate() {active_ = true;};
  void deactivate() {active_ = false;};
  void group(std::string group) {group_ = group;}

 private:
  std::map<std::string, std::pair<size_t, float> > cutflow_;
  size_t npoints_;
  bool active_, verbose_;
  std::string group_;
};

#endif

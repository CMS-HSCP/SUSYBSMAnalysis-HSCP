#include "SUSYBSMAnalysis/HSCP/interface/CutFlowTracker.h"
#include "TH1F.h"
#include <boost/algorithm/string/predicate.hpp>
#include <tuple>
#include <vector>
#include <algorithm>

void CutFlowTracker::track(std::string label, float weight, std::string groupname)
{
  if(groupname.size() == 0 && group_.size() != 0) groupname = group_;
  std::string pointname = (groupname.size() != 0) ? groupname+"/"+label : label;
  if(!active_) return;
  auto point = cutflow_.find(pointname);
  if(point != cutflow_.end()){
    point->second.second += weight;
  }
  else {
    cutflow_.insert( 
		    std::make_pair(
				   pointname,
				   std::make_pair(
						  npoints_,
						  weight
						  )
				   )
		     );
    npoints_++;
  }
}

void CutFlowTracker::writeTo(TDirectory &file)
{
  file.cd();
  TH1F histo("cut_flow", "cut_flow", npoints_, 0, npoints_);
  TAxis *xax = histo.GetXaxis();
  for(auto& point : cutflow_){
    xax->SetBinLabel(point.second.first+1, point.first.c_str());
    histo.SetBinContent(point.second.first+1, point.second.second);
  }
  histo.Write();
}

void CutFlowTracker::writeTo(TDirectory &file, std::string group)
{
  file.cd();
  size_t npts = 0;
  std::vector< std::tuple<std::string, size_t, float> > bins;
  for(auto& entry : cutflow_) {
    if(boost::starts_with(entry.first, group)) {
      npts++;
      std::string name = entry.first.substr(group.size()+1,entry.first.size()-group.size()-1);
      size_t idx = entry.second.first;
      float  val = entry.second.second;
      bins.push_back(
		     std::make_tuple(name, idx, val)
		     );
    }
  }
  std::sort(bins.begin(), bins.end(), [](auto const &t1, auto const &t2) { //black compiler magic!
      return std::get<1>(t1) < std::get<1>(t2); 
    });

  TH1F histo("cut_flow", "cut_flow", npts, 0, npts);
  TAxis *xax = histo.GetXaxis();
  for(size_t i = 0; i < npts; ++i){
    xax->SetBinLabel(i+1, std::get<0>(bins[i]).c_str());
    histo.SetBinContent(i+1, std::get<2>(bins[i]));
  }
  histo.Write();
}

#include "RecoTotemRP/RPRomanPotResolutionService/interface/RPFitResolution.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RP2DHit.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"


#include <iostream>

RPFitResolution::RPFitResolution(const edm::ParameterSet& conf)
{
  strip_alignment_res_degradation_ = conf.getParameter<double>("StripAlignmentResolutionDegradation");
  verbosity_ = conf.getParameter<int>("Verbosity");
  var_degrad_ = strip_alignment_res_degradation_*strip_alignment_res_degradation_;
}


RP2DHit RPFitResolution::Create2DHit(unsigned int rp_id, const TotemRPLocalTrack &track)
{
  RP2DHit tmp(track.getX0(), track.getY0(), track.getX0Variance()*var_degrad_, 
      track.getY0Variance()*var_degrad_, track.getZ0());

  if(verbosity_)
  {
    std::cout<<"RPFitResolution::Create2DHit"<<std::endl;
    std::cout<<"rp_id="<<rp_id<<" pos=("<<track.getX0()<<","<<track.getY0()<<","<<track.getZ0()<<")"<<std::endl;
    std::cout<<"\t sigma x="<<track.getX0Sigma()<<" sigma y="<<track.getY0Sigma()<<std::endl;
    std::cout<<"\t oryg variance x="<<track.getX0Variance()<<" sigma y="<<track.getY0Variance()<<std::endl;
    std::cout<<"\t degrad var x="<<track.getX0Variance()*var_degrad_<<" sigma y="<<track.getY0Variance()*var_degrad_<<std::endl;
    std::cout<<"\tvar_degrad_="<<var_degrad_<<std::endl;
    std::cout<<std::endl;
    std::cout<<tmp<<std::endl;
  }
  return tmp;
}


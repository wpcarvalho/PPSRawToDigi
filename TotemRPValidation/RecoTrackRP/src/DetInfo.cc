#include "TotemRPValidation/RecoTrackRP/interface/DetInfo.h"

DetInfo::DetInfo(const std::string &path, RPDetId det_id, const edm::ParameterSet& conf)
 : BaseHistogramManager(path), det_id_(det_id)
{
  InitializeHistograms();
}


void DetInfo::InitializeHistograms()
{
  char name[1024];
  
  sprintf(name, "cluster_size_%04i", det_id_);
  cluster_size_ = TH1F(name, name, 11, -0.5, 10.5);
  cluster_size_.SetDirectory(0);
  RegisterHistogram(cluster_size_);
  
  sprintf(name, "cluster_multiplicity_%04i", det_id_);
  cluster_multiplicity_ = TH1F(name, name, 21, -0.5, 20.5);
  cluster_multiplicity_.SetDirectory(0);
  RegisterHistogram(cluster_multiplicity_);
  
  sprintf(name, "strip_profile_%04i", det_id_);
  strip_profile_ = TH1F(name, name, 513, -0.5, 512.5);
  strip_profile_.SetDirectory(0);
  RegisterHistogram(strip_profile_);
  
  sprintf(name, "track_residual_%04i", det_id_);
  track_residual_ = TH1F(name, name, 300, -0.3, 0.3); ///<[mm]
  track_residual_.SetDirectory(0);
  RegisterHistogram(track_residual_);
}


void DetInfo::FillHistograms(const clu_col_type &col, double residual)
{
  clu_col_type::const_iterator begin = col.begin();
  clu_col_type::const_iterator end = col.end();
  
  cluster_multiplicity_.Fill(col.size());
  track_residual_.Fill(residual);
  
  for(clu_col_type::const_iterator it = begin; it!=end; ++it)  ///<All clusters in the detector plane are filled in. 
  {
    cluster_size_.Fill(it->GetNumberOfStrips());
    strip_profile_.Fill(it->CentreStripPos());
  }
}

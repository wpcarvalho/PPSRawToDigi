#include "TotemRPValidation/RPReconstructedTracksValidation/interface/DetInfo.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

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
  cluster_size_.SetXTitle("Strips per cluster #");
  cluster_size_.SetYTitle("Entries");
  
  sprintf(name, "cluster_multiplicity_%04i", det_id_);
  cluster_multiplicity_ = TH1F(name, name, 21, -0.5, 20.5);
  cluster_multiplicity_.SetDirectory(0);
  RegisterHistogram(cluster_multiplicity_);
  cluster_multiplicity_.SetXTitle("Cluster multiplicity");
  cluster_multiplicity_.SetYTitle("Entries");
  
  sprintf(name, "strip_profile_%04i", det_id_);
  strip_profile_ = TH1F(name, name, 513, -0.5, 512.5);
  strip_profile_.SetDirectory(0);
  RegisterHistogram(strip_profile_);
  strip_profile_.SetXTitle("Cluster position [Strip number]");
  strip_profile_.SetYTitle("Entries");
  
  sprintf(name, "track_residual_%04i", det_id_);
  track_residual_ = TH1F(name, name, 300, -0.3, 0.3); ///<[mm]
  track_residual_.SetDirectory(0);
  RegisterHistogram(track_residual_);
  track_residual_.SetXTitle("Hit point residual [mm]");
  track_residual_.SetYTitle("Entries");
  
  sprintf(name, "trigger_sectors_%04i", det_id_);
  trigger_sectors_ = TH1F(name, name, 16, -0.5, 15.5); ///<trigger sector number
  trigger_sectors_.SetDirectory(0);
  RegisterHistogram(trigger_sectors_);
  trigger_sectors_.SetXTitle("Trigger sector");
  trigger_sectors_.SetYTitle("Entries");
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


void DetInfo::FillTriggerSector(int trigger_sector)
{
  trigger_sectors_.Fill(trigger_sector);
}

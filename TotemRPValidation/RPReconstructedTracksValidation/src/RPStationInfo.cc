#include "TotemRPValidation/RPReconstructedTracksValidation/interface/RPStationInfo.h"

RPStationInfo::RPStationInfo(const std::string &path, RPStationId rp_station_id, const edm::ParameterSet& conf)
 : BaseHistogramManager(path), rp_station_id_(rp_station_id)
{
  InitializeHistograms();
}


void RPStationInfo::InitializeHistograms()
{
  char name[1024];
  
  sprintf(name, "log10_t_dist_%04i", rp_station_id_);
  log10_t_dist_ = TH1F(name, name, 840, -4, 10);
  log10_t_dist_.SetDirectory(0);
  RegisterHistogram(log10_t_dist_);
  log10_t_dist_.SetXTitle("Log(-t/GeV^{2})");
  log10_t_dist_.SetYTitle("Entries");
  
  sprintf(name, "log10_t_dist_seen_tracks_%04i", rp_station_id_);
  log10_t_dist_seen_tracks_ = TH1F(name, name, 840, -4, 10);
  log10_t_dist_seen_tracks_.SetDirectory(0);
  RegisterHistogram(log10_t_dist_seen_tracks_);
  log10_t_dist_seen_tracks_.SetXTitle("Log(-t/GeV^{2})");
  log10_t_dist_seen_tracks_.SetYTitle("Entries");
  
  sprintf(name, "log10_t_dist_acceptance_%04i", rp_station_id_);
  log10_t_dist_acceptance_ = TH1F(name, name, 840, -4, 10);
  log10_t_dist_acceptance_.SetDirectory(0);
  RegisterHistogram(log10_t_dist_acceptance_);
  log10_t_dist_acceptance_.SetXTitle("Log(-t/GeV^{2})");
  log10_t_dist_acceptance_.SetYTitle("Acceptance");
  
  sprintf(name, "log_ksi_log_t_distribution_%04i", rp_station_id_);
  log_ksi_log_t_distribution_ = TH2F(name, name, 50, -4, 1, 35, -3.5, 0.);
  log_ksi_log_t_distribution_.SetDirectory(0);
  RegisterHistogram(log_ksi_log_t_distribution_);
  log_ksi_log_t_distribution_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_distribution_.SetYTitle("Log(-#xi)");
  
  sprintf(name, "log_ksi_log_t_seen_tracks_%04i", rp_station_id_);
  log_ksi_log_t_seen_tracks_ = TH2F(name, name, 50, -4, 1, 35, -3.5, 0.);
  log_ksi_log_t_seen_tracks_.SetDirectory(0);
  RegisterHistogram(log_ksi_log_t_seen_tracks_);
  log_ksi_log_t_seen_tracks_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_seen_tracks_.SetYTitle("Log(-#xi)");
  
  sprintf(name, "log_ksi_log_t_acceptance_%04i", rp_station_id_);
  log_ksi_log_t_acceptance_ = TH2F(name, name, 50, -4, 1, 35, -3.5, 0.);
  log_ksi_log_t_acceptance_.SetDirectory(0);
  RegisterHistogram(log_ksi_log_t_acceptance_);
  log_ksi_log_t_acceptance_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_acceptance_.SetYTitle("Log(-#xi)");
  
  sprintf(name, "track_distribution_xy_centre_%04i", rp_station_id_);
  track_distribution_xy_centre_ = TH2F(name, name, 80, -40, 40, 80, -40, 40);
  //track_distribution_xy_centre_.SetBit(TH1::kCanRebin);
  track_distribution_xy_centre_.SetDirectory(0);
  RegisterHistogram(track_distribution_xy_centre_);
  track_distribution_xy_centre_.SetXTitle("x [mm]");
  track_distribution_xy_centre_.SetYTitle("y [mm]");
   
  sprintf(name, "track_distribution_x_centre_%04i", rp_station_id_);
  track_distribution_x_centre_ = TH1F(name, name, 400, -40, 40);
  //track_distribution_x_centre_.SetBit(TH1::kCanRebin);
  track_distribution_x_centre_.SetDirectory(0);
  RegisterHistogram(track_distribution_x_centre_); 
  track_distribution_x_centre_.SetXTitle("x [mm]");
  track_distribution_x_centre_.SetYTitle("Entries");
   
  sprintf(name, "track_distribution_y_centre_%04i", rp_station_id_);
  track_distribution_y_centre_ = TH1F(name, name, 400, -40, 40);
  //track_distribution_y_centre_.SetBit(TH1::kCanRebin);
  track_distribution_y_centre_.SetDirectory(0);
  RegisterHistogram(track_distribution_y_centre_);
  track_distribution_y_centre_.SetXTitle("y [mm]");
  track_distribution_y_centre_.SetYTitle("Entries");
   
  sprintf(name, "track_distribution_thx_centre_%04i", rp_station_id_);
  track_distribution_thx_centre_ = TH1F(name, name, 400, -0.0001, 0.0001);
  track_distribution_thx_centre_.SetBit(TH1::kCanRebin);
  track_distribution_thx_centre_.SetDirectory(0);
  RegisterHistogram(track_distribution_thx_centre_); 
  track_distribution_thx_centre_.SetXTitle("#Theta_{x} [rad]");
  track_distribution_thx_centre_.SetYTitle("Entries");
   
  sprintf(name, "track_distribution_thy_centre_%04i", rp_station_id_);
  track_distribution_thy_centre_ = TH1F(name, name, 400, -0.0001, 0.0001);
  track_distribution_thy_centre_.SetBit(TH1::kCanRebin);
  track_distribution_thy_centre_.SetDirectory(0);
  RegisterHistogram(track_distribution_thy_centre_); 
  track_distribution_thy_centre_.SetXTitle("#Theta_{y} [rad]");
  track_distribution_thy_centre_.SetYTitle("Entries");
   
  sprintf(name, "track_dist_thx_vs_x_centre_%04i", rp_station_id_);
  track_dist_thx_vs_x_centre_ = TH2F(name, name, 300, -10, 10, 300, -0.0001, 0.0001);
  track_dist_thx_vs_x_centre_.SetBit(TH1::kCanRebin);
  track_dist_thx_vs_x_centre_.SetDirectory(0);
  RegisterHistogram(track_dist_thx_vs_x_centre_); 
  track_dist_thx_vs_x_centre_.SetXTitle("x [mm]");
  track_dist_thx_vs_x_centre_.SetYTitle("#Theta_{x} [rad]");
   
  sprintf(name, "track_dist_thy_vs_y_centre_%04i", rp_station_id_);
  track_dist_thy_vs_y_centre_ = TH2F(name, name, 300, -10, 10, 300, -0.0001, 0.0001);
  track_dist_thy_vs_y_centre_.SetBit(TH1::kCanRebin);
  track_dist_thy_vs_y_centre_.SetDirectory(0);
  RegisterHistogram(track_dist_thy_vs_y_centre_); 
  track_dist_thy_vs_y_centre_.SetXTitle("y [mm]");
  track_dist_thy_vs_y_centre_.SetYTitle("#Theta_{y} [rad]");
   
  sprintf(name, "track_dist_thy_vs_thx_centre_%04i", rp_station_id_);
  track_dist_thy_vs_thx_centre_ = TH2F(name, name, 300, -0.0001, 0.0001, 300, -0.0001, 0.0001);
  track_dist_thy_vs_thx_centre_.SetBit(TH1::kCanRebin);
  track_dist_thy_vs_thx_centre_.SetDirectory(0);
  RegisterHistogram(track_dist_thy_vs_thx_centre_); 
  track_dist_thy_vs_thx_centre_.SetXTitle("#Theta_{x} [rad]");
  track_dist_thy_vs_thx_centre_.SetYTitle("#Theta_{y} [rad]");
}

void RPStationInfo::Finalize()
{
  log10_t_dist_acceptance_.Divide(&log10_t_dist_seen_tracks_, &log10_t_dist_);
  log_ksi_log_t_acceptance_.Divide(&log_ksi_log_t_seen_tracks_, &log_ksi_log_t_distribution_);
}


void RPStationInfo::Fill_t_Dist(double t)
{
  log10_t_dist_.Fill(TMath::Log10(-t));
}


void RPStationInfo::FillTrackSeen(double t)
{
  log10_t_dist_seen_tracks_.Fill(TMath::Log10(-t));
}


void RPStationInfo::Fill_ksi_t_Dist(double t, double ksi)
{
  log_ksi_log_t_distribution_.Fill(TMath::Log10(-t), TMath::Log10(-ksi));
}


void RPStationInfo::FillTrackSeen(double t, double ksi)
{
  log_ksi_log_t_seen_tracks_.Fill(TMath::Log10(-t), TMath::Log10(-ksi));
}

//GetStationCentreZPosition should be used
void RPStationInfo::FillBasicTrackInfo(const RP2DHit& track_in, const RP2DHit& track_out, double station_z_centre)
{  
  double tx = (track_out.X() - track_in.X())/(track_out.Z() - track_in.Z());
  double x_centre = tx*(station_z_centre-track_in.Z()) + track_in.X();
  
  double ty = (track_out.Y() - track_in.Y())/(track_out.Z() - track_in.Z());
  double y_centre = ty*(station_z_centre-track_in.Z()) + track_in.Y();
  
  double thx = TMath::ATan(tx);
  double thy = TMath::ATan(ty);
  
  track_distribution_xy_centre_.Fill(x_centre, y_centre);
  
  track_distribution_x_centre_.Fill(x_centre);
  track_distribution_y_centre_.Fill(y_centre);
  
  track_distribution_thx_centre_.Fill(thx);
  track_distribution_thy_centre_.Fill(thy);
  
  track_dist_thx_vs_x_centre_.Fill(x_centre, thx);
  track_dist_thy_vs_y_centre_.Fill(y_centre, thy);
  track_dist_thy_vs_thx_centre_.Fill(thx, thy);
}


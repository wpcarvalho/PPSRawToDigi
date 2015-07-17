#include "TotemRPValidation/RecoTrackRP/interface/RPTrackInfo.h"

RPTrackInfo::RPTrackInfo(const std::string &path, RPId rp_id, const edm::ParameterSet&)
 : BaseHistogramManager(path), rp_id_(rp_id)
{
  InitializeHistograms();
}


void RPTrackInfo::InitializeHistograms()
{
  char name[1024];
  
  sprintf(name, "log10_t_dist_%04i", rp_id_);
  log10_t_dist_ = TH1F(name, name, 840, -4, 10);
  log10_t_dist_.SetDirectory(0);
  RegisterHistogram(log10_t_dist_);
  
  sprintf(name, "log10_t_dist_seen_tracks_%04i", rp_id_);
  log10_t_dist_seen_tracks_ = TH1F(name, name, 840, -4, 10);
  log10_t_dist_seen_tracks_.SetDirectory(0);
  RegisterHistogram(log10_t_dist_seen_tracks_);
  
  sprintf(name, "log10_t_dist_acceptance_%04i", rp_id_);
  log10_t_dist_acceptance_ = TH1F(name, name, 840, -4, 10);
  log10_t_dist_acceptance_.SetDirectory(0);
  RegisterHistogram(log10_t_dist_acceptance_);
  
  sprintf(name, "log_ksi_log_t_distribution_%04i", rp_id_);
  log_ksi_log_t_distribution_ = TH2F(name, name, 50, -4, 1, 35, -3.5, 0.);
  log_ksi_log_t_distribution_.SetDirectory(0);
  RegisterHistogram(log_ksi_log_t_distribution_);
  
  sprintf(name, "log_ksi_log_t_seen_tracks_%04i", rp_id_);
  log_ksi_log_t_seen_tracks_ = TH2F(name, name, 50, -4, 1, 35, -3.5, 0.);
  log_ksi_log_t_seen_tracks_.SetDirectory(0);
  RegisterHistogram(log_ksi_log_t_seen_tracks_);
  
  sprintf(name, "log_ksi_log_t_acceptance_%04i", rp_id_);
  log_ksi_log_t_acceptance_ = TH2F(name, name, 50, -4, 1, 35, -3.5, 0.);
  log_ksi_log_t_acceptance_.SetDirectory(0);
  RegisterHistogram(log_ksi_log_t_acceptance_);
  
  sprintf(name, "position_residuum_x_%04i", rp_id_);
  position_residuum_x_ = TH1F(name, name, 2000, -1, 1);
  //position_residuum_x_.SetBit(TH1::kCanRebin);
  position_residuum_x_.SetDirectory(0);
  RegisterHistogram(position_residuum_x_);
  
  sprintf(name, "position_residuum_y_%04i", rp_id_);
  position_residuum_y_ = TH1F(name, name, 2000, -1, 1);
  //position_residuum_y_.SetBit(TH1::kCanRebin);
  position_residuum_y_.SetDirectory(0);
  RegisterHistogram(position_residuum_y_);
  
  sprintf(name, "position_residuum_x_if_exactly_2_rps_in_150_station_%04i", rp_id_);
  position_residuum_x_if_exactly_2_rps_in_150_station_ = TH1F(name, name, 2000, -1, 1);
  //position_residuum_x_if_exactly_2_rps_in_150_station_.SetBit(TH1::kCanRebin);
  position_residuum_x_if_exactly_2_rps_in_150_station_.SetDirectory(0);
  RegisterHistogram(position_residuum_x_if_exactly_2_rps_in_150_station_);
  
  sprintf(name, "position_residuum_y_if_exactly_2_rps_in_150_station_%04i", rp_id_);
  position_residuum_y_if_exactly_2_rps_in_150_station_ = TH1F(name, name, 2000, -1, 1);
  //position_residuum_y_if_exactly_2_rps_in_150_station_.SetBit(TH1::kCanRebin);
  position_residuum_y_if_exactly_2_rps_in_150_station_.SetDirectory(0);
  RegisterHistogram(position_residuum_y_if_exactly_2_rps_in_150_station_);
  
  sprintf(name, "position_pool_x_%04i", rp_id_);
  position_pool_x_ = TH1F(name, name, 4000, -30, 30);
  //position_pool_x_.SetBit(TH1::kCanRebin);
  position_pool_x_.SetDirectory(0);
  RegisterHistogram(position_pool_x_); 
   
  sprintf(name, "position_pool_y_%04i", rp_id_);
  position_pool_y_ = TH1F(name, name, 4000, -30, 30);
  //position_pool_y_.SetBit(TH1::kCanRebin);
  position_pool_y_.SetDirectory(0);
  RegisterHistogram(position_pool_y_); 
   
  sprintf(name, "track_distribution_xy_%04i", rp_id_);
  track_distribution_xy_ = TH2F(name, name, 160, -40, 40, 160, -40, 40);
  //track_distribution_xy_.SetBit(TH1::kCanRebin);
  track_distribution_xy_.SetDirectory(0);
  RegisterHistogram(track_distribution_xy_);
   
  sprintf(name, "track_distribution_x_%04i", rp_id_);
  track_distribution_x_ = TH1F(name, name, 400, -40, 40);
  //track_distribution_x_.SetBit(TH1::kCanRebin);
  track_distribution_x_.SetDirectory(0);
  RegisterHistogram(track_distribution_x_); 
   
  sprintf(name, "track_distribution_y_%04i", rp_id_);
  track_distribution_y_ = TH1F(name, name, 400, -40, 40);
  //track_distribution_y_.SetBit(TH1::kCanRebin);
  track_distribution_y_.SetDirectory(0);
  RegisterHistogram(track_distribution_y_);
   
  sprintf(name, "track_distribution_thx_%04i", rp_id_);
  track_distribution_thx_ = TH1F(name, name, 400, -0.0001, 0.0001);
  track_distribution_thx_.SetBit(TH1::kCanRebin);
  track_distribution_thx_.SetDirectory(0);
  RegisterHistogram(track_distribution_thx_); 
   
  sprintf(name, "track_distribution_thy_%04i", rp_id_);
  track_distribution_thy_ = TH1F(name, name, 400, -0.0001, 0.0001);
  track_distribution_thy_.SetBit(TH1::kCanRebin);
  track_distribution_thy_.SetDirectory(0);
  RegisterHistogram(track_distribution_thy_); 
   
  sprintf(name, "track_reco_hit_no_%04i", rp_id_);
  track_reco_hit_no_ = TH1F(name, name, 11, -0.5, 10.5);
  //track_reco_hit_no_.SetBit(TH1::kCanRebin);
  track_reco_hit_no_.SetDirectory(0);
  RegisterHistogram(track_reco_hit_no_); 
   
  sprintf(name, "track_chi_sq_ndf_%04i", rp_id_);
  track_chi_sq_ndf_ = TH1F(name, name, 500, 0, 20);
  //track_chi_sq_ndf_.SetBit(TH1::kCanRebin);
  track_chi_sq_ndf_.SetDirectory(0);
  RegisterHistogram(track_chi_sq_ndf_);
}

void RPTrackInfo::Finalize()
{
  log10_t_dist_acceptance_.Divide(&log10_t_dist_seen_tracks_, &log10_t_dist_);
  log_ksi_log_t_acceptance_.Divide(&log_ksi_log_t_seen_tracks_, &log_ksi_log_t_distribution_);
}


void RPTrackInfo::Fill_t_Dist(double t)
{
  log10_t_dist_.Fill(TMath::Log10(-t));
}


void RPTrackInfo::FillTrackSeen(double t)
{
  log10_t_dist_seen_tracks_.Fill(TMath::Log10(-t));
}


void RPTrackInfo::Fill_ksi_t_Dist(double t, double ksi)
{
  log_ksi_log_t_distribution_.Fill(TMath::Log10(-t), TMath::Log10(-ksi));
}


void RPTrackInfo::FillTrackSeen(double t, double ksi)
{
  log_ksi_log_t_seen_tracks_.Fill(TMath::Log10(-t), TMath::Log10(-ksi));
}


void RPTrackInfo::FillRPPositionResiduals(double delta_x, double delta_y)
{
  position_residuum_x_.Fill(delta_x);
  position_residuum_y_.Fill(delta_y);
}


void RPTrackInfo::FillRPPositionResiduals2ProtonsSeenAt150Station(double delta_x, double delta_y)
{
  position_residuum_x_if_exactly_2_rps_in_150_station_.Fill(delta_x);
  position_residuum_y_if_exactly_2_rps_in_150_station_.Fill(delta_y);
}


void RPTrackInfo::FillPools(double delta_x, double delta_y)
{
  position_pool_x_.Fill(delta_x);
  position_pool_y_.Fill(delta_y);
}

void RPTrackInfo::FillBasicTrackInfo(const RPFittedTrack& track)
{
  track_distribution_xy_.Fill(track.X0(), track.Y0());
  track_distribution_x_.Fill(track.X0());
  track_distribution_y_.Fill(track.Y0());
  track_distribution_thx_.Fill(TMath::ATan(track.GetTx()));
  track_distribution_thy_.Fill(TMath::ATan(track.GetTy()));
  track_reco_hit_no_.Fill(track.GetHitEntries());
  track_chi_sq_ndf_.Fill(track.ChiSquaredOverN());
}


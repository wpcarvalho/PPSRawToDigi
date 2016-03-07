#include "TotemRPValidation/RPReconstructedTracksValidation/interface/RPTrackInfo.h"

RPTrackInfo::RPTrackInfo(const std::string &path, RPId rp_id, const edm::ParameterSet& conf)
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
  log10_t_dist_.SetXTitle("Log(-t/GeV^{2})");
  log10_t_dist_.SetYTitle("Entries");
  
  sprintf(name, "log10_t_dist_seen_tracks_%04i", rp_id_);
  log10_t_dist_seen_tracks_ = TH1F(name, name, 840, -4, 10);
  log10_t_dist_seen_tracks_.SetDirectory(0);
  RegisterHistogram(log10_t_dist_seen_tracks_);
  log10_t_dist_seen_tracks_.SetXTitle("Log(-t/GeV^{2})");
  log10_t_dist_seen_tracks_.SetYTitle("Entries");
  
  sprintf(name, "log10_t_dist_acceptance_%04i", rp_id_);
  log10_t_dist_acceptance_ = TH1F(name, name, 840, -4, 10);
  log10_t_dist_acceptance_.SetDirectory(0);
  RegisterHistogram(log10_t_dist_acceptance_);
  log10_t_dist_acceptance_.SetXTitle("Log(-t/GeV^{2})");
  log10_t_dist_acceptance_.SetYTitle("Acceptance");
  
  sprintf(name, "log10_xi_dist_%04i", rp_id_);
  log10_xi_dist_ = TH1F(name, name, 800, -8, 0);
  log10_xi_dist_.SetDirectory(0);
  RegisterHistogram(log10_xi_dist_);
  log10_xi_dist_.SetXTitle("Log(-#xi)");
  log10_xi_dist_.SetYTitle("Entries");
  
  sprintf(name, "log10_xi_dist_seen_tracks_%04i", rp_id_);
  log10_xi_dist_seen_tracks_ = TH1F(name, name, 800, -8, 0);
  log10_xi_dist_seen_tracks_.SetDirectory(0);
  RegisterHistogram(log10_xi_dist_seen_tracks_);
  log10_xi_dist_seen_tracks_.SetXTitle("Log(-#xi)");
  log10_xi_dist_seen_tracks_.SetYTitle("Entries");
  
  sprintf(name, "log10_xi_dist_acceptance_%04i", rp_id_);
  log10_xi_dist_acceptance_ = TH1F(name, name, 800, -8, 0);
  log10_xi_dist_acceptance_.SetDirectory(0);
  RegisterHistogram(log10_xi_dist_acceptance_);
  log10_xi_dist_acceptance_.SetXTitle("Log(-#xi)");
  log10_xi_dist_acceptance_.SetYTitle("Acceptance");
  
  sprintf(name, "log_ksi_log_t_distribution_%04i", rp_id_);
  log_ksi_log_t_distribution_ = TH2F(name, name, 50, -4, 1, 70, -7, 0.);
  log_ksi_log_t_distribution_.SetDirectory(0);
  RegisterHistogram(log_ksi_log_t_distribution_);
  log_ksi_log_t_distribution_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_distribution_.SetYTitle("Log(-#xi)");
  
  sprintf(name, "log_ksi_log_t_seen_tracks_%04i", rp_id_);
  log_ksi_log_t_seen_tracks_ = TH2F(name, name, 50, -4, 1, 70, -7, 0.);
  log_ksi_log_t_seen_tracks_.SetDirectory(0);
  RegisterHistogram(log_ksi_log_t_seen_tracks_);
  log_ksi_log_t_seen_tracks_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_seen_tracks_.SetYTitle("Log(-#xi)");
  
  sprintf(name, "log_ksi_log_t_acceptance_%04i", rp_id_);
  log_ksi_log_t_acceptance_ = TH2F(name, name, 50, -4, 1, 70, -7, 0.);
  log_ksi_log_t_acceptance_.SetDirectory(0);
  RegisterHistogram(log_ksi_log_t_acceptance_);
  log_ksi_log_t_acceptance_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_acceptance_.SetYTitle("Log(-#xi)");
  
  sprintf(name, "position_residuum_x_%04i", rp_id_);
  position_residuum_x_ = TH1F(name, name, 2000, -1, 1);
  position_residuum_x_.SetDirectory(0);
  RegisterHistogram(position_residuum_x_);
  position_residuum_x_.SetXTitle("x [mm]");
  position_residuum_x_.SetYTitle("Entries");
  
  sprintf(name, "position_residuum_y_%04i", rp_id_);
  position_residuum_y_ = TH1F(name, name, 2000, -1, 1);
  position_residuum_y_.SetDirectory(0);
  RegisterHistogram(position_residuum_y_);
  position_residuum_y_.SetXTitle("y [mm]");
  position_residuum_y_.SetYTitle("Entries");
  
  sprintf(name, "position_residuum_x_if_exactly_2_rps_in_150_station_%04i", rp_id_);
  position_residuum_x_if_exactly_2_rps_in_150_station_ = TH1F(name, name, 2000, -1, 1);
  position_residuum_x_if_exactly_2_rps_in_150_station_.SetDirectory(0);
  RegisterHistogram(position_residuum_x_if_exactly_2_rps_in_150_station_);
  position_residuum_x_if_exactly_2_rps_in_150_station_.SetXTitle("x [mm]");
  position_residuum_x_if_exactly_2_rps_in_150_station_.SetYTitle("Entries");
  
  sprintf(name, "position_residuum_y_if_exactly_2_rps_in_150_station_%04i", rp_id_);
  position_residuum_y_if_exactly_2_rps_in_150_station_ = TH1F(name, name, 2000, -1, 1);
  position_residuum_y_if_exactly_2_rps_in_150_station_.SetDirectory(0);
  RegisterHistogram(position_residuum_y_if_exactly_2_rps_in_150_station_);
  position_residuum_y_if_exactly_2_rps_in_150_station_.SetXTitle("y [mm]");
  position_residuum_y_if_exactly_2_rps_in_150_station_.SetYTitle("Entries");
  
  sprintf(name, "position_pool_x_%04i", rp_id_);
  position_pool_x_ = TH1F(name, name, 4000, -30, 30);
  position_pool_x_.SetDirectory(0);
  RegisterHistogram(position_pool_x_); 
  position_pool_x_.SetXTitle("Normalised position x");
  position_pool_x_.SetYTitle("Entries");
   
  sprintf(name, "position_pool_y_%04i", rp_id_);
  position_pool_y_ = TH1F(name, name, 4000, -30, 30);
  position_pool_y_.SetDirectory(0);
  RegisterHistogram(position_pool_y_); 
  position_pool_y_.SetXTitle("Normalised position y");
  position_pool_y_.SetYTitle("Entries");
   
  sprintf(name, "track_distribution_xy_%04i", rp_id_);
  track_distribution_xy_ = TH2F(name, name, 80, -40, 40, 80, -40, 40);
  track_distribution_xy_.SetDirectory(0);
  RegisterHistogram(track_distribution_xy_);
  track_distribution_xy_.SetXTitle("x [mm]");
  track_distribution_xy_.SetYTitle("y [mm]");
   
  sprintf(name, "track_distribution_x_%04i", rp_id_);
  track_distribution_x_ = TH1F(name, name, 400, -40, 40);
  track_distribution_x_.SetDirectory(0);
  RegisterHistogram(track_distribution_x_); 
  track_distribution_x_.SetXTitle("x [mm]");
  track_distribution_x_.SetYTitle("Entries");
   
  sprintf(name, "track_distribution_y_%04i", rp_id_);
  track_distribution_y_ = TH1F(name, name, 400, -40, 40);
  track_distribution_y_.SetDirectory(0);
  RegisterHistogram(track_distribution_y_);
  track_distribution_y_.SetXTitle("y [mm]");
  track_distribution_y_.SetYTitle("Entries");
   
  sprintf(name, "track_distribution_thx_%04i", rp_id_);
  track_distribution_thx_ = TH1F(name, name, 400, -0.0001, 0.0001);
  track_distribution_thx_.SetDirectory(0);
  RegisterHistogram(track_distribution_thx_); 
  track_distribution_thx_.SetXTitle("#Theta_{x} [rad]");
  track_distribution_thx_.SetYTitle("Entries");
   
  sprintf(name, "track_distribution_thy_%04i", rp_id_);
  track_distribution_thy_ = TH1F(name, name, 400, -0.0001, 0.0001);
  track_distribution_thy_.SetDirectory(0);
  RegisterHistogram(track_distribution_thy_); 
  track_distribution_thy_.SetXTitle("#Theta_{y} [rad]");
  track_distribution_thy_.SetYTitle("Entries");
   
  sprintf(name, "track_reco_hit_no_%04i", rp_id_);
  track_reco_hit_no_ = TH1F(name, name, 11, -0.5, 10.5);
  track_reco_hit_no_.SetDirectory(0);
  RegisterHistogram(track_reco_hit_no_); 
  track_reco_hit_no_.SetXTitle("Spatial points per track");
  track_reco_hit_no_.SetYTitle("Entries");
   
  sprintf(name, "track_chi_sq_ndf_%04i", rp_id_);
  track_chi_sq_ndf_ = TH1F(name, name, 500, 0, 20);
  track_chi_sq_ndf_.SetDirectory(0);
  RegisterHistogram(track_chi_sq_ndf_);
  track_chi_sq_ndf_.SetXTitle("#Chi^{2}/NDF");
  track_chi_sq_ndf_.SetYTitle("Entries");
  
  sprintf(name, "u_sectors_on_%04i", rp_id_);
  u_sectors_on_ = TH1F(name, name, 16, -0.5, 15.5);
  u_sectors_on_.SetDirectory(0);
  RegisterHistogram(u_sectors_on_);
  u_sectors_on_.SetXTitle("u coincidence sectors");
  u_sectors_on_.SetYTitle("Entries");

  sprintf(name, "v_sectors_on_%04i", rp_id_);
  v_sectors_on_ = TH1F(name, name, 16, -0.5, 15.5);
  v_sectors_on_.SetDirectory(0);
  RegisterHistogram(v_sectors_on_);
  v_sectors_on_.SetXTitle("v coincidence sectors");
  v_sectors_on_.SetYTitle("Entries");
   
  sprintf(name, "uv_sectors_on_%04i", rp_id_);
  uv_sectors_on_ = TH2F(name, name, 16, -0.5, 15.5, 16, -0.5, 15.5);
  uv_sectors_on_.SetDirectory(0);
  RegisterHistogram(uv_sectors_on_);
  uv_sectors_on_.SetXTitle("u coincidence sectors");
  uv_sectors_on_.SetYTitle("v coincidence sectors");
}

void RPTrackInfo::Finalize()
{
  log10_t_dist_acceptance_.Divide(&log10_t_dist_seen_tracks_, &log10_t_dist_);
  log_ksi_log_t_acceptance_.Divide(&log_ksi_log_t_seen_tracks_, &log_ksi_log_t_distribution_);
  log10_xi_dist_acceptance_.Divide(&log10_xi_dist_seen_tracks_, &log10_xi_dist_);
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
  log10_xi_dist_.Fill(TMath::Log10(-ksi));
}


void RPTrackInfo::FillTrackSeen(double t, double ksi)
{
  log_ksi_log_t_seen_tracks_.Fill(TMath::Log10(-t), TMath::Log10(-ksi));
  log10_xi_dist_seen_tracks_.Fill(TMath::Log10(-ksi));
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


void RPTrackInfo::FillRPCoincidenceChipInfo(const RPCCBits& u_coincidence_bits, const RPCCBits& v_coincidence_bits)
{
  const bool* u_bits = u_coincidence_bits.getRawBS();
  const bool* v_bits = v_coincidence_bits.getRawBS();
  
//  std::cout<<"RPTrackInfo::FillRPCoincidenceChipInfo RP:"<<rp_id_;
  
  std::vector<int> u_sectors, v_sectors;
  
  for(int ui=0; ui<16; ++ui)
  {
    if(u_bits[ui])
    {
      u_sectors_on_.Fill(ui);
      u_sectors.push_back(ui);
    }
  }
  
  for(int vi=0; vi<16; ++vi)
  {
    if(v_bits[vi])
    {
      v_sectors_on_.Fill(vi);
      v_sectors.push_back(vi);
    }
  }
  
  for(unsigned int ui=0; ui<u_sectors.size(); ++ui)
  {
    for(unsigned int vi=0; vi<v_sectors.size(); ++vi)
    {
      uv_sectors_on_.Fill(u_sectors[ui], v_sectors[vi]);
    }
  }
}


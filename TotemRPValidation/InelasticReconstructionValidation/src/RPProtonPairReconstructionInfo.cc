#include "TotemRPValidation/InelasticReconstructionValidation/interface/RPProtonPairReconstructionInfo.h"
#include "TProfile3D.h"
#include "TMath.h"
#include "TH2F.h"
#include "TotemRPValidation/ValidationTools/interface/ReconstructionProfile.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"



RPProtonPairReconstructionInfo::RPProtonPairReconstructionInfo(
    const std::string &path, const edm::ParameterSet& conf)
: BaseHistogramManager(path)
{
  chi2_upper_tail_prob_cut_ = conf.getParameter<double>("Chi2UpperTailProbCut");
  InitializeHistograms();
}


void RPProtonPairReconstructionInfo::InitializeHistograms()
{  
  char name[1024];

  sprintf(name, "log_ksi_log_t_distribution_left_");
  log_ksi_log_t_distribution_left_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_distribution_left_.SetDirectory(0);
  log_ksi_log_t_distribution_left_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_distribution_left_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_distribution_left_);
  
  sprintf(name, "log_ksi_log_t_distribution_right_");
  log_ksi_log_t_distribution_right_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_distribution_right_.SetDirectory(0);
  log_ksi_log_t_distribution_right_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_distribution_right_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_distribution_right_);
  
  sprintf(name, "log_ksi_log_t_reconstructed_protons_left_");
  log_ksi_log_t_reconstructed_protons_left_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_reconstructed_protons_left_.SetDirectory(0);
  log_ksi_log_t_reconstructed_protons_left_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_reconstructed_protons_left_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_reconstructed_protons_left_);
  
  sprintf(name, "log_ksi_log_t_reconstructed_protons_right_");
  log_ksi_log_t_reconstructed_protons_right_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_reconstructed_protons_right_.SetDirectory(0);
  log_ksi_log_t_reconstructed_protons_right_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_reconstructed_protons_right_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_reconstructed_protons_right_);
  
  sprintf(name, "log_ksi_log_t_reconstruction_failed_left_");
  log_ksi_log_t_reconstruction_failed_left_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_reconstruction_failed_left_.SetDirectory(0);
  log_ksi_log_t_reconstruction_failed_left_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_reconstruction_failed_left_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_reconstruction_failed_left_);
  
  sprintf(name, "log_ksi_log_t_reconstruction_failed_right_");
  log_ksi_log_t_reconstruction_failed_right_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_reconstruction_failed_right_.SetDirectory(0);
  log_ksi_log_t_reconstruction_failed_right_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_reconstruction_failed_right_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_reconstruction_failed_right_);
  
  sprintf(name, "log_ksi_log_t_reconstruction_failure_rate_left_");
  log_ksi_log_t_reconstruction_failure_rate_left_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_reconstruction_failure_rate_left_.SetDirectory(0);
  log_ksi_log_t_reconstruction_failure_rate_left_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_reconstruction_failure_rate_left_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_reconstruction_failure_rate_left_);
  
  sprintf(name, "log_ksi_log_t_reconstruction_failure_rate_right_");
  log_ksi_log_t_reconstruction_failure_rate_right_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_reconstruction_failure_rate_right_.SetDirectory(0);
  log_ksi_log_t_reconstruction_failure_rate_right_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_reconstruction_failure_rate_right_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_reconstruction_failure_rate_right_);
  
  sprintf(name, "log_ksi_log_t_total_accepted_left_");
  log_ksi_log_t_total_accepted_left_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_total_accepted_left_.SetDirectory(0);
  log_ksi_log_t_total_accepted_left_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_total_accepted_left_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_total_accepted_left_);
  
  sprintf(name, "log_ksi_log_t_total_accepted_right_");
  log_ksi_log_t_total_accepted_right_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_total_accepted_right_.SetDirectory(0);
  log_ksi_log_t_total_accepted_right_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_total_accepted_right_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_total_accepted_right_);
  
  sprintf(name, "log_ksi_log_t_reconstruction_failure_to_total_accepted_left_");
  log_ksi_log_t_reconstruction_failure_to_total_accepted_left_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_reconstruction_failure_to_total_accepted_left_.SetDirectory(0);
  log_ksi_log_t_reconstruction_failure_to_total_accepted_left_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_reconstruction_failure_to_total_accepted_left_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_reconstruction_failure_to_total_accepted_left_);
  
  sprintf(name, "log_ksi_log_t_reconstruction_failure_to_total_accepted_right_");
  log_ksi_log_t_reconstruction_failure_to_total_accepted_right_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_reconstruction_failure_to_total_accepted_right_.SetDirectory(0);
  log_ksi_log_t_reconstruction_failure_to_total_accepted_right_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_reconstruction_failure_to_total_accepted_right_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_reconstruction_failure_to_total_accepted_right_);
  
  sprintf(name, "log_ksi_log_t_acceptance_left_");
  log_ksi_log_t_acceptance_left_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_acceptance_left_.SetDirectory(0);
  log_ksi_log_t_acceptance_left_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_acceptance_left_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_acceptance_left_);
  
  sprintf(name, "log_ksi_log_t_acceptance_right_");
  log_ksi_log_t_acceptance_right_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_acceptance_right_.SetDirectory(0);
  log_ksi_log_t_acceptance_right_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_acceptance_right_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_acceptance_right_);
  
  sprintf(name, "log_ksi_log_t_total_acceptance_left_");
  log_ksi_log_t_total_acceptance_left_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_total_acceptance_left_.SetDirectory(0);
  log_ksi_log_t_total_acceptance_left_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_total_acceptance_left_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_total_acceptance_left_);
  
  sprintf(name, "log_ksi_log_t_total_acceptance_right_");
  log_ksi_log_t_total_acceptance_right_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_total_acceptance_right_.SetDirectory(0);
  log_ksi_log_t_total_acceptance_right_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_total_acceptance_right_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_total_acceptance_right_);
  
  //elastic histograms
  sprintf(name, "elast_log_t_distribution_");
  elast_log_t_distribution_ = TH1F(name, name, 800, -6, 2);
  elast_log_t_distribution_.SetDirectory(0);
  elast_log_t_distribution_.SetXTitle("Log(-t/GeV^{2})");
  elast_log_t_distribution_.SetYTitle("Entries");
  RegisterHistogram(elast_log_t_distribution_);

  sprintf(name, "elast_log_t_accepted_");
  elast_log_t_accepted_ = TH1F(name, name, 800, -6, 2);
  elast_log_t_accepted_.SetDirectory(0);
  elast_log_t_accepted_.SetXTitle("Log(-t/GeV^{2})");
  elast_log_t_accepted_.SetYTitle("Reconstructed protons");
  RegisterHistogram(elast_log_t_accepted_);

  sprintf(name, "elast_log_t_acceptance_");
  elast_log_t_acceptance_ = TH1F(name, name, 800, -6, 2);
  elast_log_t_acceptance_.SetDirectory(0);
  elast_log_t_acceptance_.SetXTitle("Log(-t/GeV^{2})");
  elast_log_t_acceptance_.SetYTitle("Acceptance");
  RegisterHistogram(elast_log_t_acceptance_);

  sprintf(name, "elast_log_t_total_accepted_");
  elast_log_t_total_accepted_ = TH1F(name, name, 800, -6, 2);
  elast_log_t_total_accepted_.SetDirectory(0);
  elast_log_t_total_accepted_.SetXTitle("Log(-t/GeV^{2})");
  elast_log_t_total_accepted_.SetYTitle("Accepted protons");
  RegisterHistogram(elast_log_t_total_accepted_);

  sprintf(name, "elast_log_t_total_acceptance_");
  elast_log_t_total_acceptance_ = TH1F(name, name, 800, -6, 2);
  elast_log_t_total_acceptance_.SetDirectory(0);
  elast_log_t_total_acceptance_.SetXTitle("Log(-t/GeV^{2})");
  elast_log_t_total_acceptance_.SetYTitle("Total acceptance");
  RegisterHistogram(elast_log_t_total_acceptance_);
  
  sprintf(name, "elast_left_thx_rec_error_");
  elast_left_thx_rec_error_ = TH1F(name, name, 500, -0.0001, 0.0001);
  elast_left_thx_rec_error_.SetBit(TH1::kCanRebin);
  elast_left_thx_rec_error_.SetDirectory(0);
  elast_left_thx_rec_error_.SetXTitle("Right proton #Delta#Theta_{x} [rad]");
  elast_left_thx_rec_error_.SetYTitle("Entries");
  RegisterHistogram(elast_left_thx_rec_error_);
  
  sprintf(name, "elast_left_thy_rec_error_");
  elast_left_thy_rec_error_ = TH1F(name, name, 500, -0.0001, 0.0001);
  elast_left_thy_rec_error_.SetBit(TH1::kCanRebin);
  elast_left_thy_rec_error_.SetDirectory(0);
  elast_left_thy_rec_error_.SetXTitle("Right proton #Delta#Theta_{y} [rad]");
  elast_left_thy_rec_error_.SetYTitle("Entries");
  RegisterHistogram(elast_left_thy_rec_error_);
  
  sprintf(name, "elast_right_thx_rec_error_");
  elast_right_thx_rec_error_ = TH1F(name, name, 500, -0.0001, 0.0001);
  elast_right_thx_rec_error_.SetBit(TH1::kCanRebin);
  elast_right_thx_rec_error_.SetDirectory(0);
  elast_right_thx_rec_error_.SetXTitle("Left proton #Delta#Theta_{x} [rad]");
  elast_right_thx_rec_error_.SetYTitle("Entries");
  RegisterHistogram(elast_right_thx_rec_error_);
  
  sprintf(name, "elast_right_thy_rec_error_");
  elast_right_thy_rec_error_ = TH1F(name, name, 500, -0.0001, 0.0001);
  elast_right_thy_rec_error_.SetBit(TH1::kCanRebin);
  elast_right_thy_rec_error_.SetDirectory(0);
  elast_right_thy_rec_error_.SetXTitle("Left proton #Delta#Theta_{y} [rad]");
  elast_right_thy_rec_error_.SetYTitle("Entries");
  RegisterHistogram(elast_right_thy_rec_error_);
  
  sprintf(name, "elast_scattering_thx_rec_error_");
  elast_scattering_thx_rec_error_ = TH1F(name, name, 500, -0.0001, 0.0001);
  elast_scattering_thx_rec_error_.SetBit(TH1::kCanRebin);
  elast_scattering_thx_rec_error_.SetDirectory(0);
  elast_scattering_thx_rec_error_.SetXTitle("#Delta#Theta_{x} [rad]");
  elast_scattering_thx_rec_error_.SetYTitle("Entries");
  RegisterHistogram(elast_scattering_thx_rec_error_);
  
  sprintf(name, "elast_scattering_thy_rec_error_");
  elast_scattering_thy_rec_error_ = TH1F(name, name, 500, -0.0001, 0.0001);
  elast_scattering_thy_rec_error_.SetBit(TH1::kCanRebin);
  elast_scattering_thy_rec_error_.SetDirectory(0);
  elast_scattering_thy_rec_error_.SetXTitle("#Delta#Theta_{y} [rad]");
  elast_scattering_thy_rec_error_.SetYTitle("Entries");
  RegisterHistogram(elast_scattering_thy_rec_error_);
  
  sprintf(name, "elast_left_xi_error_");
  elast_left_xi_error_ = TH1F(name, name, 500, -0.001, 0.001);
  elast_left_xi_error_.SetBit(TH1::kCanRebin);
  elast_left_xi_error_.SetDirectory(0);
  elast_left_xi_error_.SetXTitle("Left proton #Delta#xi");
  elast_left_xi_error_.SetYTitle("Entries");
  RegisterHistogram(elast_left_xi_error_);

  sprintf(name, "elast_right_xi_error_");
  elast_right_xi_error_ = TH1F(name, name, 500, -0.001, 0.001);
  elast_right_xi_error_.SetBit(TH1::kCanRebin);
  elast_right_xi_error_.SetDirectory(0);
  elast_right_xi_error_.SetXTitle("Right proton #Delta#xi");
  elast_right_xi_error_.SetYTitle("Entries");
  RegisterHistogram(elast_right_xi_error_);

  sprintf(name, "elast_thx_left_right_difference_");
  elast_thx_left_right_difference_ = TH1F(name, name, 500, -0.0001, 0.0010);
  elast_thx_left_right_difference_.SetBit(TH1::kCanRebin);
  elast_thx_left_right_difference_.SetDirectory(0);
  elast_thx_left_right_difference_.SetXTitle("Rec. #Theta_{x,L}-#Theta_{x,R} [rad]");
  elast_thx_left_right_difference_.SetYTitle("Entries");
  RegisterHistogram(elast_thx_left_right_difference_);

  sprintf(name, "elast_thy_left_right_difference_");
  elast_thy_left_right_difference_ = TH1F(name, name, 500, -0.0001, 0.0001);
  elast_thy_left_right_difference_.SetBit(TH1::kCanRebin);
  elast_thy_left_right_difference_.SetDirectory(0);
  elast_thy_left_right_difference_.SetXTitle("Rec. #Theta_{y,L}-#Theta_{y,R} [rad]");
  elast_thy_left_right_difference_.SetYTitle("Entries");
  RegisterHistogram(elast_thy_left_right_difference_);

  sprintf(name, "elast_x0_rec_error_");
  elast_x0_rec_error_ = TH1F(name, name, 500, -0.1, 0.1);
  elast_x0_rec_error_.SetBit(TH1::kCanRebin);
  elast_x0_rec_error_.SetDirectory(0);
  elast_x0_rec_error_.SetXTitle("#Delta X_{0} [mm]");
  elast_x0_rec_error_.SetYTitle("Entries");
  RegisterHistogram(elast_x0_rec_error_);

  sprintf(name, "elast_y0_rec_error_");
  elast_y0_rec_error_ = TH1F(name, name, 500, -0.1, 0.1);
  elast_y0_rec_error_.SetBit(TH1::kCanRebin);
  elast_y0_rec_error_.SetDirectory(0);
  elast_y0_rec_error_.SetXTitle("#Delta Y_{0} [mm]");
  elast_y0_rec_error_.SetYTitle("Entries");
  RegisterHistogram(elast_y0_rec_error_);

  sprintf(name, "elast_z0_rec_error_");
  elast_z0_rec_error_ = TH1F(name, name, 500, -0.1, 0.1);
  elast_z0_rec_error_.SetBit(TH1::kCanRebin);
  elast_z0_rec_error_.SetDirectory(0);
  elast_z0_rec_error_.SetXTitle("#Delta Z_{0} [mm]");
  elast_z0_rec_error_.SetYTitle("Entries");
  RegisterHistogram(elast_z0_rec_error_);

  sprintf(name, "elast_chi2ndf_dist_");
  elast_chi2ndf_dist_ = TH1F(name, name, 500, 0, 10);
  elast_chi2ndf_dist_.SetBit(TH1::kCanRebin);
  elast_chi2ndf_dist_.SetDirectory(0);
  elast_chi2ndf_dist_.SetXTitle("#chi^{2}/NDF");
  elast_chi2ndf_dist_.SetYTitle("Entries");
  RegisterHistogram(elast_chi2ndf_dist_);

  sprintf(name, "elast_probe_dist_");
  elast_probe_dist_ = TH1F(name, name, 500, 0, 1);
  elast_probe_dist_.SetBit(TH1::kCanRebin);
  elast_probe_dist_.SetDirectory(0);
  elast_probe_dist_.SetXTitle("#chi^{2} upper tail dist.");
  elast_probe_dist_.SetYTitle("Entries");
  RegisterHistogram(elast_probe_dist_);
  
  
  //left side
  sprintf(name, "ksi_error_vs_log_t_log_ksi_phi_left_");
  ksi_error_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s"); 
  RegisterHistogram(ksi_error_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "smeared_ksi_error_vs_log_t_log_ksi_phi_left_");
  smeared_ksi_error_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s"); 
  RegisterHistogram(smeared_ksi_error_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "t_error_vs_log_t_log_ksi_phi_left_");
  t_error_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //t_error_vs_log_t_log_ksi_phi_left_.SetDirectory(0);
  RegisterHistogram(t_error_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "phi_error_vs_log_t_log_ksi_phi_left_");
  phi_error_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(phi_error_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "vx_error_vs_log_t_log_ksi_phi_left_");
  vx_error_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //vx_error_vs_log_t_log_ksi_phi_left_.SetDirectory(0);
  RegisterHistogram(vx_error_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "vy_error_vs_log_t_log_ksi_phi_left_");
  vy_error_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //vy_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(vy_error_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "vz_error_vs_log_t_log_ksi_phi_left_");
  vz_error_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //vy_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(vz_error_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "thetax_error_vs_log_t_log_ksi_phi_left_");
  thetax_error_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //thetax_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(thetax_error_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "thetay_error_vs_log_t_log_ksi_phi_left_");
  thetay_error_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //thetay_error_vs_log_t_log_ksi_phi_left_.SetDirectory(0);
  RegisterHistogram(thetay_error_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "smeared_thetax_error_vs_log_t_log_ksi_phi_left_");
  smeared_thetax_error_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //smeared_thetax_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(smeared_thetax_error_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "smeared_thetay_error_vs_log_t_log_ksi_phi_left_");
  smeared_thetay_error_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //smeared_thetay_error_vs_log_t_log_ksi_phi_left_.SetDirectory(0);
  RegisterHistogram(smeared_thetay_error_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "chi2_error_overn_N_vs_log_t_log_ksi_phi_left_");
  chi2_error_overn_N_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //chi2_error_overn_N_vs_log_t_log_ksi_phi_left_.SetDirectory(0);
  RegisterHistogram(chi2_error_overn_N_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "ksi_rel_error_vs_log_t_log_ksi_phi_left_");
  ksi_rel_error_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(ksi_rel_error_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "t_rel_error_vs_log_t_log_ksi_phi_left_");
  t_rel_error_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_rel_error_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "reconstruction_failure_vs_log_t_log_ksi_phi_left_");
  reconstruction_failure_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(reconstruction_failure_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "prob_function_vs_log_t_log_ksi_phi_left_");
  prob_function_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(prob_function_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "t_x_error_vs_log_t_log_ksi_phi_left_");
  t_x_error_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_x_error_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "t_y_error_vs_log_t_log_ksi_phi_left_");
  t_y_error_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_y_error_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "t_x_rel_error_vs_log_t_log_ksi_phi_left_");
  t_x_rel_error_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_x_rel_error_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "t_y_rel_error_vs_log_t_log_ksi_phi_left_");
  t_y_rel_error_vs_log_t_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_y_rel_error_vs_log_t_log_ksi_phi_left_);
  
  sprintf(name, "t_x_error_vs_log_tx_log_ksi_phi_left_");
  t_x_error_vs_log_tx_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_x_error_vs_log_tx_log_ksi_phi_left_);
  
  sprintf(name, "t_y_error_vs_log_ty_log_ksi_phi_left_");
  t_y_error_vs_log_ty_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_y_error_vs_log_ty_log_ksi_phi_left_);
  
  sprintf(name, "t_x_rel_error_vs_log_tx_log_ksi_phi_left_");
  t_x_rel_error_vs_log_tx_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_x_rel_error_vs_log_tx_log_ksi_phi_left_);
  
  sprintf(name, "t_y_rel_error_vs_log_ty_log_ksi_phi_left_");
  t_y_rel_error_vs_log_ty_log_ksi_phi_left_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_y_rel_error_vs_log_ty_log_ksi_phi_left_);
  
  
  
  // right side
  sprintf(name, "ksi_error_vs_log_t_log_ksi_phi_right_");
  ksi_error_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s"); 
  //ksi_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(ksi_error_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "smeared_ksi_error_vs_log_t_log_ksi_phi_right_");
  smeared_ksi_error_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s"); 
  //smeared_ksi_error_vs_log_t_log_ksi_phi_right_.SetDirectory(0);
  RegisterHistogram(smeared_ksi_error_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "t_error_vs_log_t_log_ksi_phi_right_");
  t_error_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //t_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(t_error_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "phi_error_vs_log_t_log_ksi_phi_right_");
  phi_error_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //phi_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(phi_error_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "vx_error_vs_log_t_log_ksi_phi_right_");
  vx_error_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //vx_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(vx_error_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "vy_error_vs_log_t_log_ksi_phi_right_");
  vy_error_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //vy_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(vy_error_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "vz_error_vs_log_t_log_ksi_phi_right_");
  vz_error_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //vy_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(vz_error_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "thetax_error_vs_log_t_log_ksi_phi_right_");
  thetax_error_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //thetax_error_vs_log_t_log_ksi_phi_right_.SetDirectory(0);
  RegisterHistogram(thetax_error_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "thetay_error_vs_log_t_log_ksi_phi_right_");
  thetay_error_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //thetay_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(thetay_error_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "smeared_thetax_error_vs_log_t_log_ksi_phi_right_");
  smeared_thetax_error_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //smeared_thetax_error_vs_log_t_log_ksi_phi_right_.SetDirectory(0);
  RegisterHistogram(smeared_thetax_error_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "smeared_thetay_error_vs_log_t_log_ksi_phi_right_");
  smeared_thetay_error_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //smeared_thetay_error_vs_log_t_log_ksi_phi_right_.SetDirectory(0);
  RegisterHistogram(smeared_thetay_error_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "chi2_error_overn_N_vs_log_t_log_ksi_phi_right_");
  chi2_error_overn_N_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //chi2_error_overn_N_vs_log_t_log_ksi_phi_right_.SetDirectory(0);
  RegisterHistogram(chi2_error_overn_N_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "ksi_rel_error_vs_log_t_log_ksi_phi_right_");
  ksi_rel_error_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(ksi_rel_error_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "t_rel_error_vs_log_t_log_ksi_phi_right_");
  t_rel_error_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_rel_error_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "reconstruction_failure_vs_log_t_log_ksi_phi_right_");
  reconstruction_failure_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(reconstruction_failure_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "prob_function_vs_log_t_log_ksi_phi_right_");
  prob_function_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(prob_function_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "t_x_error_vs_log_t_log_ksi_phi_right_");
  t_x_error_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_x_error_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "t_y_error_vs_log_t_log_ksi_phi_right_");
  t_y_error_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_y_error_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "t_x_rel_error_vs_log_t_log_ksi_phi_right_");
  t_x_rel_error_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_x_rel_error_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "t_y_rel_error_vs_log_t_log_ksi_phi_right_");
  t_y_rel_error_vs_log_t_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_y_rel_error_vs_log_t_log_ksi_phi_right_);
  
  sprintf(name, "t_x_error_vs_log_tx_log_ksi_phi_right_");
  t_x_error_vs_log_tx_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_x_error_vs_log_tx_log_ksi_phi_right_);
  
  sprintf(name, "t_y_error_vs_log_ty_log_ksi_phi_right_");
  t_y_error_vs_log_ty_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_y_error_vs_log_ty_log_ksi_phi_right_);
  
  sprintf(name, "t_x_rel_error_vs_log_tx_log_ksi_phi_right_");
  t_x_rel_error_vs_log_tx_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_x_rel_error_vs_log_tx_log_ksi_phi_right_);
  
  sprintf(name, "t_y_rel_error_vs_log_ty_log_ksi_phi_right_");
  t_y_rel_error_vs_log_ty_log_ksi_phi_right_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_y_rel_error_vs_log_ty_log_ksi_phi_right_);
  
  
  
  //common histograms
  sprintf(name, "M_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  M_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(M_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "M_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  M_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(M_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "ksi_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  ksi_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(ksi_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "ksi_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  ksi_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(ksi_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "ksi_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  ksi_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(ksi_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "ksi_higher_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  ksi_higher_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(ksi_higher_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "ksi_lower_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  ksi_lower_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(ksi_lower_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "ksi_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  ksi_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(ksi_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "t_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  t_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(t_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "t_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  t_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(t_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "t_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  t_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(t_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "t_higher_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  t_higher_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(t_higher_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "t_lower_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  t_lower_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(t_lower_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "t_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  t_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(t_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "phi_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  phi_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(phi_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "phi_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  phi_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(phi_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "phi_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  phi_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(phi_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "x_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  x_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(x_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "y_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  y_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(y_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "z_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  z_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(z_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "probe_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  probe_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(probe_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);

  sprintf(name, "chsqndf_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  chsqndf_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_ = 
      ReconstructionProfile(name, name, 50, 0, 5, 60, -5, 1, 70, -6, 1, "s");
  RegisterHistogram(chsqndf_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_);
  
  log_dpe_mass_total_accepted_dist_ = TH1F("log_dpe_mass_total_accepted_dist_","log_dpe_mass_total_accepted_dist_",100,0,4);
  RegisterHistogram(log_dpe_mass_total_accepted_dist_);
  log_dpe_mass_reconstructed_dist_ = TH1F("log_dpe_mass_reconstructed_dist_","log_dpe_mass_reconstructed_dist_",100,0,4);
  RegisterHistogram(log_dpe_mass_reconstructed_dist_);
  log_dpe_mass_generated_dist_ = TH1F("log_dpe_mass_generated_dist_","log_dpe_mass_generated_dist_",100,0,4);
  RegisterHistogram(log_dpe_mass_generated_dist_);
  log_dpe_mass_total_acceptance_ = TH1F("log_dpe_mass_total_acceptance_","log_dpe_mass_total_acceptance_",100,0,4);
  RegisterHistogram(log_dpe_mass_total_acceptance_);
  log_dpe_mass_reconstructed_acceptance_ = TH1F("log_dpe_mass_reconstructed_acceptance_","log_dpe_mass_reconstructed_acceptance_",100,0,4);
  RegisterHistogram(log_dpe_mass_reconstructed_acceptance_);
  log_dpe_mass_reconstruction_failures_ = TH1F("log_dpe_mass_reconstruction_failures_","log_dpe_mass_reconstruction_failures_",100,0,4);
  RegisterHistogram(log_dpe_mass_reconstruction_failures_);
  log_dpe_mass_reconstruction_failure_probability_ = TH1F("log_dpe_mass_reconstruction_failure_probability_","log_dpe_mass_reconstruction_failure_probability_",100,0,4);
  RegisterHistogram(log_dpe_mass_reconstruction_failure_probability_);
}


void RPProtonPairReconstructionInfo::Finalize()
{
  log_ksi_log_t_acceptance_left_
      .Divide(&log_ksi_log_t_reconstructed_protons_left_, &log_ksi_log_t_distribution_left_);
  log_ksi_log_t_reconstruction_failure_rate_left_
      .Divide(&log_ksi_log_t_reconstruction_failed_left_, &log_ksi_log_t_distribution_left_);
  log_ksi_log_t_reconstruction_failure_to_total_accepted_left_
      .Divide(&log_ksi_log_t_reconstruction_failed_left_, &log_ksi_log_t_total_accepted_left_);
  log_ksi_log_t_total_acceptance_left_
      .Divide(&log_ksi_log_t_total_accepted_left_, &log_ksi_log_t_distribution_left_);
      
      
  log_ksi_log_t_acceptance_right_
      .Divide(&log_ksi_log_t_reconstructed_protons_right_, &log_ksi_log_t_distribution_right_);
  log_ksi_log_t_reconstruction_failure_rate_right_
      .Divide(&log_ksi_log_t_reconstruction_failed_right_, &log_ksi_log_t_distribution_right_);
  log_ksi_log_t_reconstruction_failure_to_total_accepted_right_
      .Divide(&log_ksi_log_t_reconstruction_failed_right_, &log_ksi_log_t_total_accepted_right_);
  log_ksi_log_t_total_acceptance_right_
      .Divide(&log_ksi_log_t_total_accepted_right_, &log_ksi_log_t_distribution_right_);
  

  log_dpe_mass_total_acceptance_
      .Divide(&log_dpe_mass_total_accepted_dist_, &log_dpe_mass_generated_dist_);
  log_dpe_mass_reconstructed_acceptance_
      .Divide(&log_dpe_mass_reconstructed_dist_, &log_dpe_mass_generated_dist_);
  log_dpe_mass_reconstruction_failure_probability_
      .Divide(&log_dpe_mass_reconstruction_failures_, &log_dpe_mass_generated_dist_);
  
  elast_log_t_total_acceptance_
  .Divide(&elast_log_t_total_accepted_, &elast_log_t_distribution_);
  
  elast_log_t_acceptance_
  .Divide(&elast_log_t_accepted_, &elast_log_t_distribution_);
}


void RPProtonPairReconstructionInfo::FillReferenceHistograms(
      const BeamOpticsParams &BOPar, 
      const HepMC::GenParticle *pr_prot_left, 
      const HepMC::GenParticle *pr_prot_right, 
      const HepMC::GenParticle *pr_vert_left, const HepMC::GenParticle *pr_vert_right, 
      bool verbosity)
{
  if(!pr_prot_left || !pr_prot_right)
    return;
  
  //left
  double pr_px0 = pr_prot_left->momentum().x();
  double pr_py0 = pr_prot_left->momentum().y();
  double pr_pz0 = pr_prot_left->momentum().z();
  
  double pr_ksi0 = BOPar.IPNonSmearedProtonMomentumToXi(pr_px0, pr_py0, pr_pz0);
  double pr_t0 = BOPar.IPNonSmearedProtonMomentumTot(pr_px0, pr_py0, pr_pz0);
        
  double log_minus_ksi0 = TMath::Log10(-pr_ksi0);
  double log_minus_t0 = TMath::Log10(-pr_t0);
        
  log_ksi_log_t_distribution_left_.Fill(log_minus_t0, log_minus_ksi0);
  
  //right
  double pr_px1 = pr_prot_right->momentum().x();
  double pr_py1 = pr_prot_right->momentum().y();
  double pr_pz1 = pr_prot_right->momentum().z();
  
  double pr_ksi1 = BOPar.IPNonSmearedProtonMomentumToXi(pr_px1, pr_py1, pr_pz1);
  double pr_t1 = BOPar.IPNonSmearedProtonMomentumTot(pr_px1, pr_py1, pr_pz1);
        
  double log_minus_ksi1 = TMath::Log10(-pr_ksi1);
  double log_minus_t1 = TMath::Log10(-pr_t1);
        
  log_ksi_log_t_distribution_right_.Fill(log_minus_t1, log_minus_ksi1);
  
  double log_generated_dpe_mass = TMath::Log10(BOPar.HEPMCNonSmearedProtonstoM(pr_prot_right->momentum(), pr_prot_left->momentum()));
  log_dpe_mass_generated_dist_.Fill(log_generated_dpe_mass);
  
  //elastically scaterred protons
  double mean_t = (pr_t0+pr_t1)/2.0;
  double log_minus_mean_t = TMath::Log10(-mean_t);
  elast_log_t_distribution_.Fill(log_minus_mean_t);
}


void RPProtonPairReconstructionInfo::FillResidualHistograms(
    const BeamOpticsParams &BOPar, 
    const HepMC::GenParticle *pr_prot_left,
    const HepMC::GenParticle *pr_prot_right,
    const HepMC::GenParticle *pr_vert_left, const HepMC::GenParticle *pr_vert_right, 
    const RPReconstructedProtonPair* rec_prot_pair, bool verbosity)
{
  if(!pr_prot_left || !pr_prot_right || !rec_prot_pair)
    return;
    
  //left
  double pr_px0 = pr_prot_left->momentum().x();
  double pr_py0 = pr_prot_left->momentum().y();
  double pr_pz0 = pr_prot_left->momentum().z();
  
  double pr_vx0 = pr_vert_left->production_vertex()->position().x();//+BOPar.GetBeamDisplacementX()*1000.0;
  double pr_vy0 = pr_vert_left->production_vertex()->position().y();//+BOPar.GetBeamDisplacementY()*1000.0;
  double pr_vz0 = pr_vert_left->production_vertex()->position().z();//+BOPar.GetBeamDisplacementZ()*1000.0;
  
  double pr_ksi0 = BOPar.IPNonSmearedProtonMomentumToXi(pr_px0, pr_py0, pr_pz0);
  double pr_t0 = BOPar.IPNonSmearedProtonMomentumTot(pr_px0, pr_py0, pr_pz0);
  double log_minus_ksi0 = TMath::Log10(-pr_ksi0);
  double log_minus_t0 = TMath::Log10(-pr_t0);
  
  double pr_phi0 = ReduceAngleTo_0_2Pi(
      BOPar.ComputeNonSmearedProtonPhi(pr_prot_left->momentum()));
  
  double pr_theta_x0_no_cross_angle = pr_px0/TMath::Abs(pr_pz0);
  double pr_theta_y0_no_cross_angle = pr_py0/TMath::Abs(pr_pz0);
  double pr_theta_x0 = pr_theta_x0_no_cross_angle+BOPar.GetCrossingAngleX();  //the crossing angle !!
  double pr_theta_y0 = pr_theta_y0_no_cross_angle+BOPar.GetCrossingAngleY();  //the crossing angle !!
  
  double smeared_theta_x0 = pr_vert_left->momentum().x()/TMath::Abs(pr_vert_left->momentum().z());  //the crossing angle !!
  double smeared_theta_y0 = pr_vert_left->momentum().y()/TMath::Abs(pr_vert_left->momentum().z());  //the crossing angle !!
  double smeared_xi0 = BOPar.IPSmearedProtonMomentumToXi(
      pr_vert_left->momentum().x(), pr_vert_left->momentum().y(), 
      pr_vert_left->momentum().z());
  
  double cos_pr_phi0 = TMath::Cos(pr_phi0);
  double sin_pr_phi0 = TMath::Sin(pr_phi0);
  double pr_t_x0 = pr_t0*cos_pr_phi0*cos_pr_phi0;
  double pr_t_y0 = pr_t0*sin_pr_phi0*sin_pr_phi0;
  double log_minus_tx0 = TMath::Log10(-pr_t_x0);
  double log_minus_ty0 = TMath::Log10(-pr_t_y0);
  
  //right
  double pr_px1 = pr_prot_right->momentum().x();
  double pr_py1 = pr_prot_right->momentum().y();
  double pr_pz1 = pr_prot_right->momentum().z();
  
  double pr_vx1 = pr_vert_right->production_vertex()->position().x();//+BOPar.GetBeamDisplacementX()*1000.0;
  double pr_vy1 = pr_vert_right->production_vertex()->position().y();//+BOPar.GetBeamDisplacementY()*1000.0;
  double pr_vz1 = pr_vert_right->production_vertex()->position().z();//+BOPar.GetBeamDisplacementZ()*1000.0;
  
  double pr_ksi1 = BOPar.IPNonSmearedProtonMomentumToXi(pr_px1, pr_py1, pr_pz1);
  double pr_t1 = BOPar.IPNonSmearedProtonMomentumTot(pr_px1, pr_py1, pr_pz1);
  double log_minus_ksi1 = TMath::Log10(-pr_ksi1);
  double log_minus_t1 = TMath::Log10(-pr_t1);
  double pr_phi1 = ReduceAngleTo_0_2Pi(
      BOPar.ComputeNonSmearedProtonPhi(pr_prot_right->momentum()));
  
  double pr_theta_x1_no_cross_angle = pr_px1/TMath::Abs(pr_pz1);
  double pr_theta_y1_no_cross_angle = pr_py1/TMath::Abs(pr_pz1);
  double pr_theta_x1 = pr_theta_x1_no_cross_angle+BOPar.GetCrossingAngleX();  //the crossing angle !!
  double pr_theta_y1 = pr_theta_y1_no_cross_angle+BOPar.GetCrossingAngleY();  //the crossing angle !!

  double smeared_theta_x1 = pr_vert_right->momentum().x()/TMath::Abs(pr_vert_right->momentum().z());  //the crossing angle !!
  double smeared_theta_y1 = pr_vert_right->momentum().y()/TMath::Abs(pr_vert_right->momentum().z());  //the crossing angle !!
  double smeared_xi1 = BOPar.IPSmearedProtonMomentumToXi(
      pr_vert_right->momentum().x(), pr_vert_right->momentum().y(), 
      pr_vert_right->momentum().z());
  
  double cos_pr_phi1 = TMath::Cos(pr_phi1);
  double sin_pr_phi1 = TMath::Sin(pr_phi1);
  double pr_t_x1 = pr_t1*cos_pr_phi1*cos_pr_phi1;
  double pr_t_y1 = pr_t1*sin_pr_phi1*sin_pr_phi1;
  double log_minus_tx1 = TMath::Log10(-pr_t_x1);
  double log_minus_ty1 = TMath::Log10(-pr_t_y1);
  
  log_ksi_log_t_total_accepted_left_.Fill(log_minus_t0, log_minus_ksi0);
  log_ksi_log_t_total_accepted_right_.Fill(log_minus_t1, log_minus_ksi1);

  double mean_t = (pr_t0+pr_t1)/2.0;
  double log_minus_mean_t = TMath::Log10(-mean_t);
  elast_log_t_total_accepted_.Fill(log_minus_mean_t);
  
  double log_generated_dpe_mass = TMath::Log10(BOPar.HEPMCNonSmearedProtonstoM(pr_prot_right->momentum(), pr_prot_left->momentum()));
  log_dpe_mass_total_accepted_dist_.Fill(log_generated_dpe_mass);
  
  if(!rec_prot_pair->Valid()
//      || (rec_prot_pair->KsiLeft()>0.1) || (rec_prot_pair->KsiRight()>0.1)
      )
  {
    log_ksi_log_t_reconstruction_failed_left_.Fill(log_minus_t0, log_minus_ksi0);
    reconstruction_failure_vs_log_t_log_ksi_phi_left_.
          Fill(log_minus_t0, log_minus_ksi0, pr_phi0, 1.0);
    log_ksi_log_t_reconstruction_failed_right_.Fill(log_minus_t1, log_minus_ksi1);
    reconstruction_failure_vs_log_t_log_ksi_phi_right_.
          Fill(log_minus_t1, log_minus_ksi1, pr_phi1, 1.0);
    log_dpe_mass_reconstruction_failures_.Fill(log_generated_dpe_mass);
    return;
  }
  
  log_dpe_mass_reconstructed_dist_.Fill(log_generated_dpe_mass);
  
  log_ksi_log_t_reconstructed_protons_left_.Fill(log_minus_t0, log_minus_ksi0);
  log_ksi_log_t_reconstructed_protons_right_.Fill(log_minus_t1, log_minus_ksi1);
  
  
  double chi2overN = rec_prot_pair->Chi2Norm();
  double chi2 = rec_prot_pair->Chi2();
  int deg_of_freedom = rec_prot_pair->DegreesOfFreedom();
  double prob = TMath::Prob(chi2, deg_of_freedom);
  
  double x_err0 = rec_prot_pair->X3D() - pr_vx0;
  double y_err0 = rec_prot_pair->Y3D() - pr_vy0;
  double z_err0 = rec_prot_pair->Z3D() - pr_vz0;
  double thetax_err0 = rec_prot_pair->ThetaXAngleLeft() - pr_theta_x0;
  double thetay_err0 = rec_prot_pair->ThetaYAngleLeft() - pr_theta_y0;
  
  double smeared_thetax_err0 = rec_prot_pair->ThetaXAngleLeft() - smeared_theta_x0;
  double smeared_thetay_err0 = rec_prot_pair->ThetaYAngleLeft() - smeared_theta_y0;
  
  RPRecoProtMADXVariables mad_var_left = rec_prot_pair->GetMADXVariablesLeft();
  double t_reconst_left = BOPar.MADXCanonicalVariablesTot(mad_var_left);
  double t_err0 = t_reconst_left - pr_t0;
  double ksi_err0 = rec_prot_pair->KsiLeft() - pr_ksi0;
  double smeared_ksi_err0 = rec_prot_pair->KsiLeft() - smeared_xi0;
  double phi_reconst_left = 
    BOPar.MADXCanonicalVariablesToCrossingAngleCorrectedPhi(mad_var_left);
  double phi_err0 = AngleDiff(phi_reconst_left, pr_phi0);
  double cos_rec_phi0 = TMath::Cos(phi_reconst_left);
  double sin_rec_phi0 = TMath::Sin(phi_reconst_left);
  double rec_t_x0 = t_reconst_left*cos_rec_phi0*cos_rec_phi0;
  double rec_t_y0 = t_reconst_left*sin_rec_phi0*sin_rec_phi0;
  double t_x_err0 = rec_t_x0 - pr_t_x0;
  double t_y_err0 = rec_t_y0 - pr_t_y0;
  
  double x_err1 = rec_prot_pair->X3D() - pr_vx1;
  double y_err1 = rec_prot_pair->Y3D() - pr_vy1;
  double z_err1 = rec_prot_pair->Z3D() - pr_vz1;
  double thetax_err1 = rec_prot_pair->ThetaXAngleRight() - pr_theta_x1;
  double thetay_err1 = rec_prot_pair->ThetaYAngleRight() - pr_theta_y1;
  
  double smeared_thetax_err1 = rec_prot_pair->ThetaXAngleRight() - smeared_theta_x1;
  double smeared_thetay_err1 = rec_prot_pair->ThetaYAngleRight() - smeared_theta_y1;
  
  RPRecoProtMADXVariables mad_var_right = rec_prot_pair->GetMADXVariablesRight();
  double t_reconst_right = BOPar.MADXCanonicalVariablesTot(mad_var_right);
  double t_err1 = t_reconst_right - pr_t1;
  double ksi_err1 = rec_prot_pair->KsiRight() - pr_ksi1;
  double smeared_ksi_err1 = rec_prot_pair->KsiRight() - smeared_xi1;
  double phi_reconst_right = 
    BOPar.MADXCanonicalVariablesToCrossingAngleCorrectedPhi(mad_var_right);
  double phi_err1 = AngleDiff(phi_reconst_right, pr_phi1);
  double cos_rec_phi1 = TMath::Cos(phi_reconst_right);
  double sin_rec_phi1 = TMath::Sin(phi_reconst_right);
  double rec_t_x1 = t_reconst_right*cos_rec_phi1*cos_rec_phi1;
  double rec_t_y1 = t_reconst_right*sin_rec_phi1*sin_rec_phi1;
  double t_x_err1 = rec_t_x1 - pr_t_x1;
  double t_y_err1 = rec_t_y1 - pr_t_y1;
  
  //fill uppertail probability cut
  if(prob<chi2_upper_tail_prob_cut_)
    return;
  
  //fill elastic scattering histograms
//  double mean_t = BOPar.MADXCanonicalVariablesToElasticPhysics_t(mad_var_right, mad_var_left);
//  double log_minus_mean_t = TMath::Log10(-mean_t);
  elast_log_t_accepted_.Fill(log_minus_mean_t);
  
  elast_left_thx_rec_error_.Fill(thetax_err0);
  elast_left_thy_rec_error_.Fill(thetay_err0);
  elast_right_thx_rec_error_.Fill(thetax_err1);
  elast_right_thy_rec_error_.Fill(thetay_err1);
  
  double elastic_thx = BOPar.MADXCanonicalVariablesToElasticPhysicsThx(mad_var_right, mad_var_left);
  double elastic_thy = BOPar.MADXCanonicalVariablesToElasticPhysicsThy(mad_var_right, mad_var_left);
  double elastic_thx_generated_mean = (pr_theta_x1_no_cross_angle - pr_theta_x0_no_cross_angle)/2.0;
  double elastic_thy_generated_mean = (pr_theta_y1_no_cross_angle - pr_theta_y0_no_cross_angle)/2.0;
  double elastic_thx_rec_error = elastic_thx - elastic_thx_generated_mean;
  double elastic_thy_rec_error = elastic_thy - elastic_thy_generated_mean;
  elast_scattering_thx_rec_error_.Fill(elastic_thx_rec_error);
  elast_scattering_thy_rec_error_.Fill(elastic_thy_rec_error);
  
  double elastic_thx_right_left_rec_err = (rec_prot_pair->ThetaXAngleRight() - BOPar.GetCrossingAngleX()) + (rec_prot_pair->ThetaXAngleLeft() - BOPar.GetCrossingAngleX()); 
  double elastic_thy_right_left_rec_err = (rec_prot_pair->ThetaYAngleRight() - BOPar.GetCrossingAngleY()) + (rec_prot_pair->ThetaYAngleLeft() - BOPar.GetCrossingAngleY());
  elast_thx_left_right_difference_.Fill(-elastic_thx_right_left_rec_err);
  elast_thy_left_right_difference_.Fill(-elastic_thy_right_left_rec_err);
  
  elast_left_xi_error_.Fill(ksi_err0);
  elast_right_xi_error_.Fill(ksi_err1);
  
  elast_x0_rec_error_.Fill(x_err0);
  elast_y0_rec_error_.Fill(y_err0);
  elast_z0_rec_error_.Fill(z_err0);
    
  elast_chi2ndf_dist_.Fill(chi2overN);
  elast_probe_dist_.Fill(prob);

  //fill diffractive scattering histograms
  ksi_error_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, ksi_err0);
  smeared_ksi_error_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, smeared_ksi_err0);
  t_error_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, t_err0);
  phi_error_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, phi_err0);
  vx_error_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, x_err0);
  vy_error_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, y_err0);
  vz_error_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, z_err0);
  thetax_error_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, thetax_err0);
  thetay_error_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, thetay_err0);
  smeared_thetax_error_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, smeared_thetax_err0);
  smeared_thetay_error_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, smeared_thetay_err0);
  ksi_rel_error_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, ksi_err0/pr_ksi0);
  t_rel_error_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, t_err0/pr_t0);
  chi2_error_overn_N_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, chi2overN);
  prob_function_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, prob);
  t_x_error_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, t_x_err0);
  t_y_error_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, t_y_err0);
  t_x_rel_error_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, t_x_err0/pr_t_x0);
  t_y_rel_error_vs_log_t_log_ksi_phi_left_.Fill(log_minus_t0, log_minus_ksi0, pr_phi0, t_y_err0/pr_t_y0);

  t_x_error_vs_log_tx_log_ksi_phi_left_.Fill(log_minus_tx0, log_minus_ksi0, pr_phi0, t_x_err0);
  t_y_error_vs_log_ty_log_ksi_phi_left_.Fill(log_minus_ty0, log_minus_ksi0, pr_phi0, t_y_err0);
  t_x_rel_error_vs_log_tx_log_ksi_phi_left_.Fill(log_minus_tx0, log_minus_ksi0, pr_phi0, t_x_err0/pr_t_x0);
  t_y_rel_error_vs_log_ty_log_ksi_phi_left_.Fill(log_minus_ty0, log_minus_ksi0, pr_phi0, t_y_err0/pr_t_y0);

  
  ksi_error_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, ksi_err1);
  smeared_ksi_error_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, smeared_ksi_err1);
  t_error_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, t_err1);
  phi_error_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, phi_err1);
  vx_error_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, x_err1);
  vy_error_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, y_err1);
  vz_error_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, z_err1);
  thetax_error_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, thetax_err1);
  thetay_error_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, thetay_err1);
  smeared_thetax_error_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, smeared_thetax_err1);
  smeared_thetay_error_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, smeared_thetay_err1);
  ksi_rel_error_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, ksi_err1/pr_ksi1);
  t_rel_error_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, t_err1/pr_t1);
  chi2_error_overn_N_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, chi2overN);
  prob_function_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, prob);
  t_x_error_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, t_x_err1);
  t_y_error_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, t_y_err1);
  t_x_rel_error_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, t_x_err1/pr_t_x1);
  t_y_rel_error_vs_log_t_log_ksi_phi_right_.Fill(log_minus_t1, log_minus_ksi1, pr_phi1, t_y_err1/pr_t_y1);

  t_x_error_vs_log_tx_log_ksi_phi_right_.Fill(log_minus_tx1, log_minus_ksi1, pr_phi1, t_x_err1);
  t_y_error_vs_log_ty_log_ksi_phi_right_.Fill(log_minus_ty1, log_minus_ksi1, pr_phi1, t_y_err1);
  t_x_rel_error_vs_log_tx_log_ksi_phi_right_.Fill(log_minus_tx1, log_minus_ksi1, pr_phi1, t_x_err1/pr_t_x1);
  t_y_rel_error_vs_log_ty_log_ksi_phi_right_.Fill(log_minus_ty1, log_minus_ksi1, pr_phi1, t_y_err1/pr_t_y1);
  
  double ksi_rec_lo, ksi_rec_hi;
//  double t_rec_lo, t_rec_hi, phi_rec_lo, phi_rec_hi, phi_pr_lo, phi_pr_hi;
  double ksi_pr_lo, ksi_pr_hi, t_pr_lo, t_pr_hi;
  double ksi_err_lo, ksi_err_hi, t_err_lo, t_err_hi, phi_err_lo, phi_err_hi;
  
  if(TMath::Abs(pr_ksi0)>TMath::Abs(pr_ksi1))
  {
    ksi_pr_hi = pr_ksi0;
    ksi_pr_lo = pr_ksi1;
    t_pr_hi = pr_t0;
    t_pr_lo = pr_t1;
//    phi_pr_hi = pr_phi0;
//    phi_pr_lo = pr_phi1;
    
    ksi_rec_hi = rec_prot_pair->KsiLeft();
    ksi_rec_lo = rec_prot_pair->KsiRight();
//    t_rec_hi = t_reconst_left;
//    t_rec_lo = t_reconst_right;
//    phi_rec_hi = phi_reconst_left;
//    phi_rec_lo = phi_reconst_right;
    
    ksi_err_hi = ksi_err0; 
    ksi_err_lo = ksi_err1;
    t_err_hi = t_err0;
    t_err_lo = t_err1;
    phi_err_hi = phi_err0;
    phi_err_lo = phi_err1;
  }
  else
  {
    ksi_pr_hi = pr_ksi1;
    ksi_pr_lo = pr_ksi0;
    t_pr_hi = pr_t1;
    t_pr_lo = pr_t0;
//    phi_pr_hi = pr_phi1;
//    phi_pr_lo = pr_phi0;
    
    ksi_rec_hi = rec_prot_pair->KsiRight();
    ksi_rec_lo = rec_prot_pair->KsiLeft();
//    t_rec_hi = t_reconst_right;
//    t_rec_lo = t_reconst_left;
//    phi_rec_hi = phi_reconst_right;
//    phi_rec_lo = phi_reconst_left;
    
    ksi_err_hi = ksi_err1;
    ksi_err_lo = ksi_err0;
    t_err_hi = t_err1;
    t_err_lo = t_err0;
    phi_err_hi = phi_err1;
    phi_err_lo = phi_err0;
  }
  
  double M_rec = BOPar.MADXCanonicalVariablesToM(mad_var_right, mad_var_left);
  double M_pr = BOPar.HEPMCNonSmearedProtonstoM(pr_prot_right->momentum(), pr_prot_left->momentum());
  double M_error = M_rec - M_pr;
  double Log_M = TMath::Log10(M_pr);
  double Log_ksi_lo_div_ksi_hi = TMath::Log10(ksi_pr_lo/ksi_pr_hi);
  double Log_ksi_hi = TMath::Log10(-ksi_pr_hi);
  
  if(ksi_rec_hi*ksi_rec_lo<=0.0 || ksi_pr_hi*ksi_pr_lo<=0.0)
    return;
  
  if(verbosity)
  {
    std::cout<<" M_pr="<<M_pr<<" ksi_pr_hi="<<ksi_pr_hi<<" ksi_pr_lo="<<ksi_pr_lo<<" t_pr_hi="<<t_pr_hi<<" t_pr_lo="<<t_pr_lo<<std::endl;
    std::cout<<" M_error="<<M_error<<" ksi_err_hi="<<ksi_err_hi<<" ksi_err_lo="<<ksi_err_lo<<std::endl;
    std::cout<<" t_err_hi="<<t_err_hi<<" t_err_lo="<<t_err_lo<<std::endl;
    std::cout<<" x_err0="<<x_err0<<" y_err0="<<y_err0<<" z_err0="<<z_err0<<std::endl;
  }
  
  M_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, M_error);
  M_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, M_error/M_pr);
  ksi_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, ksi_err_hi);
  ksi_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, ksi_err_lo);
  ksi_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, ksi_err_hi);
  ksi_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, ksi_err_lo);
  
  ksi_higher_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, ksi_err_hi/ksi_pr_hi);
  ksi_lower_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, ksi_err_lo/ksi_pr_lo);
  ksi_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, ksi_err_hi/ksi_pr_hi);
  ksi_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, ksi_err_lo/ksi_pr_lo);
  
  t_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, t_err_hi);
  t_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, t_err_lo);
  t_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, t_err_hi);
  t_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, t_err_lo);

  t_higher_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, t_err_hi/t_pr_hi);
  t_lower_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, t_err_lo/t_pr_lo);
  t_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, t_err_hi/t_pr_hi);
  t_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, t_err_lo/t_pr_lo);
  
  phi_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, phi_err_hi);
  phi_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, phi_err_lo);
  phi_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, phi_err_hi);
  phi_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, phi_err_lo);
    
  x_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, x_err0);
  y_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, y_err0);
  z_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, z_err0);
    
  probe_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, prob);
  chsqndf_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_.Fill(Log_M, Log_ksi_lo_div_ksi_hi, Log_ksi_hi, chi2overN);
}


double RPProtonPairReconstructionInfo::ReduceAngleTo_0_2Pi(double angle)
{
  double TwoPi = 2*TMath::Pi();
  while(angle>=TwoPi)
  {
    angle -= TwoPi;
  }
  while(angle<0)
  {
    angle += TwoPi;
  }
  return angle;
}


double RPProtonPairReconstructionInfo::ReduceAngleTo_Pi_Pi(double angle)
{
  double TwoPi = 2*TMath::Pi();
  double Pi = TMath::Pi();
  while(angle>=Pi)
  {
    angle -= TwoPi;
  }
  while(angle<-Pi)
  {
    angle += TwoPi;
  }
  return angle;
}


double RPProtonPairReconstructionInfo::AngleDiff(double angle1, double angle2)
{
  double diff = ReduceAngleTo_Pi_Pi(angle2-angle1);
  
  return diff;
}


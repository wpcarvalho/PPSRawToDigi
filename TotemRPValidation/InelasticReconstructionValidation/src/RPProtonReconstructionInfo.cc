#include "TotemRPValidation/InelasticReconstructionValidation/interface/RPProtonReconstructionInfo.h"
#include "TProfile3D.h"
#include "TMath.h"
#include "TH2F.h"
#include "TotemRPValidation/ValidationTools/interface/ReconstructionProfile.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


RPProtonReconstructionInfo::RPProtonReconstructionInfo(const std::string &path, 
      int arm_id, const edm::ParameterSet& conf)
 : BaseHistogramManager(path), arm_id_(arm_id)
{
  chi2_upper_tail_prob_cut_ = conf.getParameter<double>("Chi2UpperTailProbCut");
  InitializeHistograms();
}


void RPProtonReconstructionInfo::InitializeHistograms()
{  
  char name[1024];
  
  sprintf(name, "log_ksi_log_t_distribution_%i", arm_id_);
  log_ksi_log_t_distribution_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_distribution_.SetDirectory(0);
  log_ksi_log_t_distribution_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_distribution_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_distribution_);
  
  sprintf(name, "log_ksi_log_t_reconstructed_protons_%i", arm_id_);
  log_ksi_log_t_reconstructed_protons_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_reconstructed_protons_.SetDirectory(0);
  log_ksi_log_t_reconstructed_protons_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_reconstructed_protons_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_reconstructed_protons_);
  
  sprintf(name, "log_ksi_log_t_acceptance_%i", arm_id_);
  log_ksi_log_t_acceptance_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_acceptance_.SetDirectory(0);
  log_ksi_log_t_acceptance_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_acceptance_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_acceptance_);
  
  sprintf(name, "log_ksi_log_t_reconstruction_failed_%i", arm_id_);
  log_ksi_log_t_reconstruction_failed_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_reconstruction_failed_.SetDirectory(0);
  log_ksi_log_t_reconstruction_failed_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_reconstruction_failed_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_reconstruction_failed_);
  
  sprintf(name, "log_ksi_log_t_reconstruction_failure_rate_%i", arm_id_);
  log_ksi_log_t_reconstruction_failure_rate_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_reconstruction_failure_rate_.SetDirectory(0);
  log_ksi_log_t_reconstruction_failure_rate_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_reconstruction_failure_rate_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_reconstruction_failure_rate_);
  
  sprintf(name, "log_ksi_log_t_total_accepted_%i", arm_id_);
  log_ksi_log_t_total_accepted_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_total_accepted_.SetDirectory(0);
  log_ksi_log_t_total_accepted_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_total_accepted_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_total_accepted_);
  
  sprintf(name, "log_ksi_log_t_total_acceptance_%i", arm_id_);
  log_ksi_log_t_total_acceptance_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_total_acceptance_.SetDirectory(0);
  log_ksi_log_t_total_acceptance_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_total_acceptance_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_total_acceptance_);
  
  sprintf(name, "log_ksi_log_t_reconstruction_failure_to_total_accepted_%i", arm_id_);
  log_ksi_log_t_reconstruction_failure_to_total_accepted_ = TH2F(name, name, 60, -4, 2, 70, -7, 0.);
  log_ksi_log_t_reconstruction_failure_to_total_accepted_.SetDirectory(0);
  log_ksi_log_t_reconstruction_failure_to_total_accepted_.SetXTitle("Log(-t/GeV^{2})");
  log_ksi_log_t_reconstruction_failure_to_total_accepted_.SetYTitle("Log(-#xi)");
  RegisterHistogram(log_ksi_log_t_reconstruction_failure_to_total_accepted_);

  sprintf(name, "log_ksi_distribution_%i", arm_id_);
  log_ksi_distribution_ = TH1F(name, name, 1000, -8, 0.);
  log_ksi_distribution_.SetDirectory(0);
  log_ksi_distribution_.SetXTitle("Log(-#xi)");
  log_ksi_distribution_.SetYTitle("Entries");
  RegisterHistogram(log_ksi_distribution_);

  sprintf(name, "log_ksi_reconstructed_protons_%i", arm_id_);
  log_ksi_reconstructed_protons_ = TH1F(name, name, 1000, -8, 0.);
  log_ksi_reconstructed_protons_.SetDirectory(0);
  log_ksi_reconstructed_protons_.SetXTitle("Log(-#xi)");
  log_ksi_reconstructed_protons_.SetYTitle("Rec. protons");
  RegisterHistogram(log_ksi_reconstructed_protons_);

  sprintf(name, "log_ksi_acceptance_%i", arm_id_);
  log_ksi_acceptance_ = TH1F(name, name, 1000, -8, 0.);
  log_ksi_acceptance_.SetDirectory(0);
  log_ksi_acceptance_.SetXTitle("Log(-#xi)");
  log_ksi_acceptance_.SetYTitle("Rec. acceptance");
  RegisterHistogram(log_ksi_acceptance_);

  sprintf(name, "log_ksi_total_accepted_%i", arm_id_);
  log_ksi_total_accepted_ = TH1F(name, name, 1000, -8, 0.);
  log_ksi_total_accepted_.SetDirectory(0);
  log_ksi_total_accepted_.SetXTitle("Log(-#xi)");
  log_ksi_total_accepted_.SetYTitle("Total accepted");
  RegisterHistogram(log_ksi_total_accepted_);

  sprintf(name, "log_ksi_total_acceptance_%i", arm_id_);
  log_ksi_total_acceptance_ = TH1F(name, name, 1000, -8, 0.);
  log_ksi_total_acceptance_.SetDirectory(0);
  log_ksi_total_acceptance_.SetXTitle("Log(-#xi)");
  log_ksi_total_acceptance_.SetYTitle("Total acceptance");
  RegisterHistogram(log_ksi_total_acceptance_);
  
  sprintf(name, "ksi_error_vs_log_t_log_ksi_phi_%i", arm_id_);
  ksi_error_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s"); 
  //ksi_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(ksi_error_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "smeared_ksi_error_vs_log_t_log_ksi_phi_%i", arm_id_);
  smeared_ksi_error_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s"); 
  //ksi_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(smeared_ksi_error_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "t_error_vs_log_t_log_ksi_phi_%i", arm_id_);
  t_error_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //t_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(t_error_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "phi_error_vs_log_t_log_ksi_phi_%i", arm_id_);
  phi_error_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //phi_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(phi_error_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "vx_error_vs_log_t_log_ksi_phi_%i", arm_id_);
  vx_error_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //vx_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(vx_error_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "vy_error_vs_log_t_log_ksi_phi_%i", arm_id_);
  vy_error_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //vy_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(vy_error_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "thetax_error_vs_log_t_log_ksi_phi_%i", arm_id_);
  thetax_error_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //thetax_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(thetax_error_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "smeared_thetax_error_vs_log_t_log_ksi_phi_%i", arm_id_);
  smeared_thetax_error_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //smeared_thetax_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(smeared_thetax_error_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "thetay_error_vs_log_t_log_ksi_phi_%i", arm_id_);
  thetay_error_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //thetay_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(thetay_error_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "smeared_thetay_error_vs_log_t_log_ksi_phi_%i", arm_id_);
  smeared_thetay_error_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //smeared_thetay_error_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(smeared_thetay_error_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "chi2_error_overn_N_vs_log_t_log_ksi_phi_%i", arm_id_);
  chi2_error_overn_N_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  //chi2_error_overn_N_vs_log_t_log_ksi_phi_.SetDirectory(0);
  RegisterHistogram(chi2_error_overn_N_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "ksi_rel_error_vs_log_t_log_ksi_phi_%i", arm_id_);
  ksi_rel_error_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(ksi_rel_error_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "t_rel_error_vs_log_t_log_ksi_phi_%i", arm_id_);
  t_rel_error_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_rel_error_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "t_x_error_vs_log_t_log_ksi_phi_%i", arm_id_);
  t_x_error_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_x_error_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "t_x_error_vs_log_tx_log_ksi_phi_%i", arm_id_);
  t_x_error_vs_log_tx_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_x_error_vs_log_tx_log_ksi_phi_);
  
  sprintf(name, "t_y_error_vs_log_t_log_ksi_phi_%i", arm_id_);
  t_y_error_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_y_error_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "t_y_error_vs_log_ty_log_ksi_phi_%i", arm_id_);
  t_y_error_vs_log_ty_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_y_error_vs_log_ty_log_ksi_phi_);
  
  sprintf(name, "t_x_rel_error_vs_log_t_log_ksi_phi_%i", arm_id_);
  t_x_rel_error_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_x_rel_error_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "t_x_rel_error_vs_log_tx_log_ksi_phi_%i", arm_id_);
  t_x_rel_error_vs_log_tx_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_x_rel_error_vs_log_tx_log_ksi_phi_);
  
  sprintf(name, "t_y_rel_error_vs_log_t_log_ksi_phi_%i", arm_id_);
  t_y_rel_error_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_y_rel_error_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "t_y_rel_error_vs_log_ty_log_ksi_phi_%i", arm_id_);
  t_y_rel_error_vs_log_ty_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(t_y_rel_error_vs_log_ty_log_ksi_phi_);
  
  sprintf(name, "reconstruction_failure_vs_log_t_log_ksi_phi_%i", arm_id_);
  reconstruction_failure_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(reconstruction_failure_vs_log_t_log_ksi_phi_);
  
  sprintf(name, "prob_function_vs_log_t_log_ksi_phi_%i", arm_id_);
  prob_function_vs_log_t_log_ksi_phi_ = 
      ReconstructionProfile(name, name, 75, -5.5, 2., 60, -5, 1, 50, 0, 2*TMath::Pi(), "s");
  RegisterHistogram(prob_function_vs_log_t_log_ksi_phi_);
}


void RPProtonReconstructionInfo::Finalize()
{
  log_ksi_log_t_acceptance_.Divide(&log_ksi_log_t_reconstructed_protons_, &log_ksi_log_t_distribution_);
  log_ksi_log_t_reconstruction_failure_rate_.Divide(&log_ksi_log_t_reconstruction_failed_, &log_ksi_log_t_distribution_);
  log_ksi_log_t_reconstruction_failure_to_total_accepted_.Divide(&log_ksi_log_t_reconstruction_failed_, &log_ksi_log_t_total_accepted_);
  log_ksi_log_t_total_acceptance_.Divide(&log_ksi_log_t_total_accepted_, &log_ksi_log_t_distribution_);
  log_ksi_total_acceptance_.Divide(&log_ksi_total_accepted_, &log_ksi_distribution_);
  log_ksi_acceptance_.Divide(&log_ksi_reconstructed_protons_, &log_ksi_distribution_);
}


void RPProtonReconstructionInfo::FillReferenceHistograms(
    const BeamOpticsParams &BOPar, const HepMC::GenParticle *pr_prot, 
    const HepMC::GenParticle *pr_vert, bool verbosity)
{
  if(!pr_prot)
    return;
  
//  std::cout<<"RPProtonReconstructionInfo::FillReferenceHistograms"<<std::endl;
  
  double pr_px = pr_prot->momentum().x();
  double pr_py = pr_prot->momentum().y();
  double pr_pz = pr_prot->momentum().z();
  
  double pr_ksi = BOPar.IPNonSmearedProtonMomentumToXi(pr_px, pr_py, pr_pz);
  double pr_t = BOPar.IPNonSmearedProtonMomentumTot(pr_px, pr_py, pr_pz);
  
//  std::cout<<"primary proton:"<<" pr_px="<<pr_px<<" pr_py="<<pr_py<<" pr_pz="<<pr_pz<<std::endl;
//  std::cout<<"mad_thx="<<pr_px/7000<<" mad_thy="<<pr_py/7000<<std::endl;
//
//  std::cout<<"pr_ksi = "<<pr_ksi<<std::endl;
//  std::cout<<"pr_t = "<<pr_t<<std::endl;
        
  double log_minus_ksi = TMath::Log10(-pr_ksi);
  double log_minus_t = TMath::Log10(-pr_t);

//  std::cout<<"log_minus_ksi = "<<log_minus_ksi<<std::endl;
//  std::cout<<"log_minus_t = "<<log_minus_t<<std::endl;
        
  log_ksi_log_t_distribution_.Fill(log_minus_t, log_minus_ksi);
  log_ksi_distribution_.Fill(log_minus_ksi);
//  std::cout<<"RPProtonReconstructionInfo::FillReferenceHistograms, ended"<<std::endl;
}


void RPProtonReconstructionInfo::FillResidualHistograms(
      const BeamOpticsParams &BOPar, 
      const HepMC::GenParticle *pr_prot, const HepMC::GenParticle *pr_vert, 
      const RPReconstructedProton* rec_prot, bool verbosity)
{
  if(!pr_prot || !rec_prot)
    return;

//  std::cout<<"RPProtonReconstructionInfo::FillResidualHistograms"<<std::endl;
  
  double pr_px = pr_prot->momentum().x();
  double pr_py = pr_prot->momentum().y();
  double pr_pz = pr_prot->momentum().z();
  
  double pr_vx = pr_vert->production_vertex()->position().x();//+BOPar.GetBeamDisplacementX()*1000.0;
  double pr_vy = pr_vert->production_vertex()->position().y();//+BOPar.GetBeamDisplacementY()*1000.0;
  double pr_vz = pr_vert->production_vertex()->position().z();//+BOPar.GetBeamDisplacementZ()*1000.0;
  
  //project the primary vertex on z=0 along the primary proton momentum direction
  pr_vx -= pr_px/pr_pz*pr_vz;
  pr_vy -= pr_py/pr_pz*pr_vz;
  
  double pr_ksi = BOPar.IPNonSmearedProtonMomentumToXi(pr_px, pr_py, pr_pz);
  double pr_t = BOPar.IPNonSmearedProtonMomentumTot(pr_px, pr_py, pr_pz);

//  std::cout<<"pr_ksi = "<<pr_ksi<<std::endl;
//  std::cout<<"pr_t = "<<pr_t<<std::endl;
  
  double log_minus_ksi = TMath::Log10(-pr_ksi);
  double log_minus_t = TMath::Log10(-pr_t);

//  std::cout<<"log_minus_ksi = "<<log_minus_ksi<<std::endl;
//  std::cout<<"log_minus_t = "<<log_minus_t<<std::endl;
  
  double pr_phi = ReduceAngleTo_0_2Pi(
      BOPar.ComputeNonSmearedProtonPhi(pr_prot->momentum()));
  double cos_pr_phi = TMath::Cos(pr_phi);
  double sin_pr_phi = TMath::Sin(pr_phi);
  double pr_t_x = pr_t*cos_pr_phi*cos_pr_phi;
  double pr_t_y = pr_t*sin_pr_phi*sin_pr_phi;
  double pr_theta_x = pr_px/TMath::Abs(pr_pz)+BOPar.GetCrossingAngleX();  //the angle does not include the crossing angle !!
  double pr_theta_y = pr_py/TMath::Abs(pr_pz)+BOPar.GetCrossingAngleY();  //the angle does not include the crossing angle !!
  
  double smeared_theta_x = pr_vert->momentum().x()/TMath::Abs(pr_vert->momentum().z()); //crossing angle included
  double smeared_theta_y = pr_vert->momentum().y()/TMath::Abs(pr_vert->momentum().z());  //crossing angle included
  double smeared_xi = BOPar.IPSmearedProtonMomentumToXi(
      pr_vert->momentum().x(), pr_vert->momentum().y(), pr_vert->momentum().z());
  
  double smeared_theta_x_error = rec_prot->Theta_x_angle() - smeared_theta_x;
  double smeared_theta_y_error = rec_prot->Theta_y_angle() - smeared_theta_y;
  double smeared_xi_error = rec_prot->Ksi() - smeared_xi;

  double log_minus_tx = TMath::Log10(-pr_t_x);
  double log_minus_ty = TMath::Log10(-pr_t_y);
  
  log_ksi_log_t_total_accepted_.Fill(log_minus_t, log_minus_ksi);
  log_ksi_total_accepted_.Fill(log_minus_ksi);
  
  if(!rec_prot->Valid()
//      || !(rec_prot->Ksi()<0.01) 
  )
  {
    log_ksi_log_t_reconstruction_failed_.Fill(log_minus_t, log_minus_ksi);
    reconstruction_failure_vs_log_t_log_ksi_phi_.
          Fill(log_minus_t, log_minus_ksi, pr_phi, 1.0);
    return;
  }
  
  log_ksi_log_t_reconstructed_protons_.Fill(log_minus_t, log_minus_ksi);
  log_ksi_reconstructed_protons_.Fill(log_minus_ksi);
  
  double x_err = rec_prot->X() - pr_vx;
  double y_err = rec_prot->Y() - pr_vy;
  double thetax_err = rec_prot->Theta_x_angle() - pr_theta_x;  //correct for the crossing angle
  double thetay_err = rec_prot->Theta_y_angle() - pr_theta_y;  //correct for the crossing angle
  
  RPRecoProtMADXVariables mad_var = rec_prot->GetMADXVariables();
  double reconstructed_t = BOPar.MADXCanonicalVariablesTot(mad_var);
  double t_err = reconstructed_t - pr_t;
  double ksi_err = rec_prot->Ksi() - pr_ksi;
  double reconstructed_phi = BOPar.MADXCanonicalVariablesToCrossingAngleCorrectedPhi(mad_var);
  double phi_err = AngleDiff(reconstructed_phi, pr_phi);
  
  double cos_rec_phi = TMath::Cos(reconstructed_phi);
  double sin_rec_phi = TMath::Sin(reconstructed_phi);
  double rec_t_x = reconstructed_t*cos_rec_phi*cos_rec_phi;
  double rec_t_y = reconstructed_t*sin_rec_phi*sin_rec_phi;
  double t_x_err = rec_t_x - pr_t_x;
  double t_y_err = rec_t_y - pr_t_y;
  
  double chi2overN = rec_prot->Chi2Norm();
  double chi2 = rec_prot->Chi2();
  int deg_of_freedom = rec_prot->DegreesOfFreedom();
  double prob = TMath::Prob(chi2, deg_of_freedom);
  
  if(verbosity)
  {
    std::cout<<"Reconstruction validation summary"<<std::endl;
    std::cout<<"================================="<<std::endl;
    std::cout<<"Smeared proton - reconstructed:"<<std::endl;
    std::cout<<"  x_err="<<x_err<<" rec="<<rec_prot->X()<<" true="<<pr_vx<<std::endl;
    std::cout<<"  y_err="<<y_err<<" rec="<<rec_prot->Y()<<" true="<<pr_vy<<std::endl;
    std::cout<<"  smeared_theta_x_error="<<smeared_theta_x_error<<" rec="<<rec_prot->Theta_x_angle()<<" true="<<smeared_theta_x<<std::endl;
    std::cout<<"  smeared_theta_y_error="<<smeared_theta_y_error<<" rec="<<rec_prot->Theta_y_angle()<<" true="<<smeared_theta_y<<std::endl;
    std::cout<<"  smeared_xi_error="<<smeared_xi_error<<" rec="<<rec_prot->Ksi()<<" true="<<smeared_xi<<std::endl;
    std::cout<<"Original proton - reconstructed:"<<std::endl;
    std::cout<<"  thetax_err="<<thetax_err<<" rec="<<rec_prot->Theta_x_angle()<<" true+cra="<<pr_theta_x<<std::endl;
    std::cout<<"  thetay_err="<<thetay_err<<" rec="<<rec_prot->Theta_y_angle()<<" true+cra="<<pr_theta_y<<std::endl;
    std::cout<<"  t_err="<<t_err<<" rec="<<reconstructed_t<<" true="<<pr_t<<std::endl;
    std::cout<<"  ksi_err="<<ksi_err<<" rec="<<rec_prot->Ksi()<<" true="<<pr_ksi<<std::endl;
    std::cout<<"  phi_err="<<phi_err<<" rec="<<reconstructed_phi<<" true="<<pr_phi<<std::endl;
    std::cout<<"  chi2overN="<<chi2overN<<std::endl;
  }
  
  if(prob<chi2_upper_tail_prob_cut_)
    return;
  
  ksi_error_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, ksi_err);
  smeared_ksi_error_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, smeared_xi_error);
  t_error_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, t_err);
  phi_error_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, phi_err);
  vx_error_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, x_err);
  vy_error_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, y_err);
  thetax_error_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, thetax_err);
  thetay_error_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, thetay_err);
  smeared_thetax_error_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, smeared_theta_x_error);
  smeared_thetay_error_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, smeared_theta_y_error);
  chi2_error_overn_N_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, chi2overN);
  prob_function_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, prob);
  
  ksi_rel_error_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, ksi_err/pr_ksi);
  t_rel_error_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, t_err/pr_t);
  
  t_x_error_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, t_x_err);
  t_y_error_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, t_y_err);
  t_x_rel_error_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, t_x_err/pr_t_x);
  t_y_rel_error_vs_log_t_log_ksi_phi_.Fill(log_minus_t, log_minus_ksi, pr_phi, t_y_err/pr_t_y);
  
  t_x_error_vs_log_tx_log_ksi_phi_.Fill(log_minus_tx, log_minus_ksi, pr_phi, t_x_err);
  t_y_error_vs_log_ty_log_ksi_phi_.Fill(log_minus_ty, log_minus_ksi, pr_phi, t_y_err);
  t_x_rel_error_vs_log_tx_log_ksi_phi_.Fill(log_minus_tx, log_minus_ksi, pr_phi, t_x_err/pr_t_x);
  t_y_rel_error_vs_log_ty_log_ksi_phi_.Fill(log_minus_ty, log_minus_ksi, pr_phi, t_y_err/pr_t_y);

//  std::cout<<"RPProtonReconstructionInfo::FillResidualHistograms, ended"<<std::endl;
}


double RPProtonReconstructionInfo::ReduceAngleTo_0_2Pi(double angle)
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


double RPProtonReconstructionInfo::ReduceAngleTo_Pi_Pi(double angle)
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


double RPProtonReconstructionInfo::AngleDiff(double angle1, double angle2)
{
  double diff = ReduceAngleTo_Pi_Pi(angle2-angle1);
  
  return diff;
}


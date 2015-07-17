#ifndef TotemRPValidation_RPReconstructedTracksValidation_RPProtonPairReconstructionInfo_h
#define TotemRPValidation_RPReconstructedTracksValidation_RPProtonPairReconstructionInfo_h

#include "TotemRPValidation/BaseValidationClasses/interface/BaseHistogramManager.h"
#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile3D.h"
#include "TotemRPValidation/ValidationTools/interface/ReconstructionProfile.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProtonPair.h"
#include "HepMC/GenEvent.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class RPProtonPairReconstructionInfo : public BaseHistogramManager
{
  public:
    RPProtonPairReconstructionInfo(const std::string &path, const edm::ParameterSet& conf);
    
    void FillReferenceHistograms(const BeamOpticsParams &BOPar, 
          const HepMC::GenParticle *pr_prot_left, 
          const HepMC::GenParticle *pr_prot_right, 
          const HepMC::GenParticle *pr_vert_left, 
          const HepMC::GenParticle *pr_vert_right, bool verbosity=false);
    void FillResidualHistograms(const BeamOpticsParams &BOPar, 
          const HepMC::GenParticle *pr_prot_left,
          const HepMC::GenParticle *pr_prot_right,
          const HepMC::GenParticle *pr_vert_left,
          const HepMC::GenParticle *pr_vert_right,
          const RPReconstructedProtonPair* rec_prot_pair, bool verbosity=false);
    
  private:
    void InitializeHistograms();
    virtual void Finalize();
    double ReduceAngleTo_0_2Pi(double angle);
    double ReduceAngleTo_Pi_Pi(double angle);
    double AngleDiff(double angle1, double angle2);
    double chi2_upper_tail_prob_cut_;
    
//    double mp_;
//    double p0_;
//    double E1_;
    
    TH2F log_ksi_log_t_distribution_left_;
    TH2F log_ksi_log_t_distribution_right_;
    TH2F log_ksi_log_t_reconstructed_protons_left_;
    TH2F log_ksi_log_t_reconstructed_protons_right_;
    TH2F log_ksi_log_t_reconstruction_failed_left_;
    TH2F log_ksi_log_t_reconstruction_failed_right_;
    TH2F log_ksi_log_t_reconstruction_failure_rate_left_;
    TH2F log_ksi_log_t_reconstruction_failure_rate_right_;
    TH2F log_ksi_log_t_total_accepted_left_;
    TH2F log_ksi_log_t_total_accepted_right_;
    TH2F log_ksi_log_t_reconstruction_failure_to_total_accepted_left_;
    TH2F log_ksi_log_t_reconstruction_failure_to_total_accepted_right_;
    TH2F log_ksi_log_t_acceptance_left_;
    TH2F log_ksi_log_t_acceptance_right_;
    TH2F log_ksi_log_t_total_acceptance_left_;
    TH2F log_ksi_log_t_total_acceptance_right_;
    TH1F log_dpe_mass_total_accepted_dist_;
    TH1F log_dpe_mass_reconstructed_dist_;
    TH1F log_dpe_mass_generated_dist_;
    TH1F log_dpe_mass_total_acceptance_;
    TH1F log_dpe_mass_reconstructed_acceptance_;
    TH1F log_dpe_mass_reconstruction_failures_;
    TH1F log_dpe_mass_reconstruction_failure_probability_;
    
    TH1F elast_log_t_distribution_;
    TH1F elast_log_t_accepted_;
    TH1F elast_log_t_acceptance_;
    TH1F elast_log_t_total_accepted_;
    TH1F elast_log_t_total_acceptance_;
    
    
    TH1F elast_left_thx_rec_error_;
    TH1F elast_left_thy_rec_error_;
    TH1F elast_right_thx_rec_error_;
    TH1F elast_right_thy_rec_error_;
    TH1F elast_scattering_thx_rec_error_;
    TH1F elast_scattering_thy_rec_error_;
    TH1F elast_left_xi_error_;
    TH1F elast_right_xi_error_;
    TH1F elast_thx_left_right_difference_;
    TH1F elast_thy_left_right_difference_;
    
    TH1F elast_x0_rec_error_;
    TH1F elast_y0_rec_error_;
    TH1F elast_z0_rec_error_;
    
    TH1F elast_chi2ndf_dist_;
    TH1F elast_probe_dist_;

//    TH1F elast_RP220_near_x_left_right_diff_;
//    TH1F elast_RP220_near_y_left_right_diff_;
//    TH1F elast_RP147_near_x_left_right_diff_;
//    TH1F elast_RP147_near_y_left_right_diff_;
//    TProfile elast_RP147_left_x_thy_dist_;
//    TProfile elast_RP147_left_y_thy_dist_;
//    TProfile elast_RP220_left_x_thy_dist_;
//    TProfile elast_RP220_left_y_thy_dist_;
//    TProfile elast_RP147_right_y_thy_dist_;
//    TProfile elast_RP147_right_x_thy_dist_;
//    TProfile elast_RP220_right_y_thy_dist_;
//    TProfile elast_RP220_right_x_thy_dist_;
    
    ReconstructionProfile ksi_error_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile smeared_ksi_error_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile t_error_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile phi_error_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile vx_error_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile vy_error_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile vz_error_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile thetax_error_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile thetay_error_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile smeared_thetax_error_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile smeared_thetay_error_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile chi2_error_overn_N_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile prob_function_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile reconstruction_failure_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile ksi_rel_error_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile t_rel_error_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile t_x_error_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile t_y_error_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile t_x_rel_error_vs_log_t_log_ksi_phi_left_;
    ReconstructionProfile t_y_rel_error_vs_log_t_log_ksi_phi_left_;
    
    ReconstructionProfile t_x_error_vs_log_tx_log_ksi_phi_left_;
    ReconstructionProfile t_y_error_vs_log_ty_log_ksi_phi_left_;
    ReconstructionProfile t_x_rel_error_vs_log_tx_log_ksi_phi_left_;
    ReconstructionProfile t_y_rel_error_vs_log_ty_log_ksi_phi_left_;
    
    ReconstructionProfile ksi_error_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile smeared_ksi_error_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile t_error_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile phi_error_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile vx_error_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile vy_error_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile vz_error_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile thetax_error_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile thetay_error_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile smeared_thetax_error_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile smeared_thetay_error_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile chi2_error_overn_N_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile prob_function_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile reconstruction_failure_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile ksi_rel_error_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile t_rel_error_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile t_x_error_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile t_y_error_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile t_x_rel_error_vs_log_t_log_ksi_phi_right_;
    ReconstructionProfile t_y_rel_error_vs_log_t_log_ksi_phi_right_;

    ReconstructionProfile t_x_error_vs_log_tx_log_ksi_phi_right_;
    ReconstructionProfile t_y_error_vs_log_ty_log_ksi_phi_right_;
    ReconstructionProfile t_x_rel_error_vs_log_tx_log_ksi_phi_right_;
    ReconstructionProfile t_y_rel_error_vs_log_ty_log_ksi_phi_right_;
    
    //lower means abs(ksi) lower
    ReconstructionProfile M_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile M_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile ksi_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile ksi_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile ksi_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile t_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile t_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile t_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile ksi_higher_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile ksi_lower_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile ksi_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile t_higher_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile t_lower_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile t_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile phi_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile phi_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile phi_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    
    ReconstructionProfile x_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile y_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile z_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    
    ReconstructionProfile probe_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
    ReconstructionProfile chsqndf_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_;
};

#endif

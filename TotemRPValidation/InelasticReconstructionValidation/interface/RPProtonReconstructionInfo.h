#ifndef TotemRPValidation_RPReconstructedTracksValidation_RPProtonReconstructionInfo_h
#define TotemRPValidation_RPReconstructedTracksValidation_RPProtonReconstructionInfo_h

#include "TotemRPValidation/BaseValidationClasses/interface/BaseHistogramManager.h"
#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile3D.h"
#include "TotemRPValidation/ValidationTools/interface/ReconstructionProfile.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProton.h"
#include "HepMC/GenEvent.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class RPProtonReconstructionInfo : public BaseHistogramManager
{
  public:
    RPProtonReconstructionInfo(const std::string &path, int arm_id, const edm::ParameterSet& conf);
    
    void FillReferenceHistograms(const BeamOpticsParams &BOPar, 
        const HepMC::GenParticle *pr_prot, const HepMC::GenParticle *pr_vert, 
        bool verbosity=false);
    void FillResidualHistograms(const BeamOpticsParams &BOPar, 
        const HepMC::GenParticle *pr_prot, const HepMC::GenParticle *pr_vert, 
        const RPReconstructedProton* rec_prot, 
        bool verbosity=false);
    
  private:
    void InitializeHistograms();
    virtual void Finalize();
    double ReduceAngleTo_0_2Pi(double angle);
    double ReduceAngleTo_Pi_Pi(double angle);
    double AngleDiff(double angle1, double angle2);
    
    double chi2_upper_tail_prob_cut_;
    
    int arm_id_;
    
    TH2F log_ksi_log_t_distribution_;
    TH2F log_ksi_log_t_reconstructed_protons_;
    TH2F log_ksi_log_t_reconstruction_failed_;
    TH2F log_ksi_log_t_reconstruction_failure_rate_;
    TH2F log_ksi_log_t_total_accepted_;
    TH2F log_ksi_log_t_reconstruction_failure_to_total_accepted_;
    TH2F log_ksi_log_t_acceptance_;
    TH2F log_ksi_log_t_total_acceptance_;

    TH1F log_ksi_distribution_;
    TH1F log_ksi_reconstructed_protons_;
    TH1F log_ksi_acceptance_;
    TH1F log_ksi_total_accepted_;
    TH1F log_ksi_total_acceptance_;
    
    ReconstructionProfile ksi_error_vs_log_t_log_ksi_phi_;
    ReconstructionProfile smeared_ksi_error_vs_log_t_log_ksi_phi_;
    ReconstructionProfile t_error_vs_log_t_log_ksi_phi_;
    ReconstructionProfile phi_error_vs_log_t_log_ksi_phi_;
    ReconstructionProfile vx_error_vs_log_t_log_ksi_phi_;
    ReconstructionProfile vy_error_vs_log_t_log_ksi_phi_;
    ReconstructionProfile thetax_error_vs_log_t_log_ksi_phi_;
    ReconstructionProfile smeared_thetax_error_vs_log_t_log_ksi_phi_;
    ReconstructionProfile thetay_error_vs_log_t_log_ksi_phi_;
    ReconstructionProfile smeared_thetay_error_vs_log_t_log_ksi_phi_;
    ReconstructionProfile chi2_error_overn_N_vs_log_t_log_ksi_phi_;
    
    ReconstructionProfile prob_function_vs_log_t_log_ksi_phi_;
    
    ReconstructionProfile reconstruction_failure_vs_log_t_log_ksi_phi_;
    
    ReconstructionProfile ksi_rel_error_vs_log_t_log_ksi_phi_;
    ReconstructionProfile t_rel_error_vs_log_t_log_ksi_phi_;
    ReconstructionProfile t_x_error_vs_log_t_log_ksi_phi_;
    ReconstructionProfile t_x_error_vs_log_tx_log_ksi_phi_;
    ReconstructionProfile t_y_error_vs_log_t_log_ksi_phi_;
    ReconstructionProfile t_y_error_vs_log_ty_log_ksi_phi_;
    ReconstructionProfile t_x_rel_error_vs_log_t_log_ksi_phi_;
    ReconstructionProfile t_x_rel_error_vs_log_tx_log_ksi_phi_;
    ReconstructionProfile t_y_rel_error_vs_log_t_log_ksi_phi_;
    ReconstructionProfile t_y_rel_error_vs_log_ty_log_ksi_phi_;
};

#endif

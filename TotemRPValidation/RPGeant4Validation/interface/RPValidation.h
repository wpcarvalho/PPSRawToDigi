#ifndef TotemRPValidation_RPGeant4Validation_RPValidation_h
#define TotemRPValidation_RPGeant4Validation_RPValidation_h

#include "TotemRPValidation/BaseValidationClasses/interface/BaseHistogramManager.h"
#include "TotemRPValidation/RPGeant4Validation/interface/RPSupplementaryInfo.h"
#include "Geometry/TotemRPDetTopology/interface/RPHepPDTWrapper.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class RPValidation : public BaseHistogramManager
{
  public:
    RPValidation(const std::string &path, RPId rp_id, const edm::ParameterSet&);
    void FillHistogramms(RPSupplementaryInfo &supp_info);
    void FillHistogramms(const PSimHit & sim_hit, RPSupplementaryInfo &supp_info);
  private:
    void InitializeHistograms();
    virtual void Finalize();
    
    RPId _rp_id;
    
    TH1F _pot_interaction_parts;
    TH1F _pot_scat_dist_x;
    TH1F _pot_scat_dist_y;
    TH1F _en_dep_for_primary_protons;
    TH1F _angle_at_RP_x;
    TH1F _angle_at_RP_y;
    TH2F _position_at_RP_xy;
    TH1F _interactions_per_rp_part;
    TH1F _front_wall_inter_track_no_dist;
    TH1F _front_wall_inter_PDG_dist;
    TH2F _front_wall_momentum_theta_dist;
    TH2F _front_wall_momentum_theta_charged_particle_dist;
    
    TH1F _front_wall_inter_track_no_dist_charged_above_1_GeV;
    TH1F _front_wall_inter_PDG_dist_charged_above_1_GeV;
    TH1F _front_wall_inter_theta_dist_charged_above_1_GeV;
    
    double mp;
    double p0;
    double E1;
    TH1F _log10_t_dist;
    TH1F _log10_t_dist_inelastic;
    TH1F _log10_t_dist_inelastic_ratio;
    TH1F _log10_t_dist_entered_rp;
    TH1F _log10_t_dist_geometrical_acceptance;
};

#endif

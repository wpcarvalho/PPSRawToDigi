#ifndef TotemRPValidation_RPGeant4Validation_StationValidation_h
#define TotemRPValidation_RPGeant4Validation_StationValidation_h

#include "TotemRPValidation/BaseValidationClasses/interface/BaseHistogramManager.h"
#include "TotemRPValidation/RPGeant4Validation/interface/RPSupplementaryInfo.h"
#include "Geometry/TotemRPDetTopology/interface/RPHepPDTWrapper.h"
#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class StationValidation : public BaseHistogramManager
{
  public:
    StationValidation(const std::string &path, RPStationId st_id, const edm::ParameterSet&);
    void FillHistogramms(RPSupplementaryInfo &supp_info);
    virtual void Finalize();
  private:
    void Init();
    
    RPStationId _station_id;
    
    TH1F _prim_prot_at_station_mom_dist;
    TH1F _prim_prot_at_station_theta_dist_x;
    TH1F _prim_prot_at_station_theta_dist_y;
    TH1F _prim_prot_at_station_pos_dist_x;
    TH1F _prim_prot_at_station_pos_dist_y;
    TH2F _prim_prot_at_station_pos_dist_xy;
    TH2F _prim_prot_pos_at_station_vs_theta_at_0_x;
    TH2F _prim_prot_pos_at_station_vs_theta_at_0_y;
    TH2F _prim_prot_pos_at_station_vs_pos_at_0_x;
    TH2F _prim_prot_pos_at_station_vs_pos_at_0_y;
    TH2F _prim_prot_pos_at_station_vs_ksi_x;
    TH2F _prim_prot_pos_at_station_vs_ksi_y;
    TH2F _prim_prot_theta_at_station_vs_theta_at_0_x;
    TH2F _prim_prot_theta_at_station_vs_theta_at_0_y;

    double mp, p0, E1;
};


#endif

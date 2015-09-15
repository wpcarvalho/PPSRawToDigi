#ifndef TotemRPValidation_RPReconstructedTracksValidation_RPTrackInfo_h
#define TotemRPValidation_RPReconstructedTracksValidation_RPTrackInfo_h

#include "TotemRPValidation/BaseValidationClasses/interface/BaseHistogramManager.h"
#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrack.h"
#include "DataFormats/TotemL1Trigger/interface/RPCCBits.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class RPTrackInfo : public BaseHistogramManager
{
  public:
    RPTrackInfo(const std::string &path, RPId rp_id, const edm::ParameterSet& conf);
    
    void Fill_t_Dist(double t);
    void FillTrackSeen(double t);
    void Fill_ksi_t_Dist(double t, double ksi);
    void FillTrackSeen(double t, double ksi);
    void FillRPPositionResiduals(double delta_x, double delta_y);
    void FillRPPositionResiduals2ProtonsSeenAt150Station(double delta_x, double delta_y);
    void FillPools(double delta_x, double delta_y);
    void FillBasicTrackInfo(const RPFittedTrack& track);
    void FillRPCoincidenceChipInfo(const RPCCBits& u_coincidence_bits, const RPCCBits& v_coincidence_bits);
    
  private:
    void InitializeHistograms();
    virtual void Finalize();
    
    RPId rp_id_;
    
//    double mp_;
//    double p0_;
//    double E1_;
    TH1F log10_t_dist_;
    TH1F log10_t_dist_seen_tracks_;
    TH1F log10_t_dist_acceptance_;

    TH1F log10_xi_dist_;
    TH1F log10_xi_dist_seen_tracks_;
    TH1F log10_xi_dist_acceptance_;
    
    TH2F log_ksi_log_t_distribution_;
    TH2F log_ksi_log_t_seen_tracks_;
    TH2F log_ksi_log_t_acceptance_;
    
    TH1F position_residuum_x_;
    TH1F position_residuum_y_;
    
    TH1F position_residuum_x_if_exactly_2_rps_in_150_station_;
    TH1F position_residuum_y_if_exactly_2_rps_in_150_station_;
    
    TH1F position_pool_x_;
    TH1F position_pool_y_;
    
    TH2F track_distribution_xy_;
    TH1F track_distribution_x_;
    TH1F track_distribution_y_;
    TH1F track_distribution_thx_;
    TH1F track_distribution_thy_;
    TH1F track_reco_hit_no_;
    TH1F track_chi_sq_ndf_;
    
    TH1F u_sectors_on_;
    TH1F v_sectors_on_;
    TH2F uv_sectors_on_;
};

#endif

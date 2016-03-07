#ifndef TotemRPValidation_RPReconstructedTracksValidation_RPStationInfo_h
#define TotemRPValidation_RPReconstructedTracksValidation_RPStationInfo_h

#include "TotemRPValidation/BaseValidationClasses/interface/BaseHistogramManager.h"
#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
//#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrack.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RP2DHit.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class RPStationInfo : public BaseHistogramManager
{
  public:
    RPStationInfo(const std::string &path, RPStationId rp_station_id, const edm::ParameterSet& conf);
    
    void Fill_t_Dist(double t);
    void FillTrackSeen(double t);
    void Fill_ksi_t_Dist(double t, double ksi);
    void FillTrackSeen(double t, double ksi);
    void FillBasicTrackInfo(const RP2DHit& track_in, const RP2DHit& track_out, double station_z_centre);
    
  private:
    void InitializeHistograms();
    virtual void Finalize();
    
    RPStationId rp_station_id_;
    
//    double mp_;
//    double p0_;
//    double E1_;
    TH1F log10_t_dist_;
    TH1F log10_t_dist_seen_tracks_;
    TH1F log10_t_dist_acceptance_;
    
    TH2F log_ksi_log_t_distribution_;
    TH2F log_ksi_log_t_seen_tracks_;
    TH2F log_ksi_log_t_acceptance_;
        
    TH2F track_distribution_xy_centre_;
    
    TH1F track_distribution_x_centre_;
    TH1F track_distribution_y_centre_;
    
    TH1F track_distribution_thx_centre_;
    TH1F track_distribution_thy_centre_;
    
    TH2F track_dist_thx_vs_x_centre_;
    TH2F track_dist_thy_vs_y_centre_;
    TH2F track_dist_thy_vs_thx_centre_;
};

#endif

#ifndef TotemRPValidation_RPReconstructedTracksValidation_DetInfo_h
#define TotemRPValidation_RPReconstructedTracksValidation_DetInfo_h

#include "TotemRPValidation/BaseValidationClasses/interface/BaseHistogramManager.h"
#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDigCluster.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class DetInfo : public BaseHistogramManager
{
  public:
    DetInfo(const std::string &path, unsigned int det_id, const edm::ParameterSet& conf);
    
    typedef edm::DetSet<RPDigCluster> clu_col_type;
    
    void FillHistograms(const clu_col_type &col, double residual);
    void FillTriggerSector(int trigger_sector);
    
  private:
    void InitializeHistograms();
    
    RPDetId det_id_;
    
    TH1F cluster_size_;
    TH1F cluster_multiplicity_;
    TH1F strip_profile_;
    TH1F track_residual_;
    TH1F trigger_sectors_;
};

#endif

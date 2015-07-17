#ifndef TotemRPValidation_InelasticReconstructionValidation_PoolsReconstructionInfo_h
#define TotemRPValidation_InelasticReconstructionValidation_PoolsReconstructionInfo_h

#include "TotemRPValidation/BaseValidationClasses/interface/BaseHistogramManager.h"
#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RP2DHit.h"
#include "TH1F.h"
#include "TVector2.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


class PoolsReconstructionInfo : public BaseHistogramManager
{
  public:
    PoolsReconstructionInfo(const std::string &path, RPId rp_id, const edm::ParameterSet&);
    void Fill(const TVector2 &theoretical_hit, const RP2DHit &simulated_hit, 
          const TVector2 &reconstructed_hit);
          
  private:
    void InitializeHistograms();
    RPId rp_id_;
    
    //theoretical: primary proton with crossing angle transformation transported to the RP
    //simulated: smeared primary proton transported to the RP
    //reconstructed: reconstructed proton transported to the RP
    TH1F simulated_vs_theoretical_residuum_x_;
    TH1F simulated_vs_reconstructed_residuum_x_;
    TH1F reconstructed_vs_theoretical_residuum_x_;
    
    TH1F simulated_vs_theoretical_residuum_y_;
    TH1F simulated_vs_reconstructed_residuum_y_;
    TH1F reconstructed_vs_theoretical_residuum_y_;
    
    //TH1F reconstructed_vs_theoretical_pool_x_;
    //TH1F reconstructed_vs_theoretical_pool_y_;
};


#endif


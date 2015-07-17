#ifndef TotemRPValidation_RPGeant4Validation_RPGeant4Validation_h
#define TotemRPValidation_RPGeant4Validation_RPGeant4Validation_h

#include "TotemRPValidation/RPGeant4Validation/interface/RPValidation.h"
#include "TotemRPValidation/RPGeant4Validation/interface/RPSupplementaryInfo.h"
#include "TotemRPValidation/BaseValidationClasses/interface/BaseCollectionManager.h"
#include "TotemRPValidation/RPGeant4Validation/interface/RPPSimHitDebugInfo.h"
#include "TotemRPValidation/RPGeant4Validation/interface/StationValidation.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <vector>
#include <string>

class RPGeant4Validation
{
  public:
    RPGeant4Validation();
    void FillHistogramms(const std::vector<PSimHit>& hits, const std::vector<RPPSimHitDebugInfo>& debug_hits);
    void WriteHistograms(const std::string &root_file_nam);
  private:
    BaseCollectionManager<RPValidation, unsigned int, edm::ParameterSet> rp_coll;
    BaseCollectionManager<StationValidation, unsigned int, edm::ParameterSet> station_coll;
    
    std::vector<RPId> rp_ids;
    std::vector<RPStationId> station_ids;
    RPDetSpaceGeometry det_sp_geom;
};

#endif

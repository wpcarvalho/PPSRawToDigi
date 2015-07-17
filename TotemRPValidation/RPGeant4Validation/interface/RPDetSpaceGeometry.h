#ifndef TotemRPValidation_RPGeant4Validation_RP_DET_SPACE_GEOMETRY_H
#define TotemRPValidation_RPGeant4Validation_RP_DET_SPACE_GEOMETRY_H

#include <vector>
#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"


class RPDetSpaceGeometry
{
  public:
    RPDetSpaceGeometry();
    std::vector<RPDetId> GetRPDetIdList();
    std::vector<RPDetId> GetRPDetIdList(RPId rp_id);
    std::vector<RPStationId> GetRPStationIdList();
    std::vector<RPId> GetRPRomanPotsIdList(RPStationId station_id);
    std::vector<RPId> GetAllRPRomanPotsIdList();
    static inline RPStationId GetStationId(RPDetId rp_det_id) {return (rp_det_id/100)*100;}
    inline RPId GetRPId(RPDetId rp_det_id) {return (rp_det_id/10)*10;}
    inline int GetDetSeqNoInRP(RPDetId rp_det_id) {return rp_det_id%10;}  //returns number from 0 to 9
    static inline bool IsOnTheRight(RPStationId station_id) {return station_id/1000;}
//    static inline bool IsOnTheRight(RPDetId rp_det_id) {return rp_det_id/1000;}
//    static inline bool IsOnTheRight(RPId rp_id) {return rp_id/1000;}
    static bool IsInsideStation(double z);
    static RPStationId GetStationIdByVertex(double z);
    static RPStationId GetStationIdBySeqNo(int seq_no);
  private:
}
;

#endif

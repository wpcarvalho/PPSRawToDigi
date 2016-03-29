#ifndef RecoTotemRP_RPRecoDataFormats_RPTrackCandidateCollection_h
#define RecoTotemRP_RPRecoDataFormats_RPTrackCandidateCollection_h
 
#include "RecoTotemRP/RPRecoDataFormats/interface/RPTrackCandidate.h"
#include "DataFormats/TotemRPDetId/interface/TotemRPIdTypes.h"

#include <map>
 
class RPTrackCandidateCollection : public std::map<unsigned int, RPTrackCandidate>
{
  public:
    //inline unsigned int RomanPotId() {return rp_id_;}
    //inline void RomanPotId(unsigned int rp_id) {rp_id_ = rp_id;}
  private:
    //unsigned int rp_id_;
};

#endif

#ifndef RecoTotemRP_RPRecoDataFormats_RPTrackCandidateDistinctCollectionsSet_h
#define RecoTotemRP_RPRecoDataFormats_RPTrackCandidateDistinctCollectionsSet_h
 
#include "RecoTotemRP/RPRecoDataFormats/interface/RPTrackCandidateCollection.h"
#include <vector>
#include <map>


class RPTrackCandidateDistinctCollectionsSet
   : public std::map<unsigned int, std::vector<RPTrackCandidateCollection> >
{
  public:
    //inline unsigned int RomanPotId() {return rp_id_;}
    //inline void RomanPotId(unsigned int rp_id) {rp_id_ = rp_id;}
  private:
    //unsigned int rp_id_;
};


#endif

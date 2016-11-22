#include "DataFormats/CTPPSReco/interface/DiamondRecHit.h"


DiamondRecHit::DiamondRecHit() : id_(0), time_(0), time_error_(0), 
                                x_(0), x_error_(0), y_(0), y_error_(0), flags_(0) {
}

DiamondRecHit::DiamondRecHit(const DetId& id, float time, float time_error, 
                        float x, float x_error, float y, float y_error, uint32_t flags) : 
	id_(id) , time_(time) , time_error_(time_error) , 
	x_(x) , x_error_(x_error) , y_(y) , y_error_(y_error) , flags_(flags) {
}

bool operator< (const DiamondRecHit& h1, const DiamondRecHit& h2) 
{
   return h1.detId() < h2.detId() ;
}

std::ostream& operator<<(std::ostream& s, const DiamondRecHit& hit) {
  s << hit.detId().rawId() << " , " << hit.time() << " +/- " << hit.timeError() << " ns , ";
  s << hit.x() << " +/- " << hit.xError() << " mm , ";
  s << hit.y() << " +/- " << hit.yError() << " mm , ";
  s << " flags=0x" << std::hex << hit.flags() << std::dec << " " ;
  return s;
}


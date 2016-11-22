#ifndef DATAFORMATS_DIAMONDRECHIT_H
#define DATAFORMATS_DIAMONDRECHIT_H 1

#include "DataFormats/DetId/interface/DetId.h"
#include <ostream>


/*
 * \class DiamondRecHit
 * 
 * \author W. Carvalho (UERJ)
 *
*/

class DiamondRecHit {

public:

  typedef DetId key_type;
  
  DiamondRecHit(); // for persistence
  explicit DiamondRecHit(const DetId& id, float time, float time_error, 
                float x, float x_error, float y, float y_error, uint32_t flags = 0);

  const DetId& detId() const { return id_; }

  float time() const { return time_ ; }
  void setTime(float time) { time_ = time ; }

  float timeError() const { return time_error_ ; }
  void setTimeError(float time_error) { time_error_ = time_error ; }
  
  uint32_t flags() const { return flags_ ; }
  void setFlags(uint32_t flags) { flags_ = flags ; }

  float x() const { return x_ ; }
  void setX(float x) { x_ = x ; }

  float xError() const { return x_error_ ; }
  void setXError(float x_error) { x_error_ = x_error ; }

  float y() const { return y_ ; }
  void setY(float y) { y_ = y ; }

  float yError() const { return y_error_ ; }
  void setYError(float y_error) { y_error_ = y_error ; }
  
private:

  DetId id_;
  float time_;
  float time_error_;
  // position and uncertainty of the hit in mm, wrt detector center 
  // (likewise TotemRPRecHit.h)
  float x_, x_error_;   
  float y_, y_error_;      
  uint32_t flags_;

};

bool operator<(const DiamondRecHit& , const DiamondRecHit&);

std::ostream& operator<<(std::ostream& s, const DiamondRecHit& hit);
  
#endif

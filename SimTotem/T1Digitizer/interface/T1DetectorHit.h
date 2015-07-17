#ifndef T1__DETECTOR_HIT_H
#define T1__DETECTOR_HIT_H

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/T1DetId/interface/T1DetId.h"
#include <iosfwd>

//class PSimHit;
//class T1DetId;
class T1DetectorHit
{
 public:
  T1DetectorHit(int element, float charge, float position, float time,
		const PSimHit * hitp = 0, uint32_t id=0)
    : theElement(element), theCharge(charge),
    thePosition(position),   theTime(time), theHitp(hitp), theDetId(id) {}

  int   getElement()  const {return theElement;}
  float getCharge()   const {return theCharge;}
  float getPosition() const {return thePosition;}
  float getTime()     const {return theTime;}
  const PSimHit * getSimHit() const {return theHitp;}
  uint32_t getDetId()  const {return theDetId;}
  friend std::ostream & operator<<(std::ostream &, const T1DetectorHit &);
 private:
  /// strip or wire number
  int   theElement;
  float theCharge;
  /// the position is along the element, with (0,0) the center of the chamber
  float thePosition; 
  /// start counting time at the beam crossing
  float theTime;
  /// theSimHit that created this hit
  const PSimHit * theHitp;

  uint32_t theDetId;
};

#endif


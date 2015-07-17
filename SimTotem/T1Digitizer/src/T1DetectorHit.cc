#include "SimTotem/T1Digitizer/interface/T1DetectorHit.h"
#include <iostream>

std::ostream & operator<<(std::ostream & stream, const T1DetectorHit & hit) {
  stream << "element: " << hit.getElement()
         << "  charge: " << hit.getCharge()
         << "   pos:  " << hit.getPosition()
         << "   time: " << hit.getTime()   
	 <<   "  DetId "  << hit.getDetId() << std::endl;
  return stream;
}


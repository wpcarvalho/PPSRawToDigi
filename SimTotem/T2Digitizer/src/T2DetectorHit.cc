#include "SimTotem/T2Digitizer/interface/T2DetectorHit.h"
#include <iostream>

std::ostream & operator<<(std::ostream & stream, const T2DetectorHit & hit) {
  stream << "element: " << hit.getElement()
         << "  charge: " << hit.getCharge()
         << "  row: " << hit.getRow()
         << "  column: " << hit.getCol()
         << "  time: " << hit.getTime() << std::endl; 
  return stream;
}

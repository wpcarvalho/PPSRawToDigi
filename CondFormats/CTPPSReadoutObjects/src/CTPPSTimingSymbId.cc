/****************************************************************************
 *   Seyd Mohsen Etesami
 ****************************************************************************/

#include "CondFormats/CTPPSReadoutObjects/interface/CTPPSTimingSymbId.h"


std::ostream& operator << (std::ostream& s, const CTPPSTimingSymbID &sid)
{ 
  switch (sid.subSystem) {
  case CTPPSTimingSymbID::Diamond: s << "sub-system=Diamond, "; break;
  case CTPPSTimingSymbID::FastSi: s << "sub-system=FastSi, "; break;
  case CTPPSTimingSymbID::Quartic: s << "sub-system=Quartic, "; break;


  } 
  s << "symb. id=" << sid.symbolicID;
  return s;
}


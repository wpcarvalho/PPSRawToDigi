/****************************************************************************
*   Seyed Mohsen Etesami
****************************************************************************/

#ifndef CondFormats_CTPPSReadoutObjects_CTPPSTimingSymbId
#define CondFormats_CTPPSReadoutObjects_CTPPSTimingSymbId

#include <iostream>

/**
 *\brief Symbolic ID describing an entity of a CTPPS Timing subdetector.
 **/
class CTPPSTimingSymbID
{
  public:
    /// identifies the CTPPS Timing  subsystem
    enum {Diamond, FastSi, Quartic} subSystem;

    /// integer-encoded symbolic ID
    unsigned int symbolicID;

    bool operator < (const CTPPSTimingSymbID &sid) const
    {
      if (subSystem == sid.subSystem)
		  return (symbolicID < sid.symbolicID);
      return (subSystem < sid.subSystem);
    }

    bool operator == (const CTPPSTimingSymbID &sid) const
    {
      return ((subSystem==sid.subSystem) && (symbolicID==sid.symbolicID));
    }
    
    friend std::ostream& operator << (std::ostream& s, const CTPPSTimingSymbID &sid);
};

#endif

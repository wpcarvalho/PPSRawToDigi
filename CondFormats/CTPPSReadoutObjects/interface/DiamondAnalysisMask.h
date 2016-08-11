#ifndef CondFormats_CTPPSReadoutObjects_DiamondAnalysisMask
#define CondFormats_CTPPSReadoutObjects_DiamondAnalysisMask

#include "CondFormats/CTPPSReadoutObjects/interface/CTPPSTimingSymbId.h"

#include <set>
#include <map>

//----------------------------------------------------------------------------------------------------

/**
 *\brief Contains data on masked channels of a VFAT.
 */
class DiamondVFATAnalysisMask
{
 public:
 DiamondVFATAnalysisMask() : fullMask(false) {}

  /// whether all channels of the VFAT shall be masked
  bool fullMask;

  /// list of channels to be masked
  std::set<unsigned char> maskedChannels;
};

//----------------------------------------------------------------------------------------------------

/**
 *\brief Channel-mask mapping.
 **/
class DiamondAnalysisMask
{
 public:
  std::map<CTPPSTimingSymbID, DiamondVFATAnalysisMask> analysisMask;

  void insert(const CTPPSTimingSymbID &sid, const DiamondVFATAnalysisMask &vam);
};

#endif

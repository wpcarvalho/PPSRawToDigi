/****************************************************************************
 *   Seyed Mohsen Etesami
 ****************************************************************************/

#ifndef CondFormats_CTPPSReadoutObjects_TotemDAQMappingDiamond
#define CondFormats_CTPPSReadoutObjects_TotemDAQMappingDiamond

#include "CondFormats/CTPPSReadoutObjects/interface/DiamondFramePosition.h"

#include "CondFormats/CTPPSReadoutObjects/interface/CTPPSTimingSymbId.h"

#include <map>
//----------------------------------------------------------------------------------------------------

/**
 *\brief Contains mappind data related to a VFAT.
 */
class DiamondVFATInfo
{
 public:
  /// is data of coincidence-chip VFAT
  enum {data, CC} type;

  /// the symbolic id
  CTPPSTimingSymbID symbolicID;

  /// the hardware ID (16 bit)
  unsigned int hwID;
    
  friend std::ostream& operator << (std::ostream& s, const DiamondVFATInfo &fp);
};

//----------------------------------------------------------------------------------------------------

/**
 *\brief The mapping between FramePosition and VFATInfo.
 */
class TotemDAQMappingDiamond
{
 public:
  std::map<DiamondFramePosition, DiamondVFATInfo> VFATMapping;
    
  void insert(const DiamondFramePosition &fp, const DiamondVFATInfo &vi);
};

#endif

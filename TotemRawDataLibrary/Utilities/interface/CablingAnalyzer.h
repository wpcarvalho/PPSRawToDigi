/****************************************************************************
*
* This is a part of the TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/


#ifndef _Totem_CablingAnalyzer_h_
#define _Totem_CablingAnalyzer_h_

#include "TotemRawDataLibrary/DataFormats/interface/FramePosition.h"

namespace Totem {
  class DataFile;
}

class VFAT2RegistersXML;

/**
\ingroup TotemRawDataLibrary
\brief Creates map (VFAT ID --> Slink position) from a data file.
**/

class CablingAnalyzer
{
  public:
    struct record
    {
      unsigned short ID;
      Totem::FramePosition index;
      record(unsigned short _ID = 0, Totem::FramePosition i = 0) : ID(_ID), index(i) {}
    };

#ifdef _MONITOR_
    unsigned short run(Totem::DataFile *, VFAT2RegistersXML *) const;
#endif
};

inline bool operator == (const CablingAnalyzer::record &a, const CablingAnalyzer::record &b)
{
  return (a.ID == b.ID && a.index == b.index) ? true : false ;
}

#endif

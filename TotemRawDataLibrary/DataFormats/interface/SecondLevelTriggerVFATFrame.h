/**********************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* This is a part of TOTEM offline software.
* Authors:
*   Michal Zmuda (m.zmuda@aol.com)
*
**********************************************************/

#pragma once

#include <vector>
#include "TotemRawDataLibrary/DataFormats/interface/OptoRxSupplementalData.h"
#include "TotemRawDataLibrary/DataFormats/interface/VFATFrame.h"


namespace Totem {

/**
 * \ingroup TotemRawDataLibrary
 * TODO: describe.
 **/
class SecondLevelTriggerVFATFrame
{
  private:
    unsigned short header;
    unsigned short tailer;
    unsigned dataSize;
    unsigned short* data;

  public:
    SecondLevelTriggerVFATFrame(unsigned short* data);
    ~SecondLevelTriggerVFATFrame();

    unsigned short getSize() const;
    bool checkFootprint() const;
    std::vector< std::pair<unsigned char, unsigned char> >* getData() const;
    void Print() const;
    void PrintBinary() const;
};

}

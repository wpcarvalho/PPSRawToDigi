/**********************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* This is a part of TOTEM offline software.
* Authors:
*   Michal Zmuda (m.zmuda@aol.com)
*
**********************************************************/

#ifndef _Totem_RawVFATFrame_h_
#define _Totem_RawVFATFrame_h_

#include <vector>
#include <bitset>
#include "TotemRawDataLibrary/DataFormats/interface/PositionedVFATFrame.h"


namespace Totem {

/**
 * \ingroup TotemRawDataLibrary
 * TODO: describe.
 **/
class RawVFATFrame: public PositionedVFATFrame
{
private:
    word infoHolder[5];
    std::vector<unsigned char> activeChannels;

public:
    // TODO: required by CMSSW, define behaviour !
    RawVFATFrame() {}

    RawVFATFrame(word data[13]);
    ~RawVFATFrame();

    word getRxFlags() const;
    word getFiberIdx() const;
    word getGohIdx() const;
    word getSize() const;
    bool passedHardwareSynchronisationChecks() const;

    word getBC() const;
    word getEC() const;
    word getFlags() const;
    word getChipID() const;

    bool checkFootprint() const;
    bool checkCRC() const;
    bool channelActive(unsigned char channel) const;
    std::vector<unsigned char> getActiveChannels() const;
    void Print(bool binary = false) const;
};
#endif
}

/**********************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* This is a part of TOTEM offline software.
* Authors:
*   Michal Zmuda (m.zmuda@aol.com)
*
**********************************************************/

#ifndef _Totem_ClusterizationVFATFrame_h_
#define _Totem_ClusterizationVFATFrame_h_

#include <vector>
#include <bitset>
#include "TotemRawDataLibrary/DataFormats/interface/PositionedVFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/OptoRxSupplementalData.h"


namespace Totem {

/**
 * \ingroup TotemRawDataLibrary
 * TODO: describe.
 **/
class ClusterizationVFATFrame: public PositionedVFATFrame
{
private:
    word infoHolder[5];
    word* dataHolder;
    unsigned dataSize;
    std::pair<unsigned char, unsigned char> getClustersAt(unsigned i) const;

public:
    // TODO: required by CMSSW, define behaviour !
    ClusterizationVFATFrame() {}

    ClusterizationVFATFrame(word* data);
    ~ClusterizationVFATFrame();

    word getRxFlags() const;
    word getFiberIdx() const;
    word getGohIdx() const;
    word getSize() const;
    bool passedHardwareSynchronisationChecks() const;

    word hasBC() const;
    word hasEC() const;
    word hasChipID() const;

    word getBC() const;
    word getEC() const;
    word getChipID() const;
    word getFlags() const;

    bool checkFootprint() const;
    bool checkCRC() const;
    bool channelActive(unsigned char channel) const;
    std::vector<unsigned char> getActiveChannels() const;
    std::vector< std::pair<unsigned char, unsigned char> > getActiveClusters() const;
    void Print(bool binary = false) const;
};
#endif
}

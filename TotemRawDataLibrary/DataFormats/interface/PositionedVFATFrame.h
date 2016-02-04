/**********************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* This is a part of TOTEM offline software.
* Authors:
*   Michal Zmuda (m.zmuda@aol.com)
*
**********************************************************/

#ifndef _Totem_PositionedVFATFrame_h_
#define _Totem_PositionedVFATFrame_h_

#include "TotemRawDataLibrary/DataFormats/interface/VFATFrame.h"


namespace Totem {

/**
 * \ingroup TotemRawDataLibrary
 * TODO: describe.
 **/
class PositionedVFATFrame: public VFATFrame
{
  public:
    virtual ~PositionedVFATFrame(){}

    virtual word getRxFlags() const = 0; // todo replace with dedicated getters for each flags when their definitions are known
    virtual word getFiberIdx() const = 0;
    virtual word getGohIdx() const = 0;
    virtual word getSize() const = 0;

    virtual bool passedHardwareSynchronisationChecks() const = 0;
};

}
#endif

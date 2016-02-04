/**********************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* This is a part of TOTEM offline software.
* Authors:
*   Michal Zmuda (m.zmuda@aol.com)
*
**********************************************************/


#ifndef _Totem_PositionedVFATFrameCollection_h_
#define _Totem_PositionedVFATFrameCollection_h_

#include "TotemRawDataLibrary/DataFormats/interface/VFATFrameCollection.h"
#include "TotemRawDataLibrary/DataFormats/interface/PositionedVFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/SecondLevelTriggerVFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/FramePosition.h"
#include <map>


namespace Totem {

/**
 * \ingroup TotemRawDataLibrary
 * TODO: describe.
 **/
class PositionedVFATFrameCollection : public VFATFrameCollection
{
private:
    std::map<FramePosition, PositionedVFATFrame*> positionToFrame;
    std::map<unsigned, PositionedVFATFrame*> idToFrame;
    std::vector<SecondLevelTriggerVFATFrame*> triggerFrames;

public:
    ~PositionedVFATFrameCollection();
    std::string GetClassName() const { return "PositionedVFATFrameCollection"; }
    const PositionedVFATFrame* GetFrameByID(unsigned int ID) const;
    const PositionedVFATFrame* GetFrameByIndex(FramePosition index) const;
    const PositionedVFATFrame* GetFrameByIndexID(FramePosition index, unsigned int ID) const;
    unsigned int Size() const;
    bool Empty() const;

    virtual void AddNewFrame(FramePosition position, PositionedVFATFrame* frame);

    //todo - find better container/application than this one
    virtual void AddTriggerFrame(SecondLevelTriggerVFATFrame* frame);
    virtual std::vector<SecondLevelTriggerVFATFrame*>* GetTriggerFrames() const;
protected:
    value_type BeginIterator() const;
    value_type NextIterator(const value_type&) const;
    bool IsEndIterator(const value_type&) const;
};

} // namespace
#endif


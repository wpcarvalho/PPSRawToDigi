/**********************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors:
*   Michal Zmuda (m.zmuda@aol.com)
*
**********************************************************/

#include "TotemRawDataLibrary/DataFormats/interface/PositionedVFATFrameCollection.h"
#include "TotemRawDataLibrary/DataFormats/interface/CommonDef.h"
#include <map>
#include <cstdio>

using namespace std;

namespace Totem {


PositionedVFATFrameCollection::~PositionedVFATFrameCollection()
{
    map<FramePosition, PositionedVFATFrame*>::iterator i = positionToFrame.begin();
    for(; i != positionToFrame.end(); i++)
        delete i->second;

    vector<SecondLevelTriggerVFATFrame*>::iterator j = triggerFrames.begin();
    for(; j != triggerFrames.end(); j++){
//        delete (*j);
    }
}

void PositionedVFATFrameCollection::AddNewFrame(FramePosition position, PositionedVFATFrame* frame)
{
    if(positionToFrame[position] != NULL) WARN("PositionedVFATFrameCollection::AddNewFrame,Adding second frame with that position (" << position << ") - overriding.");
    if(idToFrame[frame->getChipID()] != NULL) WARN("PositionedVFATFrameCollection::AddNewFrame,Adding second frame with that ID (" << frame->getChipID() << ") - overriding.");

    positionToFrame[position] = frame;
    idToFrame[frame->getChipID()] = frame;
}

const PositionedVFATFrame* PositionedVFATFrameCollection::GetFrameByID(unsigned int ID) const {
    map<unsigned, PositionedVFATFrame*>::const_iterator result = idToFrame.find(ID);
    return result != idToFrame.end() ? result->second : NULL;
}

const PositionedVFATFrame* PositionedVFATFrameCollection::GetFrameByIndex(FramePosition position) const {
    map<FramePosition, PositionedVFATFrame*>::const_iterator result = positionToFrame.find(position);
    return result != positionToFrame.end() ? result->second : NULL;
}

const PositionedVFATFrame* PositionedVFATFrameCollection::GetFrameByIndexID(FramePosition position, unsigned int ID) const{
    const PositionedVFATFrame* res = GetFrameByIndex(position);
    return res != NULL && res->getChipID() == ID ? res : NULL;
}

unsigned int PositionedVFATFrameCollection::Size() const {
    return positionToFrame.size();
}

bool PositionedVFATFrameCollection::Empty() const {
    return Size() == 0;
}

VFATFrameCollection::value_type PositionedVFATFrameCollection::BeginIterator() const
{
    if (positionToFrame.size() > 0)
      return *(positionToFrame.begin());
    else 
      return value_type(FramePosition(), NULL);
}

VFATFrameCollection::value_type PositionedVFATFrameCollection::NextIterator(const value_type& pair) const
{
    map<FramePosition, PositionedVFATFrame*>::const_iterator current = positionToFrame.find(pair.first);
    if (current == positionToFrame.end() || ++current == positionToFrame.end())
      return value_type(FramePosition(), NULL);
    else 
      return *(current++);
}

bool PositionedVFATFrameCollection::IsEndIterator(const value_type& pair) const {
   return positionToFrame.find(pair.first) == positionToFrame.end();
}

void PositionedVFATFrameCollection::AddTriggerFrame(SecondLevelTriggerVFATFrame* frame){
    triggerFrames.push_back(frame);
}
vector<SecondLevelTriggerVFATFrame*>* PositionedVFATFrameCollection::GetTriggerFrames() const{
    return new vector<SecondLevelTriggerVFATFrame*>(triggerFrames.begin(),triggerFrames.end());
}

}

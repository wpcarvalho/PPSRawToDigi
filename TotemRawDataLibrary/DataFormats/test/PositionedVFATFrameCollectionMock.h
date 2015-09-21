/**********************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* Authors:
*   Michal Zmuda (m.zmuda@aol.com)
*
**********************************************************/

#ifndef _MONITOR_QT3_

#include "TotemRawData/DataFormats/interface/PositionedVFATFrameCollection.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"


using namespace std;
using namespace Totem;
using namespace testing;


/**
 * \ingroup TotemRawData
 * TODO: describe.
 **/
class PositionedVFATFrameCollectionMock : public PositionedVFATFrameCollection
{
  public:
    MOCK_CONST_METHOD0(GetClassName, std::string());
    MOCK_CONST_METHOD1(GetFrameByID, const PositionedVFATFrame*(unsigned int ID));
    MOCK_CONST_METHOD1(GetFrameByIndex, const PositionedVFATFrame*(FramePosition index));
    MOCK_CONST_METHOD2(GetFrameByIndexID, const PositionedVFATFrame*(FramePosition index, unsigned int ID));
    MOCK_CONST_METHOD0(Size, unsigned int());
    MOCK_CONST_METHOD0(Empty, bool());
    MOCK_METHOD2(AddNewFrame, void(FramePosition position, PositionedVFATFrame* frame));
    MOCK_METHOD1(AddTriggerFrame, void(SecondLevelTriggerVFATFrame* frame));
    MOCK_CONST_METHOD0(GetTriggerFrames, vector<SecondLevelTriggerVFATFrame*>*());
};

#endif

/**********************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* Authors:
*   Michal Zmuda (m.zmuda@aol.com)
*
**********************************************************/

#ifndef _MONITOR_QT3_

#include "TotemRawData/DataFormats/interface/PositionedVFATFrame.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"


using namespace std;
using namespace Totem;
using namespace testing;

/**
 * \ingroup TotemRawData
 * TODO: describe.
 **/
class PositionedVFATFrameMock : public PositionedVFATFrame
{
  public:
    MOCK_CONST_METHOD0(getRxFlags, PositionedVFATFrame::word());
    MOCK_CONST_METHOD0(getFiberIdx, PositionedVFATFrame::word());
    MOCK_CONST_METHOD0(getGohIdx, PositionedVFATFrame::word());
    MOCK_CONST_METHOD0(getSize, PositionedVFATFrame::word());
    MOCK_CONST_METHOD0(passedHardwareSynchronisationChecks, bool());
    MOCK_CONST_METHOD0(getBC, PositionedVFATFrame::word());
    MOCK_CONST_METHOD0(getEC, PositionedVFATFrame::word());
    MOCK_CONST_METHOD0(getFlags, PositionedVFATFrame::word());
    MOCK_CONST_METHOD0(getChipID, PositionedVFATFrame::word());
    MOCK_CONST_METHOD0(checkFootprint, bool());
    MOCK_CONST_METHOD0(checkCRC, bool());
    MOCK_CONST_METHOD1(channelActive, bool(unsigned char channel));
    MOCK_CONST_METHOD0(getActiveChannels, vector<unsigned char>());
    MOCK_CONST_METHOD1(Print, void(bool));
};

#endif

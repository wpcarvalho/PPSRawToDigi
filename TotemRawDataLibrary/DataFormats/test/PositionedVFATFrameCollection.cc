/**********************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* Authors:
*   Michal Zmuda (m.zmuda@aol.com)
*
**********************************************************/

#ifndef _MONITOR_QT3_

#include "TotemRawData/DataFormats/interface/PositionedVFATFrameCollection.h"
#include "TotemRawData/DataFormats/test/PositionedVFATFrameMock.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <iostream>
#include <algorithm>

using namespace std;
using namespace Totem;
using namespace testing;

TEST(PositionedVFATFrameCollection, test_class_name)
{
    //given
    PositionedVFATFrameCollection sample;

    //then
    EXPECT_EQ(sample.GetClassName(), "PositionedVFATFrameCollection");
}

TEST(PositionedVFATFrameCollection, test_empty_on_create)
{
    //given
    PositionedVFATFrameCollection sample;

    //then
    EXPECT_EQ(sample.Size(), 0);
    EXPECT_TRUE(sample.Empty());
}

TEST(PositionedVFATFrameCollection, test_empty_after_add)
{
    //given
    PositionedVFATFrameCollection sample;
    PositionedVFATFrameMock* mock = new PositionedVFATFrameMock;
    FramePosition position;

    //when
    EXPECT_CALL(*mock, getChipID()).WillRepeatedly(Return(0));
    sample.AddNewFrame(position, mock);

    //then
    EXPECT_EQ(sample.Size(), 1);
    EXPECT_FALSE(sample.Empty());
}

TEST(PositionedVFATFrameCollection, test_find_by_id__only_one_frame_present){
    //given
    PositionedVFATFrameCollection sample;
    PositionedVFATFrameMock* mock = new PositionedVFATFrameMock;
    FramePosition position;

    //when
    EXPECT_CALL(*mock, getChipID()).WillRepeatedly(Return(300));
    sample.AddNewFrame(position, mock);

    //then
    EXPECT_EQ(sample.GetFrameByID(300), mock);
}

TEST(PositionedVFATFrameCollection, test_size_when_multple_added__no_overlaps){
    //given
    PositionedVFATFrameCollection sample;
    PositionedVFATFrameMock* mock1 = new PositionedVFATFrameMock;
    FramePosition position1(0,0,0,0,1);
    PositionedVFATFrameMock* mock2 = new PositionedVFATFrameMock;
    FramePosition position2(0,0,0,0,2);
    PositionedVFATFrameMock* mock3 = new PositionedVFATFrameMock;
    FramePosition position3(0,0,0,1,1);


    //when
    EXPECT_CALL(*mock1, getChipID()).WillRepeatedly(Return(1));
    EXPECT_CALL(*mock2, getChipID()).WillRepeatedly(Return(2));
    EXPECT_CALL(*mock3, getChipID()).WillRepeatedly(Return(3));
    sample.AddNewFrame(position1, mock1);
    sample.AddNewFrame(position2, mock2);
    sample.AddNewFrame(position3, mock3);

    //then
    EXPECT_EQ(sample.Size(), 3);
}

TEST(PositionedVFATFrameCollection, test_size_when_multple_added__position_overlaps){
    //given
    PositionedVFATFrameCollection sample;
    PositionedVFATFrameMock* mock1 = new PositionedVFATFrameMock;
    FramePosition position1(0,0,0,0,1);
    PositionedVFATFrameMock* mock2 = new PositionedVFATFrameMock;
    FramePosition position2(0,0,0,0,2);
    PositionedVFATFrameMock* mock3 = new PositionedVFATFrameMock;
    FramePosition position3(0,0,0,0,1);

    //when
    EXPECT_CALL(*mock1, getChipID()).WillRepeatedly(Return(1));
    EXPECT_CALL(*mock2, getChipID()).WillRepeatedly(Return(2));
    EXPECT_CALL(*mock3, getChipID()).WillRepeatedly(Return(3));
    sample.AddNewFrame(position1, mock1);
    sample.AddNewFrame(position2, mock2);
    sample.AddNewFrame(position3, mock3);

    //then
    EXPECT_EQ(sample.Size(), 2);

    //one of those is not part of collection and not deleted with it
    Mock::AllowLeak(mock1);
    Mock::AllowLeak(mock3);
}

TEST(PositionedVFATFrameCollection, test_size_when_multple_added__id_overlaps_only){
    //given
    PositionedVFATFrameCollection sample;
    PositionedVFATFrameMock* mock1 = new PositionedVFATFrameMock;
    FramePosition position1(0,0,0,0,1);
    PositionedVFATFrameMock* mock2 = new PositionedVFATFrameMock;
    FramePosition position2(0,0,0,0,2);
    PositionedVFATFrameMock* mock3 = new PositionedVFATFrameMock;
    FramePosition position3(0,0,0,0,3);


    //when
    EXPECT_CALL(*mock1, getChipID()).WillRepeatedly(Return(1));
    EXPECT_CALL(*mock2, getChipID()).WillRepeatedly(Return(1));
    EXPECT_CALL(*mock3, getChipID()).WillRepeatedly(Return(1));
    sample.AddNewFrame(position1, mock1);
    sample.AddNewFrame(position2, mock2);
    sample.AddNewFrame(position3, mock3);

    //then
    EXPECT_EQ(sample.Size(), 3);
}

TEST(PositionedVFATFrameCollection, test_getting_by_position__no_overlaps){
    //given
    PositionedVFATFrameCollection sample;
    PositionedVFATFrameMock* mock1 = new PositionedVFATFrameMock;
    FramePosition position1(0,0,0,0,1);
    PositionedVFATFrameMock* mock2 = new PositionedVFATFrameMock;
    FramePosition position2(0,0,0,0,2);
    PositionedVFATFrameMock* mock3 = new PositionedVFATFrameMock;
    FramePosition position3(0,0,0,1,1);


    //when
    EXPECT_CALL(*mock1, getChipID()).WillRepeatedly(Return(1));
    EXPECT_CALL(*mock2, getChipID()).WillRepeatedly(Return(2));
    EXPECT_CALL(*mock3, getChipID()).WillRepeatedly(Return(3));
    sample.AddNewFrame(position1, mock1);
    sample.AddNewFrame(position2, mock2);
    sample.AddNewFrame(position3, mock3);

    //then
    EXPECT_EQ(sample.GetFrameByIndex(position1), mock1);
    EXPECT_EQ(sample.GetFrameByIndex(position2), mock2);
    EXPECT_EQ(sample.GetFrameByIndex(position3), mock3);
}

TEST(PositionedVFATFrameCollection, test_getting_by_position_and_id__no_overlaps){
    //given
    PositionedVFATFrameCollection sample;
    PositionedVFATFrameMock* mock1 = new PositionedVFATFrameMock;
    FramePosition position1(0,0,0,0,1);
    PositionedVFATFrameMock* mock2 = new PositionedVFATFrameMock;
    FramePosition position2(0,0,0,0,2);
    PositionedVFATFrameMock* mock3 = new PositionedVFATFrameMock;
    FramePosition position3(0,0,0,1,1);


    //when
    EXPECT_CALL(*mock1, getChipID()).WillRepeatedly(Return(1));
    EXPECT_CALL(*mock2, getChipID()).WillRepeatedly(Return(2));
    EXPECT_CALL(*mock3, getChipID()).WillRepeatedly(Return(3));
    sample.AddNewFrame(position1, mock1);
    sample.AddNewFrame(position2, mock2);
    sample.AddNewFrame(position3, mock3);

    //then
    EXPECT_EQ(sample.GetFrameByIndexID(position1, 1), mock1);
    EXPECT_EQ(sample.GetFrameByIndexID(position2, 2), mock2);
    EXPECT_EQ(sample.GetFrameByIndexID(position3, 3), mock3);
}


TEST(PositionedVFATFrameCollection, test_getting_by_position_and_id__id_overlaps){
    //given
    PositionedVFATFrameCollection sample;
    PositionedVFATFrameMock* mock1 = new PositionedVFATFrameMock;
    FramePosition position1(0,0,0,0,1);
    PositionedVFATFrameMock* mock2 = new PositionedVFATFrameMock;
    FramePosition position2(0,0,0,0,2);
    PositionedVFATFrameMock* mock3 = new PositionedVFATFrameMock;
    FramePosition position3(0,0,0,1,1);


    //when
    EXPECT_CALL(*mock1, getChipID()).WillRepeatedly(Return(1));
    EXPECT_CALL(*mock2, getChipID()).WillRepeatedly(Return(2));
    EXPECT_CALL(*mock3, getChipID()).WillRepeatedly(Return(1));
    sample.AddNewFrame(position1, mock1);
    sample.AddNewFrame(position2, mock2);
    sample.AddNewFrame(position3, mock3);

    //then
    EXPECT_EQ(sample.GetFrameByIndexID(position1, 1), mock1);
    EXPECT_EQ(sample.GetFrameByIndexID(position2, 2), mock2);
    EXPECT_EQ(sample.GetFrameByIndexID(position3, 1), mock3);
}

TEST(PositionedVFATFrameCollection, test_getting_by_id__no_overlaps){
    //given
    PositionedVFATFrameCollection sample;
    PositionedVFATFrameMock* mock1 = new PositionedVFATFrameMock;
    FramePosition position1(0,0,0,0,1);
    PositionedVFATFrameMock* mock2 = new PositionedVFATFrameMock;
    FramePosition position2(0,0,0,0,2);
    PositionedVFATFrameMock* mock3 = new PositionedVFATFrameMock;
    FramePosition position3(0,0,0,1,1);


    //when
    EXPECT_CALL(*mock1, getChipID()).WillRepeatedly(Return(1));
    EXPECT_CALL(*mock2, getChipID()).WillRepeatedly(Return(2));
    EXPECT_CALL(*mock3, getChipID()).WillRepeatedly(Return(3));
    sample.AddNewFrame(position1, mock1);
    sample.AddNewFrame(position2, mock2);
    sample.AddNewFrame(position3, mock3);

    //then
    EXPECT_EQ(sample.GetFrameByID(1), mock1);
    EXPECT_EQ(sample.GetFrameByID(2), mock2);
    EXPECT_EQ(sample.GetFrameByID(3), mock3);
}

TEST(PositionedVFATFrameCollection, test_getting_by_id__id_overlaps){
    //given
    PositionedVFATFrameCollection sample;
    PositionedVFATFrameMock* mock1 = new PositionedVFATFrameMock;
    FramePosition position1(0,0,0,0,1);
    PositionedVFATFrameMock* mock2 = new PositionedVFATFrameMock;
    FramePosition position2(0,0,0,0,2);
    PositionedVFATFrameMock* mock3 = new PositionedVFATFrameMock;
    FramePosition position3(0,0,0,1,1);


    //when
    EXPECT_CALL(*mock1, getChipID()).WillRepeatedly(Return(1));
    EXPECT_CALL(*mock2, getChipID()).WillRepeatedly(Return(2));
    EXPECT_CALL(*mock3, getChipID()).WillRepeatedly(Return(1));
    sample.AddNewFrame(position1, mock1);
    sample.AddNewFrame(position2, mock2);
    sample.AddNewFrame(position3, mock3);

    //then
    EXPECT_EQ(sample.GetFrameByID(2), mock2);
    ASSERT_THAT(sample.GetFrameByID(1), AnyOf(Eq(mock1), Eq(mock3)));
}

TEST(PositionedVFATFrameCollection, test_getting_by_id__no_match){
    //given
    PositionedVFATFrameCollection sample;
    PositionedVFATFrameMock* mock1 = new PositionedVFATFrameMock;
    FramePosition position1(0,0,0,0,1);

    //when
    EXPECT_CALL(*mock1, getChipID()).WillRepeatedly(Return(1));
    sample.AddNewFrame(position1, mock1);

    //then
    EXPECT_EQ(NULL, sample.GetFrameByID(99));
}

TEST(PositionedVFATFrameCollection, test_getting_by_position__no_match){
    //given
    PositionedVFATFrameCollection sample;
    PositionedVFATFrameMock* mock1 = new PositionedVFATFrameMock;
    FramePosition position1(0,0,0,0,1);

    //when
    EXPECT_CALL(*mock1, getChipID()).WillRepeatedly(Return(1));
    sample.AddNewFrame(position1, mock1);

    //then
    EXPECT_EQ(NULL, sample.GetFrameByIndex(FramePosition(0,0,0,0,0)));
}


TEST(PositionedVFATFrameCollection, test_getting_by_position_and_id__no_match){
    //given
    PositionedVFATFrameCollection sample;
    PositionedVFATFrameMock* mock1 = new PositionedVFATFrameMock;
    FramePosition position1(0,0,0,0,1);

    //when
    EXPECT_CALL(*mock1, getChipID()).WillRepeatedly(Return(1));
    sample.AddNewFrame(position1, mock1);

    //then
    EXPECT_EQ(NULL, sample.GetFrameByIndexID(FramePosition(0,0,0,0,0), mock1->getChipID()));
    EXPECT_EQ(NULL, sample.GetFrameByIndexID(position1, 99));
}

#endif

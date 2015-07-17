/**********************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* Authors:
*   Michal Zmuda (m.zmuda@aol.com)
*
**********************************************************/

#ifndef _MONITOR_QT3_

#include "TotemRawData/DataFormats/interface/RawVFATFrame.h"
#include "gtest/gtest.h"

#include <iostream>
#include <algorithm>

using namespace std;
using namespace Totem;

// ----------------------------- footprint --------------------------

TEST(RawVFATFrame, test_footprint_ok)
{
    //given
    RawVFATFrame::word data[13] = {
       0x9000, //header
       0xA000, //BC
       0xC000, //EC, flags
       0xE000, //ID
       0x0000, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0xF000  //tailer
    };
    RawVFATFrame sample(data);

    //then
    EXPECT_TRUE(sample.checkFootprint());
}

TEST(RawVFATFrame, test_footprint_invalid__frame_header_and_tailer_only)
{
    //given
    RawVFATFrame::word data[13] = {
       0x9000, //header
       0xF000  //tailer
    };
    RawVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.checkFootprint());
    EXPECT_EQ(sample.getEC(), 0);
    EXPECT_EQ(sample.getBC(), 0);
    EXPECT_EQ(sample.getFlags(), 0);
}

TEST(RawVFATFrame, test_rx_validation__passes)
{
    //given
    RawVFATFrame::word data[13] = {
       0x9000, //header
       0xA000, //BC
       0xC000, //EC, flags
       0xE000, //ID
       0x0000, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0xF000  //tailer
    };
    RawVFATFrame sample(data);

    //then
    EXPECT_TRUE(sample.passedHardwareSynchronisationChecks());
}

TEST(RawVFATFrame, test_rx_validation__fails)
{
    //given
    RawVFATFrame::word data[13] = {
       0x9100, //header
       0xA000, //BC
       0xC000, //EC, flags
       0xE000, //ID
       0x0000, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0xF000  //tailer
    };
    RawVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.passedHardwareSynchronisationChecks());
}

TEST(RawVFATFrame, test_footprint_footer_fails)
{
    //given
    RawVFATFrame::word data[13] = {
       0x9000, //header
       0xA000, //BC
       0xC000, //EC, flags
       0xE000, //ID
       0x0000, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0xE000  //tailer
    };
    RawVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.checkFootprint());
}

TEST(RawVFATFrame, test_footprint_header_fails)
{
    //given
    RawVFATFrame::word data[13] = {
       0x8000, //header
       0xA000, //BC
       0xC000, //EC, flags
       0xE000, //ID
       0x0000, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0xF000  //tailer
    };
    RawVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.checkFootprint());
}

TEST(RawVFATFrame, test_footprint_nearBC_fails)
{
    //given
    RawVFATFrame::word data[13] = {
       0x9000, //header
       0xB000, //BC
       0xC000, //EC, flags
       0xE000, //ID
       0x0000, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0xF000  //tailer
    };
    RawVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.checkFootprint());
}

TEST(RawVFATFrame, test_footprint_nearEC_fails)
{
    //given
    RawVFATFrame::word data[13] = {
       0x9000, //header
       0xA000, //BC
       0xD000, //EC, flags
       0xE000, //ID
       0x0000, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0xF000  //tailer
    };
    RawVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.checkFootprint());
}

TEST(RawVFATFrame, test_footprint_nearID_fails)
{
    //given
    RawVFATFrame::word data[13] = {
       0x9000, //header
       0xA000, //BC
       0xD000, //EC, flags
       0xF000, //ID
       0x0000, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0xF000  //tailer
    };
    RawVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.checkFootprint());
}

// ----------------------------- getters --------------------------

TEST(RawVFATFrame, test_getRxFlags)
{
    //given
    RawVFATFrame::word data[13] = {
       0x0F00, //header
       0x0000, //BC
       0x0000, //EC, flags
       0x0000, //ID
       0x0000, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000  //tailer
    };
    RawVFATFrame sample(data);

    //then
    EXPECT_EQ(sample.getRxFlags(), 0xF);
}

TEST(RawVFATFrame, test_getFiberIdx)
{
    //given
    RawVFATFrame::word data[13] = {
       0x00F0, //header
       0x0000, //BC
       0x0000, //EC, flags
       0x0000, //ID
       0x0000, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000  //tailer
    };
    RawVFATFrame sample(data);

    //then
    EXPECT_EQ(sample.getFiberIdx(), 0xF);
}

TEST(RawVFATFrame, test_getGohIdx)
{
    //given
    RawVFATFrame::word data[13] = {
       0x000F, //header
       0x0000, //BC
       0x0000, //EC, flags
       0x0000, //ID
       0x0000, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000  //tailer
    };
    RawVFATFrame sample(data);

    //then
    EXPECT_EQ(sample.getGohIdx(), 0xF);
}

TEST(RawVFATFrame, test_getBC)
{
    //given
    RawVFATFrame::word data[13] = {
       0x0000, //header
       0x0FFF, //BC
       0x0000, //EC, flags
       0x0000, //ID
       0x0000, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000  //tailer
    };
    RawVFATFrame sample(data);

    //then
    EXPECT_EQ(sample.getBC(), 0xFFF);
}

TEST(RawVFATFrame, test_getEC)
{
    //given
    RawVFATFrame::word data[13] = {
       0x0000, //header
       0x0000, //BC
       0x0FF0, //EC, flags
       0x0000, //ID
       0x0000, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000  //tailer
    };
    RawVFATFrame sample(data);

    //then
    EXPECT_EQ(sample.getEC(), 0xFF);
}

TEST(RawVFATFrame, test_getFlags)
{
    //given
    RawVFATFrame::word data[13] = {
       0x0000, //header
       0x0000, //BC
       0x000F, //EC, flags
       0x0000, //ID
       0x0000, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000  //tailer
    };
    RawVFATFrame sample(data);

    //then
    EXPECT_EQ(sample.getFlags(), 0xF);
}

TEST(RawVFATFrame, test_getChipID)
{
    //given
    RawVFATFrame::word data[13] = {
       0x0000, //header
       0x0000, //BC
       0x0000, //EC, flags
       0x0FFF, //ID
       0x0000, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000  //tailer
    };
    RawVFATFrame sample(data);

    //then
    EXPECT_EQ(sample.getChipID(), 0xFFF);
}

TEST(RawVFATFrame, test_getSize)
{
    //given
    RawVFATFrame::word data[13] = {
       0x0000, //header
       0x0000, //BC
       0x0000, //EC, flags
       0x0000, //ID
       0x0000, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x00FF  //tailer
    };
    RawVFATFrame sample(data);

    //then
    EXPECT_EQ(sample.getSize(), 0xFF);
}

// ----------------------------- channels --------------------------

TEST(RawVFATFrame, test_activeChannels_first_active)
{
    //given
    RawVFATFrame::word data[13] = {
       0x0000, //header
       0x0000, //BC
       0x0000, //EC, flags
       0x0000, //ID
       0x0001, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000  //tailer
    };
    RawVFATFrame sample(data);

    //when
    std::vector<unsigned char> hits = sample.getActiveChannels();

    //then
    EXPECT_EQ(hits.size(), 1);
    EXPECT_NE(find(hits.begin(),hits.end(),0), hits.end());
}
TEST(RawVFATFrame, test_activeChannels_last_active)
{
    //given
    RawVFATFrame::word data[13] = {
       0x0000, //header
       0x0000, //BC
       0x0000, //EC, flags
       0x0000, //ID
       0x0000, // channels...
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x0000,
       0x8000,
       0x0000  //tailer
    };
    RawVFATFrame sample(data);

    //when
    std::vector<unsigned char> hits = sample.getActiveChannels();

    //then
    EXPECT_EQ(hits.size(), 1);
    EXPECT_NE(find(hits.begin(),hits.end(),127), hits.end());
}

TEST(RawVFATFrame, test_activeChannels_regular_active)
{
    //given
    RawVFATFrame::word data[13] = {
       0x0000, //header
       0x0000, //BC
       0x0000, //EC, flags
       0x0000, //ID
       0x0001, // channels...
       0x0001,
       0x0001,
       0x0001,
       0x0001,
       0x0001,
       0x0001,
       0x0001,
       0x0000  //tailer
    };
    RawVFATFrame sample(data);

    //when
    std::vector<unsigned char> hits = sample.getActiveChannels();

    //then
    EXPECT_EQ(hits.size(), 8);
    EXPECT_NE(find(hits.begin(),hits.end(),0), hits.end());
    EXPECT_NE(find(hits.begin(),hits.end(),16), hits.end());
    EXPECT_NE(find(hits.begin(),hits.end(),32), hits.end());
    EXPECT_NE(find(hits.begin(),hits.end(),48), hits.end());
    EXPECT_NE(find(hits.begin(),hits.end(),64), hits.end());
    EXPECT_NE(find(hits.begin(),hits.end(),80), hits.end());
    EXPECT_NE(find(hits.begin(),hits.end(),96), hits.end());
    EXPECT_NE(find(hits.begin(),hits.end(),112), hits.end());
}

TEST(RawVFATFrame, test_activeChannels_many_active)
{
    //given
    RawVFATFrame::word data[13] = {
       0x0000, //header
       0x0000, //BC
       0x0000, //EC, flags
       0x0000, //ID
       0x0000, // channels...
       0x0100,
       0x0008,
       0x0020,
       0x0100,
       0x0000,
       0x00F0,
       0x3000,
       0x0000  //tailer
    };
    RawVFATFrame sample(data);

    //when
    std::vector<unsigned char> hits = sample.getActiveChannels();

    //then
    EXPECT_EQ(hits.size(), 10);

    EXPECT_NE(find(hits.begin(),hits.end(),24), hits.end());
    EXPECT_NE(find(hits.begin(),hits.end(),35), hits.end());
    EXPECT_NE(find(hits.begin(),hits.end(),53), hits.end());
    EXPECT_NE(find(hits.begin(),hits.end(),72), hits.end());

    EXPECT_NE(find(hits.begin(),hits.end(),100), hits.end());
    EXPECT_NE(find(hits.begin(),hits.end(),101), hits.end());
    EXPECT_NE(find(hits.begin(),hits.end(),102), hits.end());
    EXPECT_NE(find(hits.begin(),hits.end(),103), hits.end());

    EXPECT_NE(find(hits.begin(),hits.end(),124), hits.end());
    EXPECT_NE(find(hits.begin(),hits.end(),125), hits.end());
}

#endif

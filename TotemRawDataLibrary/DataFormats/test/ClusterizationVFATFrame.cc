/**********************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* Authors:
*   Michal Zmuda (m.zmuda@aol.com)
*
**********************************************************/

#ifndef _MONITOR_QT3_

#include "TotemRawData/DataFormats/interface/ClusterizationVFATFrame.h"
#include "TotemRawData/DataFormats/interface/OptoRxSupplementalData.h" //todo maybe a mock?
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <iostream>
#include <algorithm>

using namespace std;
using namespace Totem;
using namespace testing;

// ----------------------------- footprint --------------------------

TEST(ClusterizationVFATFrame, test_footprint_ok)
{
    //given
    ClusterizationVFATFrame::word data[5] = {
       0x8000, //header
       0xA000, //BC
       0xC000, //EC, flags
       0xE000, //ID
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_TRUE(sample.checkFootprint());
}

TEST(ClusterizationVFATFrame, test_rx_validation__rx_other_fails)
{
    //given
    ClusterizationVFATFrame::word data[5] = {
       0x8100, //header
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.passedHardwareSynchronisationChecks());
}

TEST(ClusterizationVFATFrame, test_rx_validation__ec_fails)
{
    //given
    ClusterizationVFATFrame::word data[5] = {
       0x8200, //header
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.passedHardwareSynchronisationChecks());
}

TEST(ClusterizationVFATFrame, test_rx_validation_bc_fails)
{
    //given
    ClusterizationVFATFrame::word data[5] = {
       0x8400, //header
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.passedHardwareSynchronisationChecks());
}

TEST(ClusterizationVFATFrame, test_rx_validation_passes)
{
    //given
    ClusterizationVFATFrame::word data[5] = {
       0x8000, //header
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_TRUE(sample.passedHardwareSynchronisationChecks());
}

TEST(ClusterizationVFATFrame, test_rx_validation_passes__when_ec_bc_present)
{
//given
ClusterizationVFATFrame::word data[5] = {
        0x8600, //header
        0xA000, //BC
        0xC000, //EC, flags
        0xE000, //ID
        0xE000  //tailer
};
ClusterizationVFATFrame sample(data);

//then
EXPECT_TRUE(sample.passedHardwareSynchronisationChecks());
}

TEST(ClusterizationVFATFrame, test_footprint_footer_fails)
{
    //given
    ClusterizationVFATFrame::word data[5] = {
       0x8000, //header
       0xA000, //BC
       0xC000, //EC, flags
       0xE000, //ID
       0xE000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.checkFootprint());
}

TEST(ClusterizationVFATFrame, test_footprint_header_fails)
{
    //given
    ClusterizationVFATFrame::word data[5] = {
       0x9000, //header
       0xA000, //BC
       0xC000, //EC, flags
       0xE000, //ID
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.checkFootprint());
}

TEST(ClusterizationVFATFrame, test_footprint_nearBC_fails)
{
    //given
    ClusterizationVFATFrame::word data[5] = {
       0x9000, //header
       0xB000, //BC
       0xC000, //EC, flags
       0xE000, //ID
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.checkFootprint());
}

TEST(ClusterizationVFATFrame, test_footprint_nearEC_fails)
{
    //given
    ClusterizationVFATFrame::word data[5] = {
       0x9000, //header
       0xA000, //BC
       0xD000, //EC, flags
       0xE000, //ID
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.checkFootprint());
}

TEST(ClusterizationVFATFrame, test_footprint_nearID_fails)
{
    //given
    ClusterizationVFATFrame::word data[5] = {
       0x9000, //header
       0xA000, //BC
       0xD000, //EC, flags
       0xF000, //ID
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.checkFootprint());
}

// ----------------------------- if has field ---------------------

TEST(ClusterizationVFATFrame, test_hasBC_true)
{
    //given
    ClusterizationVFATFrame::word data[3] = {
       0x0000, //header
       0xAFFF, //BC
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_TRUE(sample.hasBC());
}

TEST(ClusterizationVFATFrame, test_hasBC_false)
{
    //given
    ClusterizationVFATFrame::word data[2] = {
       0x0000, //header
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.hasBC());
}

TEST(ClusterizationVFATFrame, test_hasEC_true)
{
    //given
    ClusterizationVFATFrame::word data[3] = {
       0x0000, //header
       0xCFF0, //EC, flags
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_TRUE(sample.hasEC());
}

TEST(ClusterizationVFATFrame, test_hasEC_false)
{
    //given
    ClusterizationVFATFrame::word data[2] = {
       0x0000, //header
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.hasEC());
}

TEST(ClusterizationVFATFrame, test_hasChipID_true)
{
    //given
    ClusterizationVFATFrame::word data[3] = {
       0x0000, //header
       0xEFFF, //ID
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_TRUE(sample.hasChipID());
}

TEST(ClusterizationVFATFrame, test_hasChipID_false)
{
    //given
    ClusterizationVFATFrame::word data[2] = {
       0x0000, //header
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.hasChipID());
}

// ----------------------------- getters --------------------------

TEST(ClusterizationVFATFrame, test_getRxFlags)
{
    //given
    ClusterizationVFATFrame::word data[2] = {
       0x0F00, //header
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_EQ(sample.getRxFlags(), 0xF);
}

TEST(ClusterizationVFATFrame, test_getFiberIdx)
{
    //given
    ClusterizationVFATFrame::word data[2] = {
       0x00F0, //header
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_EQ(sample.getFiberIdx(), 0xF);
}

TEST(ClusterizationVFATFrame, test_getGohIdx)
{
    //given
    ClusterizationVFATFrame::word data[2] = {
       0x000F, //header
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_EQ(sample.getGohIdx(), 0xF);
}

TEST(ClusterizationVFATFrame, test_getBC)
{
    //given
    ClusterizationVFATFrame::word data[3] = {
       0x0000, //header
       0xAFFF, //BC
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_EQ(sample.getBC(), 0xFFF);
}

TEST(ClusterizationVFATFrame, test_getBC_none)
{
    //given
    ClusterizationVFATFrame::word data[2] = {
       0x0000, //header
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.hasBC());
    EXPECT_EQ(sample.getBC(), 0);
}

TEST(ClusterizationVFATFrame, test_getEC)
{
    //given
    ClusterizationVFATFrame::word data[3] = {
       0x0000, //header
       0xCFF0, //EC, flags
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_EQ(sample.getEC(), 0xFF);
}

TEST(ClusterizationVFATFrame, test_getEC_none)
{
    //given
    ClusterizationVFATFrame::word data[2] = {
       0x0000, //header
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.hasEC());
    EXPECT_EQ(sample.getEC(), 0);
}

TEST(ClusterizationVFATFrame, test_getFlags)
{
    //given
    ClusterizationVFATFrame::word data[3] = {
       0x0000, //header
       0xC00F, //EC, flags
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_EQ(sample.getFlags(), 0xF);
}

TEST(ClusterizationVFATFrame, test_getFlags_none)
{
    //given
    ClusterizationVFATFrame::word data[2] = {
       0x0000, //header
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.hasEC());
    EXPECT_EQ(sample.getFlags(), 0);
}

TEST(ClusterizationVFATFrame, test_getChipID)
{
    //given
    ClusterizationVFATFrame::word data[3] = {
       0x0000, //header
       0xEFFF, //ID
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_EQ(sample.getChipID(), 0xFFF);
}

TEST(ClusterizationVFATFrame, test_getChipID_none)
{
    //given
    ClusterizationVFATFrame::word data[2] = {
       0x0000, //header
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_FALSE(sample.hasChipID());
    EXPECT_EQ(sample.getChipID(), 0);
}

TEST(ClusterizationVFATFrame, test_getSize)
{
    //given
    ClusterizationVFATFrame::word data[2] = {
       0x0000, //header
       0xF0FF  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //then
    EXPECT_EQ(sample.getSize(), 0xFF);
}

// ----------------------------- channels --------------------------

TEST(ClusterizationVFATFrame, test_activeChannels_first_active)
{
    //given
    ClusterizationVFATFrame::word data[3] = {
       0x0000, //header
       0x0100,
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //when
    std::vector<unsigned char> hits = sample.getActiveChannels();

    //then
    EXPECT_EQ(hits.size(), 1);
    EXPECT_THAT(hits, UnorderedElementsAre(0));
}
TEST(ClusterizationVFATFrame, test_activeChannels_last_active)
{
    //given
    ClusterizationVFATFrame::word data[13] = {
       0x0000, //header
       0x017F,
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //when
    std::vector<unsigned char> hits = sample.getActiveChannels();

    //then
    EXPECT_EQ(hits.size(), 1);
    EXPECT_THAT(hits, UnorderedElementsAre(127));
}

TEST(ClusterizationVFATFrame, test_activeChannels_regular_active)
{
    //given
    ClusterizationVFATFrame::word data[13] = {
       0x0000, //header
       0x0100,
       0x0110,
       0x0120,
       0x0130,
       0x0140,
       0x0150,
       0x0160,
       0x0170,
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //when
    std::vector<unsigned char> hits = sample.getActiveChannels();

    //then
    EXPECT_EQ(hits.size(), 8);
    EXPECT_THAT(hits, UnorderedElementsAre(0,16,32,48,64,80,96,112));
}

TEST(ClusterizationVFATFrame, test_activeChannels_many_active)
{
    //given
    ClusterizationVFATFrame::word data[13] = {
       0x0000, //header
       0xA000, //BC
       0xC000, //EC, flags
       0xE000, //ID
       0x0118,
       0x0123,
       0x0135,
       0x0148,
       0x0464,
       0x027C,
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //when
    std::vector<unsigned char> hits = sample.getActiveChannels();

    //then
    EXPECT_EQ(hits.size(), 10);
    EXPECT_THAT(hits, UnorderedElementsAre(24, 35, 53, 72, 100,101,102,103, 124,125));
}

TEST(ClusterizationVFATFrame, test_activeChannels__many_active__get_clusters)
{
    //given
    ClusterizationVFATFrame::word data[13] = {
       0x0000, //header
       0xA000, //BC
       0xC000, //EC, flags
       0xE000, //ID
       0x0118,
       0x0123,
       0x0135,
       0x0148,
       0x0464,
       0x027C,
       0xF000  //tailer
    };
    ClusterizationVFATFrame sample(data);

    //when
    vector<pair<unsigned char, unsigned char>> clusters = sample.getActiveClusters();

    //then
    EXPECT_EQ(clusters.size(), 6);
    EXPECT_THAT(clusters, UnorderedElementsAre(
            Pair(1,24),
            Pair(1,35),
            Pair(1,53),
            Pair(1,72),
            Pair(4,100),
            Pair(2,124)
    ));
}

#endif

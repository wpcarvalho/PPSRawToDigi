/**********************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* Authors:
*   Michal Zmuda (m.zmuda@aol.com)
*
**********************************************************/

#ifndef _MONITOR_QT3_

#include "TotemRawData/DataFormats/interface/OptoRxSupplementalData.h"
#include "gtest/gtest.h"

#include <iostream>

using namespace std;
using namespace Totem;

TEST(OptoRxSupplementalData, test_FOV)
{
    //given
    OptoRxSupplementalData sample(0x0F0, 0x0, 0x0);

    //then
    EXPECT_EQ(sample.getFOV(), 0xF);
    EXPECT_EQ(sample.getOptoRxID(), 0); //neighbor
}

TEST(OptoRxSupplementalData, test_OptoRxID)
{
    //given
    OptoRxSupplementalData sample(0x0FFF00, 0x0, 0x0);

    //then
    EXPECT_EQ(sample.getFOV(), 0); //neighbor
    EXPECT_EQ(sample.getOptoRxID(), 0xFFF);
    EXPECT_EQ(sample.getBC(), 0); //neighbor
}

TEST(OptoRxSupplementalData, test_BC)
{
    //given
    OptoRxSupplementalData sample(0x0FFF00000, 0x0, 0x0);

    //then
    EXPECT_EQ(sample.getOptoRxID(), 0); //neighbor
    EXPECT_EQ(sample.getBC(), 0xFFF);
    EXPECT_EQ(sample.getEC(), 0); //neighbor
}

TEST(OptoRxSupplementalData, test_EC)
{
    //given
    OptoRxSupplementalData sample(0x0FFFFFF00000000, 0x0, 0x0);

    //then
    EXPECT_EQ(sample.getBC(), 0); //neighbor
    EXPECT_EQ(sample.getEC(), 0xFFFFFF);
    EXPECT_EQ(sample.getEvtTy(), 0); //neighbor
}

TEST(OptoRxSupplementalData, test_EvtTy)
{
    //given
    OptoRxSupplementalData sample(0x0F00000000000000, 0x0, 0x0);

    //then
    EXPECT_EQ(sample.getEC(), 0); //neighbor
    EXPECT_EQ(sample.getEvtTy(), 0xF);
    EXPECT_EQ(sample.getBOE(), 0); //neighbor
}

TEST(OptoRxSupplementalData, test_BOE)
{
    //given
    OptoRxSupplementalData sample(0xF000000000000000, 0x0, 0x0);

    //then
    EXPECT_EQ(sample.getEvtTy(), 0); //neighbor
    EXPECT_EQ(sample.getBOE(), 0xF);
}

TEST(OptoRxSupplementalData, test_TTCstat)
{
    //given
    OptoRxSupplementalData sample(0x0, 0x0F0, 0x0);

    //then
    EXPECT_EQ(sample.getCRC(), 0); //neighbor
    EXPECT_EQ(sample.getTTCstat(), 0xF);
}

TEST(OptoRxSupplementalData, test_CRC)
{
    //given
    OptoRxSupplementalData sample(0x0, 0x0FFFF0000, 0x0);

    //then
    EXPECT_EQ(sample.getTTCstat(), 0); //neighbor
    EXPECT_EQ(sample.getCRC(), 0xFFFF);
    EXPECT_EQ(sample.getEventSize(), 0); //neighbor
}

TEST(OptoRxSupplementalData, test_EventSize)
{
    //given
    OptoRxSupplementalData sample(0x0, 0x0FF00000000, 0x0);

    //then
    EXPECT_EQ(sample.getCRC(), 0); //neighbor
    EXPECT_EQ(sample.getEventSize(), 0xFF);
    EXPECT_EQ(sample.getEOE(), 0); //neighbor
}

TEST(OptoRxSupplementalData, test_EOE)
{
    //given
    OptoRxSupplementalData sample(0x0, 0xF000000000000000, 0x0);

    //then
    EXPECT_EQ(sample.getEventSize(), 0); //neighbor
    EXPECT_EQ(sample.getEOE(), 0xF);
}

TEST(OptoRxSupplementalData, test_OrbitCounter)
{
    //given
    OptoRxSupplementalData sample(0x0, 0x0, 0xFFFFFFFF);

    //then
    EXPECT_EQ(sample.getOrbitCounter(), 0xFFFFFFFF);
}

TEST(OptoRxSupplementalData, test_OrbitCounter_overflow)
{
    //given
    OptoRxSupplementalData full(0x0, 0x0, 0xFFFFFFFFFFFFFFFF);
    OptoRxSupplementalData empty(0x0, 0x0, 0xFFFFFFFF00000000);

    //then
    EXPECT_EQ(full.getOrbitCounter(), 0xFFFFFFFF);
    EXPECT_EQ(empty.getOrbitCounter(), 0x0);
}

#endif

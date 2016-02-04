/**********************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* Authors:
*   Michal Zmuda (m.zmuda@aol.com)
*
**********************************************************/

#ifndef _MONITOR_QT3_

#include "TotemRawData/Readers/interface/DynamicOptoRxVMEBFile.h"
#include "TotemRawData/DataFormats/interface/RawEvent.h"
#include "TotemRawData/DataFormats/test/PositionedVFATFrameMock.h"
#include "TotemRawData/DataFormats/test/PositionedVFATFrameCollectionMock.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <algorithm>

using namespace std;
using namespace Totem;
using namespace testing;

struct DynamicOptoRxVMEBFileExposed: public DynamicOptoRxVMEBFile{
    using DynamicOptoRxVMEBFile::ProcessVMEBEvent;
    using DynamicOptoRxVMEBFile::ProcessSubEvent;
    using DynamicOptoRxVMEBFile::ProcessOptoRxFrame;
    using DynamicOptoRxVMEBFile::IsCrcForOptoRxFrameValid;
    using DynamicOptoRxVMEBFile::ProcessRawVFATFrames;
    using DynamicOptoRxVMEBFile::ProcessClusterizationVFATFrames;
    using DynamicOptoRxVMEBFile::ProcessMixedRawAndClusterizationVFATFrames;
};

TEST(Integration_DynamicOptoRxVMEBFile, test_GetClassName){
    //given
    DynamicOptoRxVMEBFileExposed file;

    //then
    EXPECT_EQ("DynamicOptoRxVMEBFile", file.GetClassName());
}

TEST(Integration_DynamicOptoRxVMEBFile, test_empty_optorx){
    //given
    unsigned long long buffer[3] =
            {
                  //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)                            | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                  //OrbitCounter(32)                   | padding
                    0b00000000000000000000000000000llu | 0b00000000000000000000000000000000000llu<<32,
                  //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, _)).Times(0);
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 3*8, &event);
}

TEST(Integration_DynamicOptoRxVMEBFile, test_empty_optorx_crc){
    //given
    //crc generated for 0x8005 poly with http://depa.usst.edu.cn/chenjq/www2/software/crc/CRC_Javascript/CRCcalculation.htm
    unsigned long long buffer[3] =
            {
                  //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)                            | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                  //OrbitCounter(32)                   | padding
                    0b00000000000000000000000000000llu | 0b00000000000000000000000000000000000llu<<32,
                  //0(4)      | TCC(4)       | 0(8)             | CRC(16)    | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0x0000<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;

    //then
    EXPECT_TRUE(DynamicOptoRxVMEBFileExposed::IsCrcForOptoRxFrameValid(buffer, 3*8));
}

TEST(Integration_DynamicOptoRxVMEBFile, test_artificial_optorx_crc){
    //given
    //crc generated for 0x8005 poly with http://depa.usst.edu.cn/chenjq/www2/software/crc/CRC_Javascript/CRCcalculation.htm
    unsigned long long buffer[3] =
            {
                  //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)                            | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                  //OrbitCounter(32)                   | padding
                    0b00000000000000000000000000000llu | 0b00000000000000000000000000000000000llu<<32,
                  //0(4)      | TCC(4)       | 0(8)             | CRC(16)       | EvtSize(10)       | 0(14)                     | x(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0x8005llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0x01llu<<56 /*CHANGE HERE*/
            };
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;

    //then
    EXPECT_TRUE(DynamicOptoRxVMEBFileExposed::IsCrcForOptoRxFrameValid(buffer, 3*8));
}

TEST(Integration_DynamicOptoRxVMEBFile, test_artificial_optorx_crc_2){
    //given
    //crc generated for 0x8005 poly with http://depa.usst.edu.cn/chenjq/www2/software/crc/CRC_Javascript/CRCcalculation.htm
    unsigned long long buffer[3] =
            {
                  //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)                            | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                  //OrbitCounter(32)                   | CHANGE
                    0b00000000000000000000000000000llu | 0xFFFFFFFFllu<<32,
                  //0(4)      | TCC(4)       | 0(8)             | CRC(16)       | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0x28A0llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;

    //then
    EXPECT_TRUE(DynamicOptoRxVMEBFileExposed::IsCrcForOptoRxFrameValid(buffer, 3*8));
}

TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__one_clusterization_frame){
    //given
    unsigned long long buffer[3] =
            {
                  //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)                            | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                  //OrbitCounter(32)                   | header1(16)   | tailer1(16)
                    0b00000000000000000000000000000llu | 0x8000llu<<32 | 0xF002llu<<48,
                  //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, _)).Times(1);
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 3*8, &event);
}

TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__one_clusterization_frame__one_empty_trigger){
    //given
    unsigned long long buffer[4] =

            {
                  //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)                            | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                  //OrbitCounter(32)                   | header1(16)   | tailer1(16)
                    0b00000000000000000000000000000llu | 0x8000llu<<32 | 0xF002llu<<48,
                  // headerTrigger1(16) | tailerTrigger1(16) | padding
                     0xB000llu          | 0xF002llu<<16      | 0llu<<32,
                  //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, _)).Times(1);
    EXPECT_CALL(*collection, AddTriggerFrame(_)).Times(1);
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 4*8, &event);
}

TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__one_clusterization_frame__one_nonempty_trigger){
    //given
    unsigned long long buffer[4] =

            {
                  //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)                            | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                  //OrbitCounter(32)                   | header1(16)   | tailer1(16)
                    0b00000000000000000000000000000llu | 0x8000llu<<32 | 0xF002llu<<48,
                  // headerTrigger1(16) | payload    | tailerTrigger1(16) | padding
                     0xB000llu          | 0x7FFF<<16 | 0xF003llu<<16      | 0llu<<48,
                  //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, _)).Times(1);
    EXPECT_CALL(*collection, AddTriggerFrame(_)).Times(1);
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 4*8, &event);
}

TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__one_empty_trigger__one_clusterization){
    //given
    unsigned long long buffer[4] =

            {
                  //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)                            | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                  //OrbitCounter(32)                   | headerTrigger1(16)   | tailerTrigger1(16)
                    0b00000000000000000000000000000llu | 0xB000llu<<32        | 0xF002llu<<48,
                  // header1(16) | tailer1(16)   | padding
                     0x8000llu   | 0xF002llu<<16 | 0llu<<32,
                  //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, _)).Times(1);
    EXPECT_CALL(*collection, AddTriggerFrame(_)).Times(1);
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 4*8, &event);
}

TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__many_clusterization_frames){
    //given
    unsigned long long buffer[4] =
            {
                  //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)                            | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                  //OrbitCounter(32)                   | header1(16)   | tailer1(16)
                    0b00000000000000000000000000000llu | 0x8000llu<<32 | 0xF002llu<<48,
                  //header2(16)| tailer2(16)   | header3(16)   | tailer3(16)
                    0x8001llu  | 0xF002llu<<16 | 0x8002llu<<32 | 0xF002llu<<48,
                  //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, _)).Times(3);
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 4*8, &event);
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__many_clusterization_frames__padded_zeros){
    //given
    unsigned long long buffer[4] =
            {
                  //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)                            | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                  //OrbitCounter(32)                   | header1(16)   | tailer1(16)
                    0b00000000000000000000000000000llu | 0x8000llu<<32 | 0xF002llu<<48,
                  //header2(16)| tailer2(16)   | padding
                    0x8001llu  | 0xF002llu<<16 | 0b00000000000000000000000000000llu<<32,
                  //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, _)).Times(2);
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 4*8, &event);
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__many_clusterization_frames__fail_rx){
    //given
    unsigned long long buffer[3] =
            {
                  //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)                            | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                  //OrbitCounter(32)                   | header1(16)   | tailer1(16)
                    0b00000000000000000000000000000llu | 0x8100llu<<32 | 0xF002llu<<48,
                  //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, _)).Times(0);
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 3*8, &event);
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__many_clusterization_frames__fail_ec__ec_bc_absent){
    //given
    unsigned long long buffer[3] =
            {
                  //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)                            | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                  //OrbitCounter(32)                   | header1(16)   | tailer1(16)
                    0b00000000000000000000000000000llu | 0x8400llu<<32 | 0xF002llu<<48,
                  //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, _)).Times(0);
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 3*8, &event);
}

TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__many_clusterization_frames__fail_ec__ec_present_so_passes){
    //given
    unsigned long long EC = 34;
    unsigned long long buffer[4] =
            {
                    //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)   | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                    //OrbitCounter(32)                   | header1(16)   | EC word (16)
                    0b00000000000000000000000000000llu | 0x8200llu<<32 | EC<<52 | 0b1100llu<<60,
                    //tailer1(16) | padding
                    0xF003llu   | 0llu<<16,
                    //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, _)).Times(1);
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 4*8, &event);
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__many_clusterization_frames__fail_bc__ec_bc_absent){
    //given
    unsigned long long buffer[3] =
            {
                  //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)                            | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                  //OrbitCounter(32)                   | header1(16)   | tailer1(16)
                    0b00000000000000000000000000000llu | 0x8200llu<<32 | 0xF002llu<<48,
                  //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, _)).Times(0);
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 3*8, &event);
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__many_clusterization_frames__fail_bc__bc_present_so_passes){
    //given
    unsigned long long BC = 34;
    unsigned long long buffer[4] =
            {
                    //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)   | EC(24)                            | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                    //OrbitCounter(32)                   | header1(16)   | BC word (16)
                    0b00000000000000000000000000000llu | 0x8400llu<<32 | BC<<48 | 0b1010llu<<60,
                    //tailer1(16) | padding
                    0xF002llu   | 0llu<<16,
                    //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, _)).Times(1);
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 4*8, &event);
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__many_clusterization_frames__padded_with_random_but_no_start_of_frame){
    //given
    unsigned long long buffer[4] =
            {
                  //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)                            | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                  //OrbitCounter(32)                   | header1(16)   | tailer1(16)
                    0b00000000000000000000000000000llu | 0x8000llu<<32 | 0xF002llu<<48,
                  //header2(16)| tailer2(16)   | padding
                    0x8001llu  | 0xF002llu<<16 | 0x532EF410llu<<32,
                  //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, _)).Times(2);
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 4*8, &event);
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__many_clusterization_frames__faulty_part_with_random_but_no_start_of_frame){
    //given
    unsigned long long buffer[4] =
            {
                  //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)                            | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                  //OrbitCounter(32)                   | header1(16)   | tailer1(16)
                    0b00000000000000000000000000000llu | 0x8000llu<<32 | 0xF002llu<<48,
                  // MESS             | header2(16)    | tailer2(16)
                    0x532EF410llu | 0x8001llu<<32  | 0xF002llu<<48 ,
                  //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, _)).Times(2);
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 4*8, &event);
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__nonempty_clusterization_frames__clusters_presence){
    //given
    unsigned long long buffer[4] =
            {
                  //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)                            | Evt(4)        | BOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                  //OrbitCounter(32)                   | header1(16)   | cluster(16)
                    0b00000000000000000000000000000llu | 0x8000llu<<32 | 0x0148llu<<48,
                  //cluster(16)  | cluster(16)    | tailer1(16)   | padding
                    0x0464llu    | 0x027Cllu<<16  | 0xF002llu<<32 | 0llu<<48,
                  //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                    0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollection* collection = new PositionedVFATFrameCollection;
    RawEvent event;
    event.frames = collection;

    //when
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 4*8, &event);

    //then
    const PositionedVFATFrame* frame = collection->GetFrameByIndex(FramePosition(0,0,0,0,0));
    EXPECT_EQ(1, collection->Size());
    EXPECT_TRUE(frame != NULL);
    EXPECT_THAT(frame->getActiveChannels(), UnorderedElementsAre(72, 100,101,102,103, 124,125));
}

TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__empty_clusterization_frame){
    //given
    unsigned long long buffer[3] =
            {
                    //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)                            | Evt(4)        | BOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                    //OrbitCounter(32)                   | header1(16)   | tailer1(16)
                      0b00000000000000000000000000000llu | 0x8000llu<<32 | 0xF002llu<<48,
                    //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollection* collection = new PositionedVFATFrameCollection;
    RawEvent event;
    event.frames = collection;

    //when
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 3*8, &event);

    //then
    const PositionedVFATFrame* frame = collection->GetFrameByIndex(FramePosition(0,0,0,0,0));
    EXPECT_EQ(1, collection->Size());
    EXPECT_TRUE(frame != NULL);
    EXPECT_THAT(frame->getActiveChannels(), UnorderedElementsAre());
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__no_data_in_clusterization_frame__check_ec){
    //given
    unsigned long long EC = 34;
    unsigned long long buffer[3] =
            {
                    //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24) | Evt(4)        | BOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | EC<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                    //OrbitCounter(32)                   | header1(16)   | tailer1(16)
                      0b00000000000000000000000000000llu | 0x8000llu<<32 | 0xF002llu<<48,
                    //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollection* collection = new PositionedVFATFrameCollection;
    RawEvent event;
    event.frames = collection;

    //when
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 3*8, &event);

    //then
    const PositionedVFATFrame* frame = collection->GetFrameByIndex(FramePosition(0,0,0,0,0));
    EXPECT_TRUE(frame != NULL);
    EXPECT_EQ(0, frame->getEC());
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__no_data_in_clusterization_frame__check_ec_present){
    //given
    unsigned long long EC = 34;
    unsigned long long buffer[4] =
            {
                    //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)                | EC(24)   | Evt(4)        | BOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0b000000000000llu<<20 | 0llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                    //OrbitCounter(32)                   | header1(16)   | EC word (16)
                      0b00000000000000000000000000000llu | 0x8000llu<<32 | 0llu<<48 | EC<<52 | 0b1100llu<<60,
                    //tailer1(16) | padding
                      0xF003llu   | 0llu<<16,
                    //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollection* collection = new PositionedVFATFrameCollection;
    RawEvent event;
    event.frames = collection;

    //when
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 4*8, &event);

    //then
    const PositionedVFATFrame* frame = collection->GetFrameByIndex(FramePosition(0,0,0,0,0));
    EXPECT_TRUE(frame != NULL);
    EXPECT_EQ(EC, frame->getEC());
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__no_data_in_clusterization_frame__check_bc){
    //given
    unsigned long long BC = 34;
    unsigned long long buffer[3] =
            {
                    //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12) | EC(24)                            | Evt(4)        | BOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | BC<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                    //OrbitCounter(32)                   | header1(16)   | tailer1(16)
                      0b00000000000000000000000000000llu | 0x8000llu<<32 | 0xF002llu<<48,
                    //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollection* collection = new PositionedVFATFrameCollection;
    RawEvent event;
    event.frames = collection;

    //when
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 3*8, &event);

    //then
    const PositionedVFATFrame* frame = collection->GetFrameByIndex(FramePosition(0,0,0,0,0));
    EXPECT_TRUE(frame != NULL);
    EXPECT_EQ(0, frame->getBC());
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__no_data_in_clusterization_frame__check_bc_present){
    //given
    unsigned long long BC = 34;
    unsigned long long buffer[4] =
            {
                    //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)   | EC(24)                            | Evt(4)        | BOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                    //OrbitCounter(32)                   | header1(16)   | BC word (16)
                      0b00000000000000000000000000000llu | 0x8000llu<<32 | BC<<48 | 0b1010llu<<60,
                    //tailer1(16) | padding
                      0xF003llu   | 0llu<<16,
                    //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollection* collection = new PositionedVFATFrameCollection;
    RawEvent event;
    event.frames = collection;

    //when
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 4*8, &event);

    //then
    const PositionedVFATFrame* frame = collection->GetFrameByIndex(FramePosition(0,0,0,0,0));
    EXPECT_TRUE(frame != NULL);
    EXPECT_EQ(BC, frame->getBC());
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__no_data_in_clusterization_frame__check_id_present){
    //given
    unsigned long long ID = 34;
    unsigned long long buffer[4] =
            {
                    //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)   | EC(24)                            | Evt(4)        | BOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                    //OrbitCounter(32)                   | header1(16)   | EC word (16)
                      0b00000000000000000000000000000llu | 0x8000llu<<32 | ID<<48 | 0b1110llu<<60,
                    //tailer1(16) | padding
                      0xF002llu   | 0llu<<16,
                    //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollection* collection = new PositionedVFATFrameCollection;
    RawEvent event;
    event.frames = collection;

    //when
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 4*8, &event);

    //then
    const PositionedVFATFrame* frame = collection->GetFrameByIndex(FramePosition(0,0,0,0,0));
    EXPECT_TRUE(frame != NULL);
    EXPECT_EQ(ID, frame->getChipID());
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__no_data_in_clusterization_frame__check_id_default){
    //given
    unsigned long long buffer[3] =
            {
                    //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)   | EC(24)                            | Evt(4)        | BOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                    //OrbitCounter(32)                   | header1(16)   | tailer1(16)
                      0b00000000000000000000000000000llu | 0x8000llu<<32 | 0xF002llu<<48,
                    //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollection* collection = new PositionedVFATFrameCollection;
    RawEvent event;
    event.frames = collection;

    //when
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 3*8, &event);

    //then
    const PositionedVFATFrame* frame = collection->GetFrameByIndex(FramePosition(0,0,0,0,0));
    EXPECT_TRUE(frame != NULL);
    EXPECT_EQ(0, frame->getChipID());
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__raw_frames__one_valid_goh){
    //given
    unsigned long long buffer[55] =
            {
                    //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)   | EC(24)                            | Evt(4)        | BOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                    //OrbitCounter(32)                   | header0(16)   | BCword(16)
                      0b00000000000000000000000000000llu | 0x9000llu<<32 | 0xA000llu<<48,
                    //ECword(16)|IDword(16)     | 32 channels
                      0xC000llu | 0xE000llu<<16 | 0llu<<32,
                    //64 channels
                      0llu,
                    //32 channels   | tailer0(16)    | header1(16)
                      0llu          | 0xE00Dllu<<32  | 0x9010llu<<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                      0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                      0llu,
                    //48 channels | tailer1(16)
                      0llu        | 0xE00Dllu<<48,
                    // header2(16) |BCword(16)    | ECword(16)    |IDword(16)
                      0x9020llu    |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                      0llu,
                    //64 channels
                      0llu,
                    //tailer2(16)    | header3(16)     | BCword(16)    | ECword(16)
                      0xE00Dllu      | 0x9030llu<<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                      0xE000llu  | 0llu<<16,
                    //64 channels
                      0llu,
                    //16 channels | tailer3(16)    | header4(16)   | BCword(16)
                      0llu        | 0xE00Dllu<<16  | 0x9040llu<<32 | 0xA000llu<<48,
                    //ECword(16)    | IDword(16)     | 32 channels
                      0xC000llu     | 0xE000llu<<16  | 0llu<<32,
                    //64 channels
                      0llu,
                    //32 channels   | tailer4(16)    | header5(16)
                      0llu          | 0xE00Dllu<<32  | 0x9050llu<<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                      0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                      0llu,
                    //48 channels | tailer5(16)
                      0llu        | 0xE00Dllu<<48,
                    // header6(16) |BCword(16)    | ECword(16)    |IDword(16)
                      0x9060llu    |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                      0llu,
                    //64 channels
                      0llu,
                    //tailer6(16) | header7(16)     | BCword(16)    | ECword(16)
                      0xE00Dllu   | 0x9070llu<<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                      0xE000llu  | 0llu<<16,
                    //64 channels
                      0llu,
                    //16 channels | tailer7(16)    | header8(16)   | BCword(16)
                      0llu        | 0xE00Dllu<<16  | 0x9080llu<<32 | 0xA000llu<<48,
                    //ECword(16)    | IDword(16)     | 32 channels
                      0xC000llu     | 0xE000llu<<16  | 0llu<<32,
                    //64 channels
                      0llu,
                    //32 channels   | tailer8(16) | header9(16)
                      0llu          | 0xE00Dllu<<32  | 0x9090llu<<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                      0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                      0llu,
                    //48 channels | tailer9(16)
                      0llu        | 0xE00Dllu<<48,
                    // headerA(16) |BCword(16)    | ECword(16)    |IDword(16)
                      0x90A0llu    |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                      0llu,
                    //64 channels
                      0llu,
                    //tailerA(16) | headerB(16)     | BCword(16)    | ECword(16)
                      0xE00Dllu   | 0x90B0llu<<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                      0xE000llu  | 0llu<<16,
                    //64 channels
                      0llu,
                    //16 channels | tailerB(16)    | headerC(16)   | BCword(16)
                      0llu        | 0xE00Dllu<<16  | 0x90C0llu<<32 | 0xA000llu<<48,
                    //ECword(16)    | IDword(16)     | 32 channels
                      0xC000llu     | 0xE000llu<<16  | 0llu<<32,
                    //64 channels
                      0llu,
                    //32 channels | tailerC(16) | headerD(16)
                      0llu        | 0xE00Dllu<<32  | 0x90D0llu<<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                      0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                      0llu,
                    //48 channels | tailerD(16)
                      0llu        | 0xE00Dllu<<48,
                    // headerE(16) |BCword(16)    | ECword(16)    |IDword(16)
                      0x90E0llu    |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                      0llu,
                    //64 channels
                      0llu,
                    //tailerE(16) | headerF(16)     | BCword(16)    | ECword(16)
                      0xE00Dllu   | 0x90F0llu<<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                      0xE000llu  | 0llu<<16,
                    //64 channels
                      0llu,
                    //16 channels | tailerF(16)    | padding
                      0llu        | 0xE00Dllu<<16  | 0llu<<32,
                    //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollection* collection = new PositionedVFATFrameCollection;
    RawEvent event;
    event.frames = collection;

    //when
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 55*8, &event);

    //then
    EXPECT_EQ(16, collection->Size());
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__raw_frames__one_almost_valid_goh__one_wrong_EC){
    //given
    unsigned long long buffer[55] =
            {
                    //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)   | EC(24)                            | Evt(4)        | BOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                    //OrbitCounter(32)                   | header0(16)   | BCword(16)
                      0b00000000000000000000000000000llu | 0x9000llu<<32 | 0xA000llu<<48,
                    //ECword(16)|IDword(16)     | 32 channels
                      0xC000llu | 0xE000llu<<16 | 0llu<<32,
                    //64 channels
                      0llu,
                    //32 channels   | tailer0(16)    | header1(16)
                      0llu          | 0xE00Dllu<<32  | 0x9010llu<<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                      0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                      0llu,
                    //48 channels | tailer1(16)
                      0llu        | 0xE00Dllu<<48,
                    // header2(16) |BCword(16)    | ECword(16)    |IDword(16)
                      0x9020llu    |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                      0llu,
                    //64 channels
                      0llu,
                    //tailer2(16)    | header3(16)     | BCword(16)    | ECword(16)
                      0xE00Dllu      | 0x9030llu<<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                      0xE000llu  | 0llu<<16,
                    //64 channels
                      0llu,
                    //16 channels | tailer3(16)    | header4(16)   | BCword(16)
                      0llu        | 0xE00Dllu<<16  | 0x9040llu<<32 | 0xA000llu<<48,
                    //ECword(16)    | IDword(16)     | 32 channels
                      0xCFF0llu     | 0xE000llu<<16  | 0llu<<32,                            /*WRONG EC */
                    //64 channels
                      0llu,
                    //32 channels   | tailer4(16)    | header5(16)
                      0llu          | 0xE00Dllu<<32  | 0x9050llu<<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                      0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                      0llu,
                    //48 channels | tailer5(16)
                      0llu        | 0xE00Dllu<<48,
                    // header6(16) |BCword(16)    | ECword(16)    |IDword(16)
                      0x9060llu    |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                      0llu,
                    //64 channels
                      0llu,
                    //tailer6(16) | header7(16)     | BCword(16)    | ECword(16)
                      0xE00Dllu   | 0x9070llu<<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                      0xE000llu  | 0llu<<16,
                    //64 channels
                      0llu,
                    //16 channels | tailer7(16)    | header8(16)   | BCword(16)
                      0llu        | 0xE00Dllu<<16  | 0x9080llu<<32 | 0xA000llu<<48,
                    //ECword(16)    | IDword(16)     | 32 channels
                      0xC000llu     | 0xE000llu<<16  | 0llu<<32,
                    //64 channels
                      0llu,
                    //32 channels   | tailer8(16) | header9(16)
                      0llu          | 0xE00Dllu<<32  | 0x9090llu<<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                      0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                      0llu,
                    //48 channels | tailer9(16)
                      0llu        | 0xE00Dllu<<48,
                    // headerA(16) |BCword(16)    | ECword(16)    |IDword(16)
                      0x90A0llu    |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                      0llu,
                    //64 channels
                      0llu,
                    //tailerA(16) | headerB(16)     | BCword(16)    | ECword(16)
                      0xE00Dllu   | 0x90B0llu<<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                      0xE000llu  | 0llu<<16,
                    //64 channels
                      0llu,
                    //16 channels | tailerB(16)    | headerC(16)   | BCword(16)
                      0llu        | 0xE00Dllu<<16  | 0x90C0llu<<32 | 0xA000llu<<48,
                    //ECword(16)    | IDword(16)     | 32 channels
                      0xC000llu     | 0xE000llu<<16  | 0llu<<32,
                    //64 channels
                      0llu,
                    //32 channels | tailerC(16) | headerD(16)
                      0llu        | 0xE00Dllu<<32  | 0x90D0llu<<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                      0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                      0llu,
                    //48 channels | tailerD(16)
                      0llu        | 0xE00Dllu<<48,
                    // headerE(16) |BCword(16)    | ECword(16)    |IDword(16)
                      0x90E0llu    |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                      0llu,
                    //64 channels
                      0llu,
                    //tailerE(16) | headerF(16)     | BCword(16)    | ECword(16)
                      0xE00Dllu   | 0x90F0llu<<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                      0xE000llu  | 0llu<<16,
                    //64 channels
                      0llu,
                    //16 channels | tailerF(16)    | padding
                      0llu        | 0xE00Dllu<<16  | 0llu<<32,
                    //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollection* collection = new PositionedVFATFrameCollection;
    RawEvent event;
    event.frames = collection;

    //when
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 55*8, &event);

    //then
    EXPECT_EQ(0, collection->Size());
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__raw_frames__one_almost_valid_goh__one_wrong_BC){
    //given
    unsigned long long buffer[55] =
            {
                    //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)   | EC(24)                            | Evt(4)        | BOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                    //OrbitCounter(32)                   | header0(16)   | BCword(16)
                      0b00000000000000000000000000000llu | 0x9000llu<<32 | 0xA000llu<<48,
                    //ECword(16)|IDword(16)     | 32 channels
                      0xC000llu | 0xE000llu<<16 | 0llu<<32,
                    //64 channels
                      0llu,
                    //32 channels   | tailer0(16)    | header1(16)
                      0llu          | 0xE00Dllu<<32  | 0x9010llu<<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                      0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                      0llu,
                    //48 channels | tailer1(16)
                      0llu        | 0xE00Dllu<<48,
                    // header2(16) |BCword(16)    | ECword(16)    |IDword(16)
                      0x9020llu    |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                      0llu,
                    //64 channels
                      0llu,
                    //tailer2(16)    | header3(16)     | BCword(16)    | ECword(16)
                      0xE00Dllu      | 0x9030llu<<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                      0xE000llu  | 0llu<<16,
                    //64 channels
                      0llu,
                    //16 channels | tailer3(16)    | header4(16)   | BCword(16)
                      0llu        | 0xE00Dllu<<16  | 0x9040llu<<32 | 0xA000llu<<48,
                    //ECword(16)    | IDword(16)     | 32 channels
                      0xC000llu     | 0xE000llu<<16  | 0llu<<32,
                    //64 channels
                      0llu,
                    //32 channels   | tailer4(16)    | header5(16)
                      0llu          | 0xE00Dllu<<32  | 0x9050llu<<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                      0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                      0llu,
                    //48 channels | tailer5(16)
                      0llu        | 0xE00Dllu<<48,
                    // header6(16) |BCword(16)    | ECword(16)    |IDword(16)
                      0x9060llu    |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                      0llu,
                    //64 channels
                      0llu,
                    //tailer6(16) | header7(16)     | BCword(16)    | ECword(16)
                      0xE00Dllu   | 0x9070llu<<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                      0xE000llu  | 0llu<<16,
                    //64 channels
                      0llu,
                    //16 channels | tailer7(16)    | header8(16)   | BCword(16)
                      0llu        | 0xE00Dllu<<16  | 0x9080llu<<32 | 0xAFFFllu<<48,             /*WRONG BC */
                    //ECword(16)    | IDword(16)     | 32 channels
                      0xC000llu     | 0xE000llu<<16  | 0llu<<32,
                    //64 channels
                      0llu,
                    //32 channels   | tailer8(16) | header9(16)
                      0llu          | 0xE00Dllu<<32  | 0x9090llu<<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                      0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                      0llu,
                    //48 channels | tailer9(16)
                      0llu        | 0xE00Dllu<<48,
                    // headerA(16) |BCword(16)    | ECword(16)    |IDword(16)
                      0x90A0llu    |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                      0llu,
                    //64 channels
                      0llu,
                    //tailerA(16) | headerB(16)     | BCword(16)    | ECword(16)
                      0xE00Dllu   | 0x90B0llu<<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                      0xE000llu  | 0llu<<16,
                    //64 channels
                      0llu,
                    //16 channels | tailerB(16)    | headerC(16)   | BCword(16)
                      0llu        | 0xE00Dllu<<16  | 0x90C0llu<<32 | 0xA000llu<<48,
                    //ECword(16)    | IDword(16)     | 32 channels
                      0xC000llu     | 0xE000llu<<16  | 0llu<<32,
                    //64 channels
                      0llu,
                    //32 channels | tailerC(16) | headerD(16)
                      0llu        | 0xE00Dllu<<32  | 0x90D0llu<<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                      0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                      0llu,
                    //48 channels | tailerD(16)
                      0llu        | 0xE00Dllu<<48,
                    // headerE(16) |BCword(16)    | ECword(16)    |IDword(16)
                      0x90E0llu    |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                      0llu,
                    //64 channels
                      0llu,
                    //tailerE(16) | headerF(16)     | BCword(16)    | ECword(16)
                      0xE00Dllu   | 0x90F0llu<<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                      0xE000llu  | 0llu<<16,
                    //64 channels
                      0llu,
                    //16 channels | tailerF(16)    | padding
                      0llu        | 0xE00Dllu<<16  | 0llu<<32,
                    //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollection* collection = new PositionedVFATFrameCollection;
    RawEvent event;
    event.frames = collection;

    //when
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 55*8, &event);

    //then
    EXPECT_EQ(0, collection->Size());
}


TEST(Integration_DynamicOptoRxVMEBFile, test_optorx__raw_frames__one_valid_goh_surrounded_with_not_complete_gohs){
    //given
    unsigned long long buffer[81] =
            {
                    //0(4)      |FOV(4)        |OptoRxID(12)          | BC(12)   | EC(24)                            | Evt(4)        | BOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b000000000000llu<<8 | 0llu<<20 | 0b000000000000000000000000llu<<32 | 0b0000llu<<56 | 0b0000llu<<60,
                    //OrbitCounter(32)                   | header0(16)   | BCword(16)
                      0b00000000000000000000000000000llu | 0x9000llu<<32 | 0xA000llu<<48,
                    //ECword(16)|IDword(16)     | 32 channels
                      0xC000llu | 0xE000llu<<16 | 0llu<<32,
                    //64 channels
                      0llu,
                    //32 channels   | tailer0(16)    | header1(16)
                      0llu          | 0xE00Dllu<<32  | 0x9010llu<<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                      0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                      0llu,
                    //48 channels | tailer1(16)
                      0llu        | 0xE00Dllu<<48,
                    // header2(16) |BCword(16)    | ECword(16)    |IDword(16)
                      0x9020llu    |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                      0llu,
                    //64 channels
                      0llu,
                    //tailer2(16)    | header3(16)     | BCword(16)    | ECword(16)
                      0xE00Dllu      | 0x9030llu<<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                      0xE000llu  | 0llu<<16,
                    //64 channels
                      0llu,
                    //16 channels | tailer3(16)    | header4(16)   | BCword(16)
                      0llu        | 0xE00Dllu<<16  | 0x9040llu<<32 | 0xA000llu<<48,
                    //ECword(16)    | IDword(16)     | 32 channels
                      0xC000llu     | 0xE000llu<<16  | 0llu<<32,
                    //64 channels
                      0llu,
                    //32 channels   | tailer4(16)    | header5(16)
                      0llu          | 0xE00Dllu<<32  | 0x9050llu<<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                      0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                      0llu,
                    //48 channels | tailer5(16)
                      0llu        | 0xE00Dllu<<48,

                    //INVALID START
                    // header2(16) |BCword(16)    | ECword(16)    |IDword(16)
                    0x9021llu     |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                    0llu,
                    //64 channels
                    0llu,
                    //tailer2(16)    | header3(16)     | BCword(16)    | ECword(16)
                    0xE00Dllu      | 0x9031llu <<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                    0xE000llu  | 0llu<<16,
                    //64 channels
                    0llu,
                    //16 channels | tailer3(16)    | header4(16)   | BCword(16)
                    0llu        | 0xE00Dllu<<16  | 0x9041llu <<32 | 0xA000llu<<48,
                    //ECword(16)    | IDword(16)     | 32 channels
                    0xC000llu     | 0xE000llu<<16  | 0llu<<32,
                    //64 channels
                    0llu,
                    //32 channels   | tailer4(16)    | header5(16)
                    0llu          | 0xE00Dllu<<32  | 0x9051llu <<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                    0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                    0llu,
                    //48 channels | tailer5(16)
                    0llu        | 0xE00Dllu<<48,
                    //INVALID END

                    // header6(16) |BCword(16)    | ECword(16)    |IDword(16)
                      0x9060llu    |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                      0llu,
                    //64 channels
                      0llu,
                    //tailer6(16) | header7(16)     | BCword(16)    | ECword(16)
                      0xE00Dllu   | 0x9070llu<<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                      0xE000llu  | 0llu<<16,
                    //64 channels
                      0llu,
                    //16 channels | tailer7(16)    | header8(16)   | BCword(16)
                      0llu        | 0xE00Dllu<<16  | 0x9080llu<<32 | 0xA000llu<<48,
                    //ECword(16)    | IDword(16)     | 32 channels
                      0xC000llu     | 0xE000llu<<16  | 0llu<<32,
                    //64 channels
                      0llu,
                    //32 channels   | tailer8(16) | header9(16)
                      0llu          | 0xE00Dllu<<32  | 0x9090llu<<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                      0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                      0llu,
                    //48 channels | tailer9(16)
                      0llu        | 0xE00Dllu<<48,
                    // headerA(16) |BCword(16)    | ECword(16)    |IDword(16)
                      0x90A0llu    |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                      0llu,
                    //64 channels
                      0llu,
                    //tailerA(16) | headerB(16)     | BCword(16)    | ECword(16)
                      0xE00Dllu   | 0x90B0llu<<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                      0xE000llu  | 0llu<<16,
                    //64 channels
                      0llu,
                    //16 channels | tailerB(16)    | headerC(16)   | BCword(16)
                      0llu        | 0xE00Dllu<<16  | 0x90C0llu<<32 | 0xA000llu<<48,
                    //ECword(16)    | IDword(16)     | 32 channels
                      0xC000llu     | 0xE000llu<<16  | 0llu<<32,
                    //64 channels
                      0llu,
                    //32 channels | tailerC(16) | headerD(16)
                      0llu        | 0xE00Dllu<<32  | 0x90D0llu<<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                      0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                      0llu,
                    //48 channels | tailerD(16)
                      0llu        | 0xE00Dllu<<48,

                    //INVALID START
                    // header2(16) |BCword(16)    | ECword(16)    |IDword(16)
                    0x9022llu      |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                    0llu,
                    //64 channels
                    0llu,
                    //tailer2(16)    | header3(16)     | BCword(16)    | ECword(16)
                    0xE00Dllu      | 0x9032llu  <<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                    0xE000llu  | 0llu<<16,
                    //64 channels
                    0llu,
                    //16 channels | tailer3(16)    | header4(16)   | BCword(16)
                    0llu        | 0xE00Dllu<<16  | 0x9042llu  <<32 | 0xA000llu<<48,
                    //ECword(16)    | IDword(16)     | 32 channels
                    0xC000llu     | 0xE000llu<<16  | 0llu<<32,
                    //64 channels
                    0llu,
                    //32 channels   | tailer4(16)    | header5(16)
                    0llu          | 0xE00Dllu<<32  | 0x9052llu  <<48,
                    //BCword(16)| ECword(16)    |IDword(16)     | 16 channels
                    0xA000llu | 0xC000llu<<16 | 0xE000llu<<32 | 0llu<<48,
                    //64 channels
                    0llu,
                    //48 channels | tailer5(16)
                    0llu        | 0xE00Dllu<<48,
                    //INVALID END

                    // headerE(16) |BCword(16)    | ECword(16)    |IDword(16)
                      0x90E0llu    |0xA000llu<<16 | 0xC000llu<<32 | 0xE000llu<<48,
                    //64 channels
                      0llu,
                    //64 channels
                      0llu,
                    //tailerE(16) | headerF(16)     | BCword(16)    | ECword(16)
                      0xE00Dllu   | 0x90F0llu<<16   | 0xA000llu<<32 | 0xC000llu<<48,
                    //IDword(16) | 48 channels
                      0xE000llu  | 0llu<<16,
                    //64 channels
                      0llu,
                    //16 channels | tailerF(16)    | padding
                      0llu        | 0xE00Dllu<<16  | 0llu<<32,
                    //0(4)      | TCC(4)       | 0(8)             | CRC(16)                   | EvtSize(10)       | 0(14)                     | x(4)          | EOE(4)
                      0b0000llu | 0b0000llu<<4 | 0b00000000llu<<8 | 0b0000000000000000llu<<16 | 0b00000000llu<<32 | 0b0000000000000000llu<<42 | 0b0000llu<<56 | 0b0000llu<<60,
            };
    PositionedVFATFrameCollection* collection = new PositionedVFATFrameCollection;
    RawEvent event;
    event.frames = collection;

    //when
    DynamicOptoRxVMEBFileExposed::ProcessOptoRxFrame(buffer, 81*8, &event);

    //then
    EXPECT_EQ(16, collection->Size());
}

#endif

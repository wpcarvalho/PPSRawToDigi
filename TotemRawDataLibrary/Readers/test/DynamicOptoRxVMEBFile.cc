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

TEST(DynamicOptoRxVMEBFile, test_GetClassName){
    //given
    DynamicOptoRxVMEBFileExposed file;

    //then
    EXPECT_EQ("DynamicOptoRxVMEBFile", file.GetClassName());
}

TEST(DynamicOptoRxVMEBFile, test_ProcessClusterizationVFATFrames_validates_clusterization_checks){
    PositionedVFATFrameMock* valid = new PositionedVFATFrameMock;
    PositionedVFATFrameMock* invalid = new PositionedVFATFrameMock;
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);

    vector<PositionedVFATFrame*> both;
    both.push_back(invalid);
    both.push_back(valid);

    //when
    EXPECT_CALL(*invalid, getEC()).WillRepeatedly(Return(0));
    EXPECT_CALL(*invalid, getBC()).WillRepeatedly(Return(0));
    EXPECT_CALL(*valid, getEC()).WillRepeatedly(Return(0));
    EXPECT_CALL(*valid, getBC()).WillRepeatedly(Return(0));
    EXPECT_CALL(*valid, getGohIdx()).WillRepeatedly(Return(0));
    EXPECT_CALL(*valid, getFiberIdx()).WillRepeatedly(Return(0));
    EXPECT_CALL(*invalid, getGohIdx()).WillRepeatedly(Return(0));
    EXPECT_CALL(*invalid, getFiberIdx()).WillRepeatedly(Return(0));
    EXPECT_CALL(*valid, passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));
    EXPECT_CALL(*invalid, passedHardwareSynchronisationChecks()).WillRepeatedly(Return(false));

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, valid)).Times(1);
    EXPECT_CALL(*collection, AddNewFrame(_, invalid)).Times(0);
    DynamicOptoRxVMEBFileExposed::ProcessClusterizationVFATFrames(both, &event, &opto);

    //cleanup
    delete invalid;
    delete valid;
}

TEST(DynamicOptoRxVMEBFile, test_ProcessClusterizationVFATFrames_validates_ec_and_bc_consistiency__passed){
    //given
    PositionedVFATFrameMock* valid = new PositionedVFATFrameMock;
    PositionedVFATFrameMock* anotherValid = new PositionedVFATFrameMock;
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);

    vector<PositionedVFATFrame*> both;
    both.push_back(anotherValid);
    both.push_back(valid);

    //when
    EXPECT_CALL(*valid, getGohIdx()).WillRepeatedly(Return(0));
    EXPECT_CALL(*valid, getFiberIdx()).WillRepeatedly(Return(0));
    EXPECT_CALL(*anotherValid, getGohIdx()).WillRepeatedly(Return(0));
    EXPECT_CALL(*anotherValid, getFiberIdx()).WillRepeatedly(Return(1));
    EXPECT_CALL(*anotherValid, getEC()).WillRepeatedly(Return(90));
    EXPECT_CALL(*anotherValid, getBC()).WillRepeatedly(Return(10));
    EXPECT_CALL(*valid, getEC()).WillRepeatedly(Return(90));
    EXPECT_CALL(*valid, getBC()).WillRepeatedly(Return(10));
    EXPECT_CALL(*valid, passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));
    EXPECT_CALL(*anotherValid, passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, valid)).Times(1);
    EXPECT_CALL(*collection, AddNewFrame(_, anotherValid)).Times(1);
    DynamicOptoRxVMEBFileExposed::ProcessClusterizationVFATFrames(both, &event, &opto);

    //cleanup
    delete anotherValid;
    delete valid;
}

TEST(DynamicOptoRxVMEBFile, test_ProcessClusterizationVFATFrames_validates_ec_and_bc_consistiency__broken){
    //given
    PositionedVFATFrameMock* anotherInvalid = new PositionedVFATFrameMock;
    PositionedVFATFrameMock* invalid = new PositionedVFATFrameMock;
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);

    vector<PositionedVFATFrame*> both;
    both.push_back(invalid);
    both.push_back(anotherInvalid);

    //when
    EXPECT_CALL(*anotherInvalid, getGohIdx()).WillRepeatedly(Return(0));
    EXPECT_CALL(*anotherInvalid, getFiberIdx()).WillRepeatedly(Return(1));
    EXPECT_CALL(*invalid, getGohIdx()).WillRepeatedly(Return(0));
    EXPECT_CALL(*invalid, getFiberIdx()).WillRepeatedly(Return(0));
    EXPECT_CALL(*invalid, getEC()).WillRepeatedly(Return(190));
    EXPECT_CALL(*invalid, getBC()).WillRepeatedly(Return(110));
    EXPECT_CALL(*anotherInvalid, getEC()).WillRepeatedly(Return(90));
    EXPECT_CALL(*anotherInvalid, getBC()).WillRepeatedly(Return(10));
    EXPECT_CALL(*anotherInvalid, passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));
    EXPECT_CALL(*invalid, passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, _)).Times(0);
    DynamicOptoRxVMEBFileExposed::ProcessClusterizationVFATFrames(both, &event, &opto);

    //cleanup
    delete invalid;
    delete anotherInvalid;
}


TEST(DynamicOptoRxVMEBFile, test_ProcessMixedRawAndClusterizationVFATFrames_no_validation_at_all){
    //given
    PositionedVFATFrameMock* valid = new PositionedVFATFrameMock;
    PositionedVFATFrameMock* invalid = new PositionedVFATFrameMock;
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);

    vector<PositionedVFATFrame*> both;
    both.push_back(invalid);
    both.push_back(valid);

    //when
    EXPECT_CALL(*valid, getGohIdx()).WillRepeatedly(Return(0));
    EXPECT_CALL(*valid, getFiberIdx()).WillRepeatedly(Return(0));
    EXPECT_CALL(*valid, getEC()).WillRepeatedly(Return(99));
    EXPECT_CALL(*valid, getBC()).WillRepeatedly(Return(199));
    EXPECT_CALL(*invalid, getGohIdx()).WillRepeatedly(Return(32));
    EXPECT_CALL(*invalid, getFiberIdx()).WillRepeatedly(Return(41));
    EXPECT_CALL(*invalid, getEC()).WillRepeatedly(Return(1));
    EXPECT_CALL(*invalid, getBC()).WillRepeatedly(Return(13));
    EXPECT_CALL(*valid, passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));
    EXPECT_CALL(*invalid, passedHardwareSynchronisationChecks()).WillRepeatedly(Return(false));

    //then
    EXPECT_CALL(*collection, AddNewFrame(_, valid)).Times(1);
    EXPECT_CALL(*collection, AddNewFrame(_, invalid)).Times(1);
    DynamicOptoRxVMEBFileExposed::ProcessMixedRawAndClusterizationVFATFrames(both, &event, &opto);

    //cleanup
    delete invalid;
    delete valid;
}

TEST(DynamicOptoRxVMEBFile, test_ProcessRawVFATFrames_validates_fiber_occupancy__no_missing){
    //given
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    PositionedVFATFrameMock frames[16];
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);

    vector<PositionedVFATFrame*> all;
    for(int i = 0; i < 16; i++){
        //when
        all.push_back(&(frames[i]));
        EXPECT_CALL(frames[i], getGohIdx()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], getFiberIdx()).WillRepeatedly(Return(i));
        EXPECT_CALL(frames[i], getEC()).WillRepeatedly(Return(1));
        EXPECT_CALL(frames[i], getBC()).WillRepeatedly(Return(1));
        EXPECT_CALL(frames[i], passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));

        //then
        EXPECT_CALL(*collection, AddNewFrame(_, &(frames[i]))).Times(1);
    }

    //when
    DynamicOptoRxVMEBFileExposed::ProcessRawVFATFrames(all, &event, &opto);
}

TEST(DynamicOptoRxVMEBFile, test_ProcessRawVFATFrames_validates_fiber_occupancy__one_missing){
    //given
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    PositionedVFATFrameMock frames[15];
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);

    vector<PositionedVFATFrame*> all;
    for(int i = 0; i < 15; i++){
        //when
        all.push_back(&(frames[i]));
        EXPECT_CALL(frames[i], getGohIdx()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], getFiberIdx()).WillRepeatedly(Return(i));
        EXPECT_CALL(frames[i], getEC()).WillRepeatedly(Return(1));
        EXPECT_CALL(frames[i], getBC()).WillRepeatedly(Return(1));
        EXPECT_CALL(frames[i], passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));

        //then
        EXPECT_CALL(*collection, AddNewFrame(_, &(frames[i]))).Times(0);
    }

    //when
    DynamicOptoRxVMEBFileExposed::ProcessRawVFATFrames(all, &event, &opto);
}


TEST(DynamicOptoRxVMEBFile, test_ProcessRawVFATFrames_validates_fiber_occupancy__all_fiber_ids_but_in_different_goh){
    //given
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    PositionedVFATFrameMock frames[16];
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);

    vector<PositionedVFATFrame*> all;
    for(int i = 0; i < 16; i++){
        //when
        all.push_back(&(frames[i]));
        EXPECT_CALL(frames[i], getGohIdx()).WillRepeatedly(Return(i%2));
        EXPECT_CALL(frames[i], getFiberIdx()).WillRepeatedly(Return(i));
        EXPECT_CALL(frames[i], getEC()).WillRepeatedly(Return(1));
        EXPECT_CALL(frames[i], getBC()).WillRepeatedly(Return(1));
        EXPECT_CALL(frames[i], passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));

        //then
        EXPECT_CALL(*collection, AddNewFrame(_, &(frames[i]))).Times(0);
    }

    //when
    DynamicOptoRxVMEBFileExposed::ProcessRawVFATFrames(all, &event, &opto);
}

TEST(DynamicOptoRxVMEBFile, test_ProcessRawVFATFrames_validates_fiber_occupancy__one_goh_fine_one_incomplete){
    //given
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    PositionedVFATFrameMock frames[19];
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);

    vector<PositionedVFATFrame*> all;
    for(int i = 0; i < 19; i++){
        //when
        all.push_back(&(frames[i]));
        EXPECT_CALL(frames[i], getGohIdx()).WillRepeatedly(Return(i/16));
        EXPECT_CALL(frames[i], getFiberIdx()).WillRepeatedly(Return(i%16));
        EXPECT_CALL(frames[i], getEC()).WillRepeatedly(Return(1));
        EXPECT_CALL(frames[i], getBC()).WillRepeatedly(Return(1));
        EXPECT_CALL(frames[i], passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));

        //then
        EXPECT_CALL(*collection, AddNewFrame(_, &(frames[i]))).Times(i<16?1:0);
    }

    //when
    DynamicOptoRxVMEBFileExposed::ProcessRawVFATFrames(all, &event, &opto);
}

//--------------------------------- test EC match ---------------------------

TEST(DynamicOptoRxVMEBFile, test_ProcessRawVFATFrames_validates_ec_and_bc_consistency__ec_all_matching){
    //given
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    PositionedVFATFrameMock frames[16];
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);//note that ec & bc from this header do not match frame values!

    vector<PositionedVFATFrame*> all;
    for(int i = 0; i < 16; i++){
        //when
        all.push_back(&(frames[i]));
        EXPECT_CALL(frames[i], getGohIdx()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], getFiberIdx()).WillRepeatedly(Return(i));
        EXPECT_CALL(frames[i], getEC()).WillRepeatedly(Return(99));
        EXPECT_CALL(frames[i], getBC()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));

        //then
        EXPECT_CALL(*collection, AddNewFrame(_, &(frames[i]))).Times(1);
    }

    //when
    DynamicOptoRxVMEBFileExposed::ProcessRawVFATFrames(all, &event, &opto);
}

TEST(DynamicOptoRxVMEBFile, test_ProcessRawVFATFrames_validates_ec_and_bc_consistency__ec_first_not_matching){
    //given
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    PositionedVFATFrameMock frames[16];
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);//note that ec & bc from this header do not match frame values!

    vector<PositionedVFATFrame*> all;
    for(int i = 0; i < 16; i++){
        //when
        all.push_back(&(frames[i]));
        EXPECT_CALL(frames[i], getGohIdx()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], getFiberIdx()).WillRepeatedly(Return(i));
        EXPECT_CALL(frames[i], getEC()).WillRepeatedly(Return(i==0?98:99));
        EXPECT_CALL(frames[i], getBC()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));

        //then
        EXPECT_CALL(*collection, AddNewFrame(_, &(frames[i]))).Times(0);
    }

    //when
    DynamicOptoRxVMEBFileExposed::ProcessRawVFATFrames(all, &event, &opto);
}

TEST(DynamicOptoRxVMEBFile, test_ProcessRawVFATFrames_validates_ec_and_bc_consistency__ec_last_not_matching){
    //given
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    PositionedVFATFrameMock frames[16];
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);//note that ec & bc from this header do not match frame values!

    vector<PositionedVFATFrame*> all;
    for(int i = 0; i < 16; i++){
        //when
        all.push_back(&(frames[i]));
        EXPECT_CALL(frames[i], getGohIdx()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], getFiberIdx()).WillRepeatedly(Return(i));
        EXPECT_CALL(frames[i], getEC()).WillRepeatedly(Return(i==15?98:99));
        EXPECT_CALL(frames[i], getBC()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));

        //then
        EXPECT_CALL(*collection, AddNewFrame(_, &(frames[i]))).Times(0);
    }

    //when
    DynamicOptoRxVMEBFileExposed::ProcessRawVFATFrames(all, &event, &opto);
}

TEST(DynamicOptoRxVMEBFile, test_ProcessRawVFATFrames_validates_ec_and_bc_consistency__ec_middle_not_matching){
    //given
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    PositionedVFATFrameMock frames[16];
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);//note that ec & bc from this header do not match frame values!

    vector<PositionedVFATFrame*> all;
    for(int i = 0; i < 16; i++){
        //when
        all.push_back(&(frames[i]));
        EXPECT_CALL(frames[i], getGohIdx()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], getFiberIdx()).WillRepeatedly(Return(i));
        EXPECT_CALL(frames[i], getEC()).WillRepeatedly(Return(i==5?98:99));
        EXPECT_CALL(frames[i], getBC()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));

        //then
        EXPECT_CALL(*collection, AddNewFrame(_, &(frames[i]))).Times(0);
    }

    //when
    DynamicOptoRxVMEBFileExposed::ProcessRawVFATFrames(all, &event, &opto);
}

TEST(DynamicOptoRxVMEBFile, test_ProcessRawVFATFrames_validates_ec_and_bc_consistency__ec_one_goh_matching_one_not){
    //given
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    PositionedVFATFrameMock frames[32];
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);//note that ec & bc from this header do not match frame values!

    vector<PositionedVFATFrame*> all;
    for(int i = 0; i < 32; i++){
        //when
        all.push_back(&(frames[i]));
        EXPECT_CALL(frames[i], getGohIdx()).WillRepeatedly(Return(i/16));
        EXPECT_CALL(frames[i], getFiberIdx()).WillRepeatedly(Return(i%16));
        //first is not ok
        EXPECT_CALL(frames[i], getEC()).WillRepeatedly(Return(i==0?98:99));
        EXPECT_CALL(frames[i], getBC()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));

        //then
        EXPECT_CALL(*collection, AddNewFrame(_, &(frames[i]))).Times(i<16?0:1);
    }

    //when
    DynamicOptoRxVMEBFileExposed::ProcessRawVFATFrames(all, &event, &opto);
}

//--------------------------------- same with BC ---------------------------

TEST(DynamicOptoRxVMEBFile, test_ProcessRawVFATFrames_validates_ec_and_bc_consistency__ec_bc_all_matching){
    //given
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    PositionedVFATFrameMock frames[16];
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);//note that ec & bc from this header do not match frame values!

    vector<PositionedVFATFrame*> all;
    for(int i = 0; i < 16; i++){
        //when
        all.push_back(&(frames[i]));
        EXPECT_CALL(frames[i], getGohIdx()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], getFiberIdx()).WillRepeatedly(Return(i));
        EXPECT_CALL(frames[i], getEC()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], getBC()).WillRepeatedly(Return(111));
        EXPECT_CALL(frames[i], passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));

        //then
        EXPECT_CALL(*collection, AddNewFrame(_, &(frames[i]))).Times(1);
    }

    //when
    DynamicOptoRxVMEBFileExposed::ProcessRawVFATFrames(all, &event, &opto);
}

TEST(DynamicOptoRxVMEBFile, test_ProcessRawVFATFrames_validates_ec_and_bc_consistency__ec_bc_first_not_matching){
    //given
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    PositionedVFATFrameMock frames[16];
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);//note that ec & bc from this header do not match frame values!

    vector<PositionedVFATFrame*> all;
    for(int i = 0; i < 16; i++){
        //when
        all.push_back(&(frames[i]));
        EXPECT_CALL(frames[i], getGohIdx()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], getFiberIdx()).WillRepeatedly(Return(i));
        EXPECT_CALL(frames[i], getEC()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], getBC()).WillRepeatedly(Return(i==0?198:199));
        EXPECT_CALL(frames[i], passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));

        //then
        EXPECT_CALL(*collection, AddNewFrame(_, &(frames[i]))).Times(0);
    }

    //when
    DynamicOptoRxVMEBFileExposed::ProcessRawVFATFrames(all, &event, &opto);
}

TEST(DynamicOptoRxVMEBFile, test_ProcessRawVFATFrames_validates_ec_and_bc_consistency__ec_bc_last_not_matching){
    //given
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    PositionedVFATFrameMock frames[16];
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);//note that ec & bc from this header do not match frame values!

    vector<PositionedVFATFrame*> all;
    for(int i = 0; i < 16; i++){
        //when
        all.push_back(&(frames[i]));
        EXPECT_CALL(frames[i], getGohIdx()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], getFiberIdx()).WillRepeatedly(Return(i));
        EXPECT_CALL(frames[i], getEC()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], getBC()).WillRepeatedly(Return(i==15?198:199));
        EXPECT_CALL(frames[i], passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));

        //then
        EXPECT_CALL(*collection, AddNewFrame(_, &(frames[i]))).Times(0);
    }

    //when
    DynamicOptoRxVMEBFileExposed::ProcessRawVFATFrames(all, &event, &opto);
}

TEST(DynamicOptoRxVMEBFile, test_ProcessRawVFATFrames_validates_ec_and_bc_consistency__ec_bc_middle_not_matching){
    //given
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    PositionedVFATFrameMock frames[16];
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);//note that ec & bc from this header do not match frame values!

    vector<PositionedVFATFrame*> all;
    for(int i = 0; i < 16; i++){
        //when
        all.push_back(&(frames[i]));
        EXPECT_CALL(frames[i], getGohIdx()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], getFiberIdx()).WillRepeatedly(Return(i));
        EXPECT_CALL(frames[i], getEC()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], getBC()).WillRepeatedly(Return(i==5?198:199));
        EXPECT_CALL(frames[i], passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));

        //then
        EXPECT_CALL(*collection, AddNewFrame(_, &(frames[i]))).Times(0);
    }

    //when
    DynamicOptoRxVMEBFileExposed::ProcessRawVFATFrames(all, &event, &opto);
}

TEST(DynamicOptoRxVMEBFile, test_ProcessRawVFATFrames_validates_ec_and_bc_consistency__ec_bc_one_goh_matching_one_not){
    //given
    PositionedVFATFrameCollectionMock* collection = new PositionedVFATFrameCollectionMock;
    PositionedVFATFrameMock frames[32];
    RawEvent event;
    event.frames = collection;
    OptoRxSupplementalData opto(0x000000000, 0x0, 0x0);//note that ec & bc from this header do not match frame values!

    vector<PositionedVFATFrame*> all;
    for(int i = 0; i < 32; i++){
        //when
        all.push_back(&(frames[i]));
        EXPECT_CALL(frames[i], getGohIdx()).WillRepeatedly(Return(i/16));
        EXPECT_CALL(frames[i], getFiberIdx()).WillRepeatedly(Return(i%16));
        //first is not ok
        EXPECT_CALL(frames[i], getEC()).WillRepeatedly(Return(0));
        EXPECT_CALL(frames[i], getBC()).WillRepeatedly(Return(i==0?198:199));
        EXPECT_CALL(frames[i], passedHardwareSynchronisationChecks()).WillRepeatedly(Return(true));

        //then
        EXPECT_CALL(*collection, AddNewFrame(_, &(frames[i]))).Times(i<16?0:1);
    }

    //when
    DynamicOptoRxVMEBFileExposed::ProcessRawVFATFrames(all, &event, &opto);
}

#endif

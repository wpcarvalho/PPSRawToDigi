/****************************************************************************
*
* This is a part of the TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors: 
*
****************************************************************************/

#include "TotemRawDataLibrary/DataFormats/interface/CommonDef.h"
#include "TotemRawDataLibrary/DataFormats/interface/FramePosition.h"
#include "TotemRawDataLibrary/Readers/interface/DynamicOptoRxVMEBFile.h"

// the RunII version
#include "TotemRawDataLibrary/DAQA/interface/event_3_14.h"

#include <cmath>
#include <set>

using namespace std;


namespace Totem {

//----------------------------------------------------------------------------------------------------

DynamicOptoRxVMEBFile::DynamicOptoRxVMEBFile() : BaseVMEBFile()
{
#ifdef DEBUG
    printf(">> DynamicOptoRxVMEBFile::DynamicOptoRxVMEBFile \n");
#endif
}

//----------------------------------------------------------------------------------------------------

DynamicOptoRxVMEBFile::~DynamicOptoRxVMEBFile()
{
#ifdef DEBUG
    printf(">> DynamicOptoRxVMEBFile::~DynamicOptoRxVMEBFile, this = %p, dataPtr = %p\n", (void *) this, dataPtr);
#endif
}

//----------------------------------------------------------------------------------------------------

unsigned char DynamicOptoRxVMEBFile::GetNextEvent(RawEvent *event)
{
	int retCode = BaseVMEBFile::GetNextEvent(event);
    if(retCode) return retCode;

    // process the buffer
    unsigned int errorCounter = ProcessVMEBEvent(dataPtr, event);

    if (errorCounter > 0) {
        ERROR("DynamicOptoRxVMEBFile::GetNextEvent") << errorCounter << " broken frame (when assumed all frames are ok) in event " << event->dataEventNumber << "." << c_endl;
        corruptedEventCounter++;
    }

    return 0;
}

//----------------------------------------------------------------------------------------------------

unsigned int DynamicOptoRxVMEBFile::ProcessVMEBEvent(char *ptr, RawEvent *event)
{
    eventHeaderStruct *eventHeader = (eventHeaderStruct *) ptr;
    bool superEvent = TEST_ANY_ATTRIBUTE(eventHeader->eventTypeAttribute, ATTR_SUPER_EVENT);

#ifdef DEBUG
    printf(">> DynamicOptoRxVMEBFile::ProcessVMEBEvent\n");

    printf("\teventSize = %i\n", eventHeader->eventSize);
    printf("\teventMagic = %i\n", eventHeader->eventMagic);
    printf("\teventHeadSize = %i\n", eventHeader->eventHeadSize);
    printf("\t* eventHeadSize (extra) = %i\n", eventHeader->eventHeadSize - EVENT_HEAD_BASE_SIZE);
    printf("\t* eventPayloadSize = %i\n", eventHeader->eventSize - eventHeader->eventHeadSize);
    printf("\teventVersion = %i\n", eventHeader->eventVersion);
    printf("\teventType = %i\n", eventHeader->eventType);
    printf("\teventRunNb = %i\n", eventHeader->eventRunNb);
    printf("\teventId = %p\n", (void*)eventHeader->eventId);
    printf("\teventTriggerPattern = %p\n", (void*)eventHeader->eventTriggerPattern);
    printf("\teventDetectorPattern = %p\n",(void*) eventHeader->eventDetectorPattern);
    printf("\teventTypeAttribute = %p\n", (void*)eventHeader->eventTypeAttribute);
    printf("\t* super event = %i\n", superEvent);
    printf("\teventLdcId = %i\n", eventHeader->eventLdcId);
    printf("\teventGdcId = %i\n", eventHeader->eventGdcId);
    printf("\teventTimestamp = %i\n", eventHeader->eventTimestamp);
#endif
  
    // store important GDC data
    event->dataEventNumber = EVENT_ID_GET_NB_IN_RUN(eventHeader->eventId) - 1;
    event->timestamp = eventHeader->eventTimestamp;

    eventSizeType eventSize = eventHeader->eventSize;
    eventHeadSizeType headSize = eventHeader->eventHeadSize;

    // process all sub-events (LDC frames)
    unsigned int errorCounter = 0;
    if (superEvent) {
        unsigned int offset = headSize;
        while (offset < eventSize) {
#ifdef DEBUG 
            printf("\t> offset before %i\n", offset);
#endif
            eventStruct *subEvPtr = (eventStruct *) (ptr + offset);
            eventSizeType subEvSize = subEvPtr->eventHeader.eventSize;
            errorCounter += ProcessSubEvent(ptr + offset, event);
            offset += subEvSize;
#ifdef DEBUG 
            printf("\t> offset after %i\n", offset);
#endif
        }
    } else
    errorCounter += ProcessSubEvent(ptr, event);
    return errorCounter;
}

//----------------------------------------------------------------------------------------------------

unsigned int DynamicOptoRxVMEBFile::ProcessSubEvent(char *ptr, RawEvent *event)
{
    eventHeaderStruct *eventHeader = (eventHeaderStruct *) ptr;

#ifdef DEBUG 
    printf("\t\t>> ProcessSubEvent\n");

    printf("\t\teventSize = %u\n", eventHeader->eventSize);
    printf("\t\teventMagic = %u\n", eventHeader->eventMagic);
    printf("\t\teventHeadSize = %u\n", eventHeader->eventHeadSize);
    printf("\t\t* eventHeadSize (extra) = %u\n", eventHeader->eventHeadSize - EVENT_HEAD_BASE_SIZE);
    printf("\t\t* eventPayloadSize = %u\n", eventHeader->eventSize - eventHeader->eventHeadSize);
    printf("\t\teventVersion = %u\n", eventHeader->eventVersion);
    printf("\t\teventType = %u\n", eventHeader->eventType);
    printf("\t\teventRunNb = %u\n", eventHeader->eventRunNb);
    printf("\t\teventId = %p\n", (void*) eventHeader->eventId);
    printf("\t\teventTriggerPattern = %p\n", (void*) eventHeader->eventTriggerPattern);
    printf("\t\teventDetectorPattern = %p\n", (void*) eventHeader->eventDetectorPattern);
    printf("\t\teventTypeAttribute = %p\n", (void*) eventHeader->eventTypeAttribute);
    printf("\t\teventLdcId = %u\n", eventHeader->eventLdcId);
    printf("\t\teventGdcId = %u\n", eventHeader->eventGdcId);
    printf("\t\teventTimestamp = %u\n", eventHeader->eventTimestamp);
#endif

    // store important LDC data
    event->ldcTimeStamps[eventHeader->eventLdcId] = eventHeader->eventTimestamp;

    unsigned long subEvSize = eventHeader->eventSize;
    unsigned long offset = eventHeader->eventHeadSize;

    // process all equipments (OptoRx frames)
    unsigned int errorCounter = 0;
    while (offset < subEvSize) {
#ifdef DEBUG 
        printf("\t\toffset (before) %lu\n", offset);
#endif
        equipmentHeaderStruct *eq = (equipmentHeaderStruct *) (ptr + offset);
        equipmentSizeType equipmentHeaderStructSize = sizeof(equipmentHeaderStruct);
        signed long payloadSize = eq->equipmentSize - equipmentHeaderStructSize;
        unsigned long long *payloadPtr = (unsigned long long *)(ptr + offset + equipmentHeaderStructSize);
#ifdef DEBUG 
        printf("\t\t\tequipmentSize = %u\n", eq->equipmentSize);
        printf("\t\t\tequipmentType = %u\n", eq->equipmentType);
        printf("\t\t\tequipmentId = %u\n", eq->equipmentId);
        printf("\t\t\tequipmentTypeAttribute = %p\n", (void*) eq->equipmentTypeAttribute);
        printf("\t\t\tequipmentBasicElementSize = %u\n", eq->equipmentBasicElementSize);
        printf("\t\t\t\t\tpayload size = %li\n", payloadSize);
        printf("\t\t\t\t\tpayload ptr = %p\n", (void*) payloadPtr);
#endif

        if (payloadSize > 0) {
            payloadSize /= 8; // bytes -> words

            switch (eq->equipmentType) {
                case etDynamicOptoRx:
                    if(!IsCrcForOptoRxFrameValid(payloadPtr, payloadSize)){
                        unsigned optoIdx = (payloadPtr[0] >> 8) & 0xFFF;
                        ERROR("DynamicOptoRxVMEBFile::ProcessSubEvent,CRC check (for OptoRxId" << optoIdx << ") has failed!");
                    } else {
                        errorCounter += ProcessOptoRxFrame(payloadPtr, payloadSize, event);
                    }
                    break;
                case etLoneG:
                    ProcessLoneGFrame(payloadPtr, payloadSize, event);
                    break;
                default:
                    ERROR("DynamicOptoRxVMEBFile::ProcessSubEvent") << "Unknown equipment type: " << eq->equipmentType << ". Skipping." << c_endl;
            }
        }
        offset += eq->equipmentSize;
#ifdef DEBUG 
        printf("\t\toffset (after) %lu\n", offset);
#endif
    }
    return errorCounter;
}

//----------------------------------------------------------------------------------------------------

unsigned int DynamicOptoRxVMEBFile::ProcessOptoRxFrame(unsigned long long* buffer, unsigned int frameSize, RawEvent *event)
{
    unsigned bufferLength = ceil((double)frameSize/sizeof(*buffer));
    OptoRxSupplementalData* optoRx = new OptoRxSupplementalData(buffer[0], buffer[bufferLength-1], buffer[1] && 0xFFFFFFFF);

    VFATFrame::word* wordsBuffer = (VFATFrame::word*) buffer;
    //skipping beginning (header data + orbit counter)
    unsigned offsetWords = 6;

    int rawFramesCount = 0;
    int clusterizationFramesCount = 0;
    vector<PositionedVFATFrame*> framesRead;

    while(offsetWords < bufferLength*sizeof(*buffer)/sizeof(VFATFrame::word)){
        int frameHeaderStamp = (wordsBuffer[offsetWords] >> 12) & 0xF;
        VFATFrame::word* shiftedBuffer =  &(wordsBuffer[offsetWords]);

        if(frameHeaderStamp == 0x9){
            PositionedVFATFrame* frame = new RawVFATFrame(shiftedBuffer);
            framesRead.push_back(frame);
            offsetWords += frame->getSize();
            rawFramesCount++;
        } else if(frameHeaderStamp == 0x8) {
            PositionedVFATFrame* frame = new ClusterizationVFATFrame(shiftedBuffer);
            framesRead.push_back(frame);
            offsetWords += frame->getSize();
            clusterizationFramesCount++;
        } else if(frameHeaderStamp == 0x0) {
            //probably padding
            offsetWords += 1;
        } else if(frameHeaderStamp == 0xB) {
            SecondLevelTriggerVFATFrame* triggerFrame = new SecondLevelTriggerVFATFrame(shiftedBuffer);
            ((PositionedVFATFrameCollection*) event->frames)->AddTriggerFrame(triggerFrame);
            offsetWords += triggerFrame->getSize();
        } else {
            WARN("DynamicOptoRxVMEBFile::ProcessOptoRxFrame,Frame header footprint (" << frameHeaderStamp << ") not reckognised - skipping one word from buffer!");
            offsetWords += 1;
        }
    }

    if(rawFramesCount > 0 && clusterizationFramesCount > 0){
        return ProcessMixedRawAndClusterizationVFATFrames(framesRead, event, optoRx);
    } else if(rawFramesCount > 0){
        return ProcessRawVFATFrames(framesRead, event, optoRx);
    } else {
        return ProcessClusterizationVFATFrames(framesRead, event, optoRx);
    }
}

//----------------------------------------------------------------------------------------------------

unsigned int DynamicOptoRxVMEBFile::ProcessRawVFATFrames(vector<PositionedVFATFrame*> framesRead, RawEvent *event, OptoRxSupplementalData* optoRx)
{
    int failedGohCount = 0;
    failedGohCount += FilterOutGohWithoutAllFibers(&framesRead, optoRx);//DAQ sends all frames but empty frames can be corrupted - this check goes first
    failedGohCount += FilterOutInvalidFrames(&framesRead, optoRx);
    failedGohCount += FilterOutInconsitientGoh(&framesRead, optoRx);

    for (unsigned int i = 0; i < framesRead.size(); i++)
      AddFrameToCollection(framesRead[i], optoRx, (PositionedVFATFrameCollection*) event->frames);

    return failedGohCount;
}

//----------------------------------------------------------------------------------------------------

unsigned int DynamicOptoRxVMEBFile::ProcessClusterizationVFATFrames(vector<PositionedVFATFrame*> framesRead,
    RawEvent *event, OptoRxSupplementalData* optoRx)
{
    int failedGohCount = 0;
    failedGohCount += FilterOutInvalidFrames(&framesRead, optoRx);
    failedGohCount += FilterOutInconsitientGoh(&framesRead, optoRx);

    for (unsigned int i = 0; i < framesRead.size(); i++)
        AddFrameToCollection(framesRead[i], optoRx, (PositionedVFATFrameCollection*) event->frames);

    return failedGohCount;
}

//----------------------------------------------------------------------------------------------------

unsigned int DynamicOptoRxVMEBFile::ProcessMixedRawAndClusterizationVFATFrames(vector<PositionedVFATFrame*> framesRead,
    RawEvent *event, OptoRxSupplementalData* optoRx)
{
    ERROR("DynamicOptoRxVMEBFile::ProcessVFATFrames,Clusterization VFAT frame when RAW frames already occured! Checks won't be performed...");
 
    for (unsigned int i = 0; i < framesRead.size(); i++)
        AddFrameToCollection(framesRead[i], optoRx, (PositionedVFATFrameCollection*) event->frames);

    return 0;
}

//----------------------------------------------------------------------------------------------------

unsigned int DynamicOptoRxVMEBFile::FilterOutInvalidFrames(vector<PositionedVFATFrame*>* framesRead, OptoRxSupplementalData* optoRx)
{
    vector<PositionedVFATFrame*>::iterator it = framesRead->begin();
    while (it != framesRead->end())
    {
        PositionedVFATFrame* frame = *it;
        if(!frame->passedHardwareSynchronisationChecks())
        {
            WARN("DynamicOptoRxVMEBFile::FilterOutInvalidFrames,For OptoRx with id " << optoRx->getOptoRxID() << " GOH (with id " << frame->getGohIdx() << ") fiber (with id " << frame->getFiberIdx() << ") - that frame was marked as failing Synchronization Checks");
            it = framesRead->erase(it);
        } else {
            it++;
        }
    }
    return 0; //any goh can fail now (only single frame can) - it will be done later
}

//----------------------------------------------------------------------------------------------------

unsigned int DynamicOptoRxVMEBFile::FilterOutGohWithoutAllFibers(vector<PositionedVFATFrame*>* framesRead, OptoRxSupplementalData* optoRx)
{
    map<pair<VFATFrame::word,VFATFrame::word>, bool> gohAndFiberToFiberPresence;
    set<VFATFrame::word> gohPresent;
    set<VFATFrame::word> invalidGoh;

    // finding present fibers
    for (unsigned int i = 0; i < framesRead->size(); i++)
    {
        PositionedVFATFrame *frame = (*framesRead)[i];
        gohAndFiberToFiberPresence[make_pair(frame->getGohIdx(), frame->getFiberIdx())] = true;
        gohPresent.insert(frame->getGohIdx());
    }

    // finding invalid gohs
    for (set<VFATFrame::word>::iterator it = gohPresent.begin(); it != gohPresent.end(); ++it)
    {
        VFATFrame::word gohId = *it;
        for(int i = 0; i < 16; i++)
        {
            if(!gohAndFiberToFiberPresence[make_pair(gohId, i)])
            {
                WARN("DynamicOptoRxVMEBFile::FilterOutGohWithoutAllFibers,For OptoRx with id " << optoRx->getOptoRxID() << " GOH (with id " << gohId << ") has failed fiber absence test; missing fiber: " << i);
                invalidGoh.insert(gohId);
                break;
            }
        }
    }

    // removing marked as invalid
    vector<PositionedVFATFrame*>::iterator it = framesRead->begin();
    while(it != framesRead->end())
    {
        PositionedVFATFrame *frame = *it;
        if(invalidGoh.find(frame->getGohIdx()) != invalidGoh.end())
        {
            it = framesRead->erase(it);
        } else {
            it++;
        }
    }

    return invalidGoh.size();
}

//----------------------------------------------------------------------------------------------------

unsigned int DynamicOptoRxVMEBFile::FilterOutInconsitientGoh(vector<PositionedVFATFrame*>* framesRead, OptoRxSupplementalData* optoRx)
{
    map<VFATFrame::word, VFATFrame::word> gohToEC;
    map<VFATFrame::word, VFATFrame::word> gohToBC;
    set<VFATFrame::word> invalidGoh;

    // finding invalid frames
    for (unsigned int i = 0; i < framesRead->size(); i++)
    {
        PositionedVFATFrame *frame = (*framesRead)[i];

        VFATFrame::word gohIdx = frame->getGohIdx();
        VFATFrame::word EC = frame->getEC();
        VFATFrame::word BC = frame->getBC();

        if(invalidGoh.find(gohIdx) == invalidGoh.end())
        {
            if (gohToEC.find(gohIdx) == gohToEC.end()) gohToEC[gohIdx] = EC;
            if (gohToBC.find(gohIdx) == gohToBC.end()) gohToBC[gohIdx] = BC;
            if (gohToEC[gohIdx] != EC || gohToBC[gohIdx] != BC) invalidGoh.insert(gohIdx);
        }
    }

    // removing marked as invalid
    vector<PositionedVFATFrame*>::iterator it = framesRead->begin();
    while (it != framesRead->end())
    {
        PositionedVFATFrame *frame = *it;
        if (invalidGoh.find(frame->getGohIdx()) != invalidGoh.end())
        {
            it = framesRead->erase(it);
        } else {
            it++;
        }
    }

    for (set<VFATFrame::word>::iterator it = invalidGoh.begin(); it != invalidGoh.end(); ++it)
        WARN("DynamicOptoRxVMEBFile::FilterOutInconsitientGoh,For OptoRx with id " <<
          optoRx->getOptoRxID() << " GOH (with id " << *it << ") has failed consistency checks");

    return invalidGoh.size();
}

//----------------------------------------------------------------------------------------------------

int32_t DynamicOptoRxVMEBFile::crcLookupTable[256] = {-1};

//todo not tested
void DynamicOptoRxVMEBFile::initializeCrcLookupTable() {
    const uint16_t polynomial = 0x8005; //todo which polynomial? backpatch
    uint16_t remainder;

    /* Compute the remainder of each possible byte dividend. */
    for (int dividend = 0; dividend < 256; ++dividend){

        /* Start with the dividend followed by zeros. */
        remainder = dividend << 8;

        /* Perform modulo-2 division, a bit at a time. */
        for (uint8_t bit = 8; bit > 0; --bit) {
            /* Try to divide the current data bit. */
            if (remainder & (1<<15)){
                remainder = (remainder << 1) ^ polynomial;
            } else {
                remainder = (remainder << 1);
            }
        }

        /*  Store the result into the table. */
        DynamicOptoRxVMEBFile::crcLookupTable[dividend] = remainder;
    }
}

//todo not tested
bool DynamicOptoRxVMEBFile::IsCrcForOptoRxFrameValid(unsigned long long* buffer, unsigned int frameSize)
{
    /* on first use lazy-initialize lookup table */
    if(crcLookupTable[0] < 0) initializeCrcLookupTable();

    /* extracting crc from buffer */
    uint16_t crc = (buffer[2] & 0x00000000FFFF0000llu) >> 16;

    /* creating copy of buffer prepared to be checked */
    uint8_t *withoutCrc = new uint8_t[frameSize];
    memcpy(withoutCrc, buffer, frameSize);

    /* erasing crc from the data */
    withoutCrc[frameSize - 5] = 0;
    withoutCrc[frameSize - 6] = 0;

    uint8_t data;
    uint16_t remainder = 0;

    /* Divide the message by the polynomial, a byte at a time. */
    for (unsigned byte = 0; byte < frameSize; ++byte)
    {
        data = withoutCrc[byte] ^ (remainder >> 8);
        remainder = crcLookupTable[data] ^ (remainder << 8);
    }

    delete [] withoutCrc;

    /* The final remainder is the CRC. */
    return remainder == crc;
}

//----------------------------------------------------------------------------------------------------

void DynamicOptoRxVMEBFile::AddFrameToCollection(PositionedVFATFrame* frame, OptoRxSupplementalData* optoRx, PositionedVFATFrameCollection* frames)
{
    //todo replace SubSystemId, TOTFEDId values (setting 0 is based on OptoRxFrameCollection)
    FramePosition position(0, 0, optoRx->getOptoRxID(), frame->getGohIdx(), frame->getFiberIdx());
    frames->AddNewFrame(position, frame);
}


} // namespace

/****************************************************************************
*
* This is a part of the TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com)
*    
****************************************************************************/

#include "TotemRawDataLibrary/DataFormats/interface/CommonDef.h"
#include "TotemRawDataLibrary/Readers/interface/VMEBFile.h"
#include "TotemRawDataLibrary/Readers/interface/VME2File.h"

// the RunII version
#include "TotemRawDataLibrary/DAQA/interface/event_3_14.h"

//#define DEBUG 0
//#define USE_CASTOR 1

using namespace std;


namespace Totem {

//----------------------------------------------------------------------------------------------------

VMEBFile::VMEBFile() : BaseVMEBFile()
{
#ifdef DEBUG
  printf(">> VMEBFile::VMEBFile \n");
#endif
}

//----------------------------------------------------------------------------------------------------

VMEBFile::~VMEBFile()
{
#ifdef DEBUG
  printf(">> VMEBFile::~VMEBFile, this = %p, dataPtr = %p\n", (void *) this, dataPtr);
#endif
}

//----------------------------------------------------------------------------------------------------

unsigned char VMEBFile::GetNextEvent(RawEvent *event)
{
  int retCode = BaseVMEBFile::GetNextEvent(event);
  if(retCode) return retCode;

  // process the buffer
  OptoRxVFATFrameCollection *oc = (OptoRxVFATFrameCollection *) event->frames;
  oc->Invalidate();
  unsigned int errorCounter = ProcessVMEBEvent(dataPtr, oc, event);
  
  if (errorCounter > 0) {
    ERROR("VMEBFile::GetNextEvent") << errorCounter << " GOH blocks have failed consistency checks in event "
      << event->dataEventNumber << "." << c_endl;
    corruptedEventCounter++;
  }

  return 0;
}

//----------------------------------------------------------------------------------------------------

unsigned int VMEBFile::ProcessVMEBEvent(char *ptr, OptoRxVFATFrameCollection *oc, RawEvent *event)
{
  eventHeaderStruct *eventHeader = (eventHeaderStruct *) ptr;
  bool superEvent = TEST_ANY_ATTRIBUTE(eventHeader->eventTypeAttribute, ATTR_SUPER_EVENT);

#ifdef DEBUG
  printf(">> VMEBFile::ProcessVMEBEvent\n");

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

      errorCounter += ProcessSubEvent(ptr + offset, oc, event);

      offset += subEvSize;
#ifdef DEBUG 
      printf("\t> offset after %i\n", offset);
#endif
    }
  } else
    errorCounter += ProcessSubEvent(ptr, oc, event);

  return errorCounter;
}

//----------------------------------------------------------------------------------------------------

unsigned int VMEBFile::ProcessSubEvent(char *ptr, OptoRxVFATFrameCollection *oc, RawEvent *event)
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
        case etOptoRxOld:
        case etOptoRx:
          errorCounter += VME2File::ProcessOptoRxFrame(payloadPtr, payloadSize, oc, event);
          break;

        case etLoneG:
          ProcessLoneGFrame(payloadPtr, payloadSize, event);
          break;

        default:
          ERROR("VMEBFile::ProcessSubEvent") << "Unknown equipment type: " << eq->equipmentType << ". Skipping." << c_endl;
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


} // namespace

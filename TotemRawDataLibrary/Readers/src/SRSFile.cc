#include "TotemRawDataLibrary/DataFormats/interface/CommonDef.h"
#include "TotemRawDataLibrary/Readers/interface/SRSFile.h"
//#include "TotemRawDataLibrary/Readers/interface/VME2File.h"

// the RunII version
#include "TotemRawDataLibrary/DAQA/interface/event_3_14.h"

//#define DEBUG 0

using namespace std;


namespace Totem {

const unsigned int SRSFile::eventHeaderSize = sizeof(eventHeaderStruct);


//----------------------------------------------------------------------------------------------------

SRSFile::SRSFile() : dataPtr(NULL), dataPtrSize(0), infile(NULL)
{
#ifdef DEBUG
  printf(">> SRSFile::SRSFile \n");
#endif
}

//----------------------------------------------------------------------------------------------------

SRSFile::~SRSFile()
{
#ifdef DEBUG
  printf(">> SRSFile::~SRSFile, this = %p, dataPtr = %p\n", (void *) this, dataPtr);
#endif

  Close();

  if(dataPtr != NULL)
    delete [] dataPtr;
}

//----------------------------------------------------------------------------------------------------

DataFile::OpenStatus SRSFile::Open(const std::string &fn)
{
  // check if the source is a SRS File
  size_t dotPos = fn.rfind('.');
  string extension = (dotPos == string::npos) ? "" : fn.substr(dotPos);
  if (extension.compare(".srs") != 0) //todo check!
    return osWrongFormat;

  infile = StorageFile::CreateInstance(fn);
  infile->OpenFile();

  if (!infile->IsOpened()) {
    infile->PrintError("Error while opening file in SRSFile::Open");
    return osCannotOpen;
  }

  // set maximum of data counters
  dataEventNumberMax = 0xFFFFFFFF;
  dataConfNumberMax = 0;

  // reset counters
  indexStatus = isNotIndexed;
  corruptedEventCounter = 0;

  filename =  fn;

  return osOK;
}

//----------------------------------------------------------------------------------------------------

DataFile::OpenStatus SRSFile::Open(StorageFile *storageFile)
{
  std::string fn = storageFile->GetURLPath();

  // check if the source is a SRS File
  size_t dotPos = fn.rfind('.');
  string extension = (dotPos == string::npos) ? "" : fn.substr(dotPos);
  if (extension.compare(".srs") != 0) //todo check!
    return osWrongFormat;

  infile = storageFile;
  infile->OpenFile();

  if (!infile->IsOpened()) {
    infile->PrintError("Error while opening file in SRSFile::Open");
    return osCannotOpen;
  }

  // set maximum of data counters
  dataEventNumberMax = 0xFFFFFFFF;
  dataConfNumberMax = 0;

  // reset counters
  indexStatus = isNotIndexed;
  corruptedEventCounter = 0;

  filename =  fn;

  return osOK;
}

//----------------------------------------------------------------------------------------------------

void SRSFile::Close()
{
#ifdef DEBUG
  printf(">> SRSFile::Close, this = %p", (void*) this);
#endif
  
  if (infile && infile->CloseFile() == EOF)
  {
    ERROR("SRSFile::Close") << "Cannot close the file." << c_endl;
    delete infile;
  } else {
	  infile=NULL;
  }
}

//----------------------------------------------------------------------------------------------------

unsigned char SRSFile::ReadToBuffer(unsigned int bytesToRead, unsigned int offset)
{
#ifdef DEBUG
  printf(">> SRSFile::ReadToBuffer(%u, %u), this = %p, dataPtr = %p, dataPtrSize = %u\n",
      bytesToRead, offset, (void*) this, dataPtr, dataPtrSize);
#endif

  // allocate new memory block if current one is too small
  if(dataPtrSize<bytesToRead+offset) {
    char *newPtr = NULL;
    try {
      newPtr = new char[bytesToRead+offset];
    }
    catch (bad_alloc& ba) {
      ERROR("SRSFile::ReadToBuffer") << "Cannot allocate buffer large enough." << c_endl;
      return 2;
    }

    if (dataPtr) {
      if (offset > 0)
        memcpy(newPtr, dataPtr, offset);
      delete [] dataPtr;
    }

    dataPtr = newPtr;
    dataPtrSize = bytesToRead+offset;
  }

  // read data at given offset
  unsigned int bytesRead;
  
  int eofFlag;
  bytesRead = infile->ReadData(dataPtr + offset, sizeof(char), bytesToRead);
  eofFlag = infile->CheckEOF();

  if (bytesRead != bytesToRead && !(bytesRead == 0 && eofFlag)) 
  {
    ERROR("SRSFile::ReadToBuffer") << "Reading from file to buffer failed. Only " << bytesRead
  			      << " B read from " << bytesToRead << " B." << c_endl;
    return 1;
  }

  return 0;
}
//----------------------------------------------------------------------------------------------------

unsigned char SRSFile::GetNextEvent(RawEvent* event)
{
#ifdef DEBUG
  printf(">> SRSFile::GetNextEvent, this = %p\n", (void*)this);
  printf("\teventHeaderSize = %u\n", eventHeaderSize);
#endif

  eventHeaderStruct *eventHeader = NULL;

  long int currentPos = infile->CurrentPosition();
  while (!infile->CheckEOF())
  {
    // read next header
    if (ReadToBuffer(eventHeaderSize, 0) != 0)
      return 10;

    eventHeader = (eventHeaderStruct *) dataPtr;

    // check the sanity of header data
    if (eventHeader->eventMagic != EVENT_MAGIC_NUMBER) {
      ERROR("SRSFile::GetNextEvent") << "Event magic check failed (" << hex << eventHeader->eventMagic << "!=" << EVENT_MAGIC_NUMBER << dec << "). Exiting." << c_endl;
      return 1;
    }

    unsigned int N = eventHeader->eventSize;
    if (N<eventHeaderSize) {
      ERROR("SRSFile::GetNextEvent") << "Event size (" << N << ") smaller than header size (" << eventHeaderSize << "). Exiting." << c_endl;
      return 1;
    }

    // get next event from the file (the header has already been read)
    if (ReadToBuffer(N-eventHeaderSize, eventHeaderSize) != 0)
      return 10;

    // because dataPtr could move, we have to be sure, that eventHeader still points the same as dataPtr
    eventHeader = (eventHeaderStruct *) dataPtr;

    // skip non physics events
    if (eventHeader->eventType != PHYSICS_EVENT)
      continue;

    break;
  }

  // check if the end of the file has been reached
  int eofFlag = infile->CheckEOF();

  if (eofFlag) {
    // if indexing, set indexed
    if (indexStatus == isIndexing)
      indexStatus = isIndexed;
    return 1;
  }

  // process the buffer
  OptoRxVFATFrameCollection *oc = (OptoRxVFATFrameCollection *) event->frames;
  oc->Invalidate();
  unsigned int errorCounter = ProcessSRSEvent(dataPtr, oc, event);

#ifdef DEBUG
  printf("* %u, %u, %u\n",
    EVENT_ID_GET_NB_IN_RUN( eventHeader->eventId ),
    EVENT_ID_GET_BURST_NB( eventHeader->eventId ),
    EVENT_ID_GET_NB_IN_BURST( eventHeader->eventId )
    );
#endif

  if (errorCounter > 0)
  {
    ERROR("SRSFile::GetNextEvent") << errorCounter << " GOH blocks have failed consistency checks in event "
      << event->dataEventNumber << "." << c_endl;
    corruptedEventCounter++;
  }

  // save event position
  if (indexStatus == isIndexing)
    positions.push_back(currentPos);

  return 0;
}

//----------------------------------------------------------------------------------------------------

unsigned int SRSFile::ProcessSRSEvent(char *ptr, OptoRxVFATFrameCollection *oc, RawEvent *event)
{
  eventHeaderStruct *eventHeader = (eventHeaderStruct *) ptr;
  bool superEvent = TEST_ANY_ATTRIBUTE(eventHeader->eventTypeAttribute, ATTR_SUPER_EVENT);

#ifdef DEBUG
  printf(">> SRSFile::ProcessVMEBEvent\n");

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

unsigned int SRSFile::ProcessSubEvent(char *ptr, OptoRxVFATFrameCollection *oc, RawEvent *event)
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
    unsigned long long *payloadPtr = (unsigned long long *)(ptr + offset + equipmentHeaderStructSize + 4); //TODO CHECK - because 32 bits - 0xFAFAFAFA

#ifdef DEBUG 
    printf("\t\t\tequipmentSize = %u\n", eq->equipmentSize);
    printf("\t\t\tequipmentType = %u\n", eq->equipmentType);
    printf("\t\t\tequipmentId = %u\n", eq->equipmentId);
    printf("\t\t\tequipmentTypeAttribute = %p\n", (void*) eq->equipmentTypeAttribute);
    printf("\t\t\tequipmentBasicElementSize = %u\n", eq->equipmentBasicElementSize);
    printf("\t\t\t\t\tpayload size = %li\n", payloadSize);
    printf("\t\t\t\t\tpayload ptr = %p\n", (void*) payloadPtr);
#endif

    if (payloadSize > 0)
    {
      payloadSize /= 8; // bytes -> words

      switch (eq->equipmentType)
      {
        case etOptoRx:
          errorCounter += SRSFile::ProcessOptoRxFrame(payloadPtr, payloadSize, oc, event);
          break;

        default:
          ERROR("SRSFile::ProcessSubEvent") << "Unknown equipment type: " << eq->equipmentType << ". Skipping." << c_endl;
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

unsigned char SRSFile::GetEvent(unsigned long n, RawEvent *event)
{
  if (n >= positions.size()) {
    ERROR("SRSFile::GetEvent") << "Requested event number (" << n
      << ") is larger than event count (" << positions.size() << ")." << c_endl;
    return 1;  
  }

  if (!infile) {
    ERROR("SRSFile::GetEvent") << "No file open." << c_endl;
    return 1;
  }

  if (infile->Seek(positions[n])) {
    ERROR("SRSFile::GetEvent") << "Seek to pos " << positions[n] << " unsuccessful (ftell="
      << infile->CurrentPosition() << ")." << c_endl;
    return 1;
  }

  return GetNextEvent(event);
}

//----------------------------------------------------------------------------------------------------

void SRSFile::Rewind()
{
  infile->Seek(0);

  corruptedEventCounter = 0;
  positions.clear();
  indexStatus = isNotIndexed;
}

//----------------------------------------------------------------------------------------------------

unsigned int SRSFile::ProcessOptoRxFrame(SRSFile::word *buf, unsigned int frameSize,
  OptoRxVFATFrameCollection *fc, RawEvent *event)
{
    // get OptoRx metadata
    unsigned long long head = buf[0];
    unsigned int OptoRxId = (head >> 8) & 0xFFF;
    unsigned long BX = (head >> 20) & 0xFFF;
    unsigned long LV1 = (head >> 32) & 0xFFFFFF;

    // save metadata to event
    if (event)
    {
        OptoRxMetaData &md = event->optoRxMetaData[OptoRxId];
        md.BX = BX;
        md.LV1 = LV1;
    }

    // TODO set appropriate value for LoneG OptoRxId
    if (OptoRxId == 0x29c)
    {
        return ProcessLoneGFrame(buf + 2, frameSize - 4, event);
    }

    // process standard (VFAT) OptoRx data
    OptoRxVFATFrameCollection::IteratorType block = fc->InsertOptoRxBlock(OptoRxId);
    unsigned int subFrames = (frameSize - 2) / 194;

    #ifdef DEBUG
      printf(">> SRSFile::ProcessOptoRxFrame > OptoRxId = %u, BX = %lu, LV1 = %lu, frameSize = %u, subFrames = %u)\n",
        OptoRxId, BX, LV1, frameSize, subFrames);
    #endif

    unsigned int errorCounter = 0;

    // process all sub-frames
    for (unsigned int r = 0; r < subFrames; ++r)
    {
        for (unsigned int c = 0; c < 4; ++c)
        {
            unsigned int head = (buf[1 + 194 * r] >> (16 * c)) & 0xFFFF;
            unsigned int foot = (buf[194 + 194 * r] >> (16 * c)) & 0xFFFF;

            #ifdef DEBUG
                printf(">>>> r = %i, c = %i: S = %i, BOF = %i, EOF = %i, ID = %i, ID' = %i\n", r, c, head & 0x1, head >> 12, foot >> 12, (head >> 8) & 0xF, (foot >> 8) & 0xF);
            #endif

            // stop if this GOH is NOT active
            if ((head & 0x1) == 0)
                continue;

            #ifdef DEBUG
                printf("\tHeader active (%04x -> %x).\n", head, head & 0x1);
            #endif

            // check structure
            if (head >> 12 != 0x4 || foot >> 12 != 0xB || ((head >> 8) & 0xF) != ((foot >> 8) & 0xF))
            {
                //stringstream ss;
                char ss[200];
                if (head >> 12 != 0x4) sprintf(ss, "\n\tHeader is not 0x4 as expected (%x).", head);
                //ss << "\n\tHeader is not 0x4 as expected (" << hex << head << dec << ").";
                if (foot >> 12 != 0xB) sprintf(ss, "\n\tFooter is not 0xB as expected (%x).", foot);
                //ss << "\n\tFooter is not 0xB as expected (" << hex << foot << dec << ").";
                if (((head >> 8) & 0xF) != ((foot >> 8) & 0xF))
                    sprintf(ss, "\n\tIncompatible GOH IDs in header (%x) and footer (%x).", ((head >> 8) & 0xF),
                            ((foot >> 8) & 0xF));
                //ss << "\n\tIncompatible GOH IDs in header ("
                //<< hex << ((head >> 8) & 0xF) << ") and footer ("
                //<< ((foot >> 8) & 0xF) << dec << ").";

                ERROR("SRSFile::ProcessOptoRxFrame") << "Wrong payload structure (in GOH block row " << r <<
                " and column " << c
                << ") in OptoRx frame ID " << OptoRxId << ". GOH block omitted." << ss << c_endl;

                errorCounter++;
                continue;
            }

            unsigned int goh = (head >> 8) & 0xF;
            VFATFrame::word **dataPtrs = fc->ValidateGOHBlock(block, goh);

            #ifdef DEBUG
                printf(">>>> transposing GOH block at prefix: %i, dataPtrs = %p\n", OptoRxId*192 + goh*16, dataPtrs);
            #endif

            // deserialization
            for (int i = 0; i < 192; i++)
            {
                int iword = 11 - i / 16;  // number of current word (11...0)
                int ibit = 15 - i % 16;   // number of current bit (15...0)
                unsigned int w = (buf[i + 2 + 194 * r] >> (16 * c)) & 0xFFFF;

                // Fill the current bit of the current word of all VFAT frames
                for (int idx = 0; idx < 16; idx++)
                    if (w & (1 << idx))
                        dataPtrs[idx][iword] |= (1 << ibit);
            }
        }
    }

    return errorCounter;
}

//----------------------------------------------------------------------------------------------------

unsigned int SRSFile::ProcessLoneGFrame(SRSFile::word *oBuf, signed long size, RawEvent *ev)
{
  if (size != 20)
  {
    ERROR("SRSFile::ProcessLoneGFrame") << "Wrong LoneG frame size: " << size << " (shall be 20)." << c_endl;
    return 1;
  }

  if (!ev)
    return 0;

  // buffer mapping: OptoRx buffer --> LoneG buffer
  SRSFile::word buf[5];
  for (unsigned int i = 0; i < 5; i++)
    buf[i] = 0;

  for (unsigned int i = 0; i < 20; i++)
  {
      int row = i / 4;
      int col = i % 4;
      buf[row] |= (oBuf[i] & 0xFFFF) << (col * 16);
  }

  ev->triggerData.type = (buf[0] >> 56) & 0xF;
  ev->triggerData.event_num = (buf[0] >> 32) & 0xFFFFFF;
  ev->triggerData.bunch_num = (buf[0] >> 20) & 0xFFF;
  ev->triggerData.src_id = (buf[0] >> 8) & 0xFFF;

  ev->triggerData.orbit_num = (buf[1] >> 32) & 0xFFFFFFFF;
  ev->triggerData.revision_num = (buf[1] >> 24) & 0xFF;

  ev->triggerData.run_num = (buf[2] >> 32) & 0xFFFFFFFF;
  ev->triggerData.trigger_num = (buf[2] >> 0) & 0xFFFFFFFF;

  ev->triggerData.inhibited_triggers_num = (buf[3] >> 32) & 0xFFFFFFFF;
  ev->triggerData.input_status_bits = (buf[3] >> 0) & 0xFFFFFFFF;

#ifdef DEBUG
  printf(">> SRSFile::ProcessLoneGFrame > size = %li\n", size);
  printf("\ttype = %x, event number = %x, bunch number = %x, id = %x\n",
    ev->triggerData.type, ev->triggerData.event_num, ev->triggerData.bunch_num, ev->triggerData.src_id);
  printf("\torbit number = %x, revision = %x\n",
    ev->triggerData.orbit_num, ev->triggerData.revision_num);
  printf("\trun number = %x, trigger number = %x\n",
    ev->triggerData.run_num, ev->triggerData.trigger_num);
  printf("\tinhibited triggers = %x, input status bits = %x\n",
    ev->triggerData.inhibited_triggers_num, ev->triggerData.input_status_bits);
#endif

  return 0;
}


} // namespace

/****************************************************************************
*
* This is a part of the TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com)
*    
****************************************************************************/

#include "TotemRawDataLibrary/DataFormats/interface/CommonDef.h"
#include "TotemRawDataLibrary/Readers/interface/BaseVMEBFile.h"
#include "TotemRawDataLibrary/Readers/interface/VME2File.h"

// the RunII version
#include "TotemRawDataLibrary/DAQA/interface/event_3_14.h"

//#define DEBUG 0
//#define USE_CASTOR 1

using namespace std;


namespace Totem {

const unsigned int BaseVMEBFile::eventHeaderSize = sizeof(eventHeaderStruct);

unsigned char BaseVMEBFile::GetNextEvent(RawEvent *event)
{
#ifdef DEBUG
  printf(">> VMEBFile::GetNextEvent, this = %p\n", (void*)this);
  printf("\teventHeaderSize = %u\n", eventHeaderSize);
#endif

  eventHeaderStruct *eventHeader = NULL;

#ifdef USE_CASTOR
  long int currentPos = rfio_ftell(infile);
  while (!rfio_feof(infile))
  {
#else
  long int currentPos = ftell(infile);
  while (!feof(infile))
  {
#endif

    // read next header
    if (ReadToBuffer(eventHeaderSize, 0) != 0)
      return 10;

    eventHeader = (eventHeaderStruct *) dataPtr;

    // check the sanity of header data
    if (eventHeader->eventMagic != EVENT_MAGIC_NUMBER) {
      ERROR("VMEBFile::GetNextEvent") << "Event magic check failed (" << hex << eventHeader->eventMagic << "!=" << EVENT_MAGIC_NUMBER << dec << "). Exiting." << c_endl;
      return 1;
    }

    unsigned int N = eventHeader->eventSize;
    if (N<eventHeaderSize) {
      ERROR("VMEBFile::GetNextEvent") << "Event size (" << N << ") smaller than header size (" << eventHeaderSize << "). Exiting." << c_endl;
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
#ifdef USE_CASTOR
  int eofFlag = rfio_feof(infile);
#else
  int eofFlag = feof(infile);
#endif

  if (eofFlag) {
    // if indexing, set indexed
    if (indexStatus == isIndexing)
      indexStatus = isIndexed;
    return 1;
  }

#ifdef DEBUG
  printf("* %u, %u, %u\n",
    EVENT_ID_GET_NB_IN_RUN( eventHeader->eventId ),
    EVENT_ID_GET_BURST_NB( eventHeader->eventId ),
    EVENT_ID_GET_NB_IN_BURST( eventHeader->eventId )
    );
#endif

  // save event position
  if (indexStatus == isIndexing)
    positions.push_back(currentPos);

  return 0;
}

//----------------------------------------------------------------------------------------------------

BaseVMEBFile::BaseVMEBFile() : dataPtr(NULL), dataPtrSize(0), infile(NULL)
{
#ifdef DEBUG
  printf(">> BaseVMEBFile::BaseVMEBFile \n");
#endif
}

//----------------------------------------------------------------------------------------------------

BaseVMEBFile::~BaseVMEBFile()
{
#ifdef DEBUG
  printf(">> BaseVMEBFile::~BaseVMEBFile, this = %p, dataPtr = %p\n", (void *) this, dataPtr);
#endif

  Close();

  if(dataPtr != NULL)
    delete [] dataPtr;
}

//----------------------------------------------------------------------------------------------------

DataFile::OpenStatus BaseVMEBFile::Open(const std::string &fn)
{
  // check if the source is a VMEB File
  size_t dotPos = fn.rfind('.');
  string extension = (dotPos == string::npos) ? "" : fn.substr(dotPos);
  if (extension.compare(getExtension()) != 0)
    return osWrongFormat;

#ifdef USE_CASTOR
  infile = rfio_fopen64(const_cast<char*>(fn.c_str()),const_cast<char*>("r"));
#else
  infile = fopen64(const_cast<char*>(fn.c_str()),const_cast<char*>("r"));
#endif
  
  if (!infile)
  {
#ifdef USE_CASTOR
	  rfio_perror(const_cast<char*>("Error while opening file in BaseVMEBFile::Open"));
#endif
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

void BaseVMEBFile::Close()
{
#ifdef DEBUG
  printf(">> BaseVMEBFile::Close, this = %p", (void*) this);
#endif
  
  if (infile) {
#ifdef USE_CASTOR    
    if (rfio_fclose(infile)==EOF)
#else
    if (fclose(infile)==EOF)
#endif
      ERROR("BaseVMEBFile::Close") << "Cannot close the file." << c_endl;
    else infile=NULL;
  }
}

//----------------------------------------------------------------------------------------------------

unsigned char BaseVMEBFile::ReadToBuffer(unsigned int bytesToRead, unsigned int offset)
{
#ifdef DEBUG
  printf(">> BaseVMEBFile::ReadToBuffer(%u, %u), this = %p, dataPtr = %p, dataPtrSize = %u\n", 
      bytesToRead, offset, (void*) this, dataPtr, dataPtrSize);
#endif

  // allocate new memory block if current one is too small
  if(dataPtrSize<bytesToRead+offset) {
    char *newPtr = NULL;
    try {
      newPtr = new char[bytesToRead+offset];
    }
    catch (bad_alloc& ba) {
      ERROR("BaseVMEBFile::ReadToBuffer") << "Cannot allocate buffer large enough." << c_endl;
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
#ifdef USE_CASTOR  
  bytesRead = rfio_fread(dataPtr+offset, sizeof(char), bytesToRead, infile);
  eofFlag = rfio_feof(infile);
#else
  bytesRead = fread(dataPtr+offset, sizeof(char), bytesToRead, infile);
  eofFlag = feof(infile);
#endif

  if (bytesRead != bytesToRead && !(bytesRead == 0 && eofFlag)) 
  {
    ERROR("BaseVMEBFile::ReadToBuffer") << "Reading from file to buffer failed. Only " << bytesRead
  			      << " B read from " << bytesToRead << " B." << c_endl;
    return 1;
  }

  return 0;
}

//----------------------------------------------------------------------------------------------------

unsigned char BaseVMEBFile::GetEvent(unsigned long n, RawEvent *event)
{
  if (n >= positions.size()) {
    ERROR("BaseVMEBFile::GetEvent") << "Requested event number (" << n
      << ") is larger than event count (" << positions.size() << ")." << c_endl;
    return 1;  
  }

  if (!infile) {
    ERROR("BaseVMEBFile::GetEvent") << "No file open." << c_endl;
    return 1;
  }

#ifdef USE_CASTOR
  // fixing bug in rfio_api
  if (rfio_feof(infile))
    if (Reopen())
      return 1;

  if (rfio_fseek(infile, positions[n], SEEK_SET)) {
    ERROR("BaseVMEBFile::GetEvent") << "Seek to pos " << positions[n] << " unsuccessful (rfio_ftell="
      << rfio_ftell(infile) << ")." << c_endl;
    return 1;
  }
#else
  if (fseek(infile, positions[n], SEEK_SET)) {
    ERROR("BaseVMEBFile::GetEvent") << "Seek to pos " << positions[n] << " unsuccessful (ftell="
      << ftell(infile) << ")." << c_endl;
    return 1;
  }
#endif  

  return GetNextEvent(event);
}

//----------------------------------------------------------------------------------------------------

void BaseVMEBFile::Rewind()
{
#ifdef USE_CASTOR
  // fixing bug in rfio_api
  if (rfio_feof(infile))
    Reopen();
  else 
    rfio_fseek(infile, 0, SEEK_SET);
#else
  fseek(infile, 0, SEEK_SET);
#endif  

  corruptedEventCounter = 0;
  positions.clear();
  indexStatus = isNotIndexed;
}

//----------------------------------------------------------------------------------------------------

#ifdef USE_CASTOR
unsigned int BaseVMEBFile::Reopen()
{
  rfio_fclose(infile);
  infile = rfio_fopen64((char *)filename.c_str(), const_cast<char*>("r"));
  if (!infile) {
    ERROR("BaseVMEBFile::Reopen") << "Could not reopen file `" << filename.c_str() << "'." << c_endl;
    return 1;
  }

  return 0;
}
#endif

//----------------------------------------------------------------------------------------------------

unsigned int BaseVMEBFile::ProcessLoneGFrame(BaseVMEBFile::word *buf, signed long size, RawEvent *ev)
{
  if (size != 5) {
    ERROR("BaseVMEBFile::ProcessLoneGFrame") << "Wrong LoneG frame size: " << size << " (shall be 5)." << c_endl;
    return 1;
  }

  if (!ev)
    return 0;

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
  printf(">> BaseVMEBFile::ProcessLoneGFrame > size = %li\n", size);
  printf("\ttype = %i, event number = %u, bunch number = %u, id = %u\n",
    ev->triggerData.type, ev->triggerData.event_num, ev->triggerData.bunch_num, ev->triggerData.src_id);
  printf("\torbit number = %u, revision = %u\n",
    ev->triggerData.orbit_num, ev->triggerData.revision_num);
  printf("\trun number = %u, trigger number = %u\n",
    ev->triggerData.run_num, ev->triggerData.trigger_num);
  printf("\tinhibited triggers = %u, input status bits = %u\n",
    ev->triggerData.inhibited_triggers_num, ev->triggerData.input_status_bits);
#endif

  return 0;
}


} // namespace

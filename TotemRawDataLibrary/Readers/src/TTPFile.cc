/****************************************************************************
*
* This is a part of theTOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors: 
*   Jan Kaspar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "TotemRawDataLibrary/DataFormats/interface/CommonDef.h"
#include "TotemRawDataLibrary/Readers/interface/TTPFile.h"

#include <cstdlib>

#ifdef USE_CASTOR
  #include "shift.h"
  // undefine these below, to be able to include standard c++ library headers.
  //#undef min
  //#undef log
#endif

namespace Totem {

TTPFile::TTPFile()
{
  indexStatus = isNotIndexed;
}

//----------------------------------------------------------------------------------------------------

TTPFile::~TTPFile()
{
  Close();
}

//----------------------------------------------------------------------------------------------------

DataFile::OpenStatus TTPFile::Open(const std::string &fn)
{
  // try to open
#ifdef USE_CASTOR
  file = rfio_fopen64(const_cast<char*>(fn.c_str()), const_cast<char*>("rb"));
#else
  file = fopen64(fn.c_str(), "rb");
#endif

  if (!file)
    return osCannotOpen;

  // check header
  char header[7]; header[6] = 0;
#ifdef USE_CASTOR
  rfio_fread(header, sizeof(char), 6, file);
#else
  fread(header, sizeof(char), 6, file);
#endif
  if (strcmp(header, "TTP001"))
    return osWrongFormat;
  
  this->filename = (char*)malloc(sizeof(char)*fn.size()+1);
  strcpy(this->filename, fn.c_str());

  return osOK;
}

//----------------------------------------------------------------------------------------------------

void TTPFile::Close()
{
  if (file) {
#ifdef USE_CASTOR
    rfio_fclose(file);
#else
    fclose(file);    
#endif
    file = NULL;
  }
}

//----------------------------------------------------------------------------------------------------

unsigned char TTPFile::GetNextEvent(RawEvent *event)
{
  // remember file position
#ifdef USE_CASTOR
  int pos = rfio_ftell(file);
#else
  int pos = ftell(file);
#endif

  unsigned char r = ReadOneEvent(event);
  if (r > 0) {
    if (indexStatus == isIndexing)
      indexStatus = isIndexed;

    if (r == 2)
      WARN("TTPFile::GetNextEvent") << "File not terminated correctly." << c_endl;

    return 1;
  }

  // add to list of indices
  if (indexStatus == isIndexing)
    positions.push_back(pos);

  return 0;
}

//----------------------------------------------------------------------------------------------------

/**
Reads one event from current position of the stream. Return values are
  - 0 => OK
  - 1 => EOF at correct place
  - 2 => EOF at wrong place
  - 3 => Forbidden number of events
**/
unsigned char TTPFile::ReadOneEvent(RawEvent *event)
{
  // read frame count
  unsigned short nFrames;  
#ifdef USE_CASTOR
  rfio_fread(&nFrames,sizeof(char), sizeof(nFrames),file);
#else
  fread(&nFrames,sizeof(char), sizeof(nFrames),file);
#endif

  int flagEof;
  int flagError;
#ifdef USE_CASTOR
  flagEof = rfio_feof(file);
  flagError = rfio_ferror(file);
#else
  flagEof = feof(file);
  flagError = ferror(file);
#endif

  if (flagEof || flagError)
    return 1;

  /*
  if (nFrames > 4) {
    fprintf(stderr, ERROR(TTPFile::ReadOneEvent,Number of frames to be read %i. It cannot be greater that 4. Exiting.), nFrames);
    return 3;
  }
  */

  SimpleVFATFrameCollection *sc = (SimpleVFATFrameCollection *) event->frames;
  sc->Clear();

  // read frames to event
  for (unsigned short i = 0; i < nFrames; i++) {
    OldVFATFrame *fp = sc->InsertEmptyFrame(i);
#ifdef USE_CASTOR
    rfio_fread(const_cast<short unsigned int*>(fp->getData()), sizeof(char), sizeof(VFATFrame), file);
#else
    fread(const_cast<short unsigned int*>(fp->getData()), sizeof(char), sizeof(VFATFrame), file);
#endif
  }

#ifdef USE_CASTOR
  flagEof = rfio_feof(file);
  flagError = rfio_ferror(file);
#else
  flagEof = feof(file);
  flagError = ferror(file);
#endif
  return (flagEof || flagError) ? 2 : 0;
}

//----------------------------------------------------------------------------------------------------

unsigned char TTPFile::GetEvent(unsigned long idx, RawEvent *event)
{
  if (indexStatus != isIndexed || idx >= positions.size())
    return 1;

  // this clears error flags, necessary to to make any operations with the file

  Myfseek(file, positions[idx], SEEK_SET);
  ReadOneEvent(event);
  return 0;
}

//----------------------------------------------------------------------------------------------------

void TTPFile::Rewind()
{
  // clear file's error flags and seek to the beginning
  Myfseek(file, 6, SEEK_SET);

  // reset status
  indexStatus = isNotIndexed;
  positions.clear();
}

//----------------------------------------------------------------------------------------------------

void TTPFile::StartIndexing()
{
  indexStatus = isIndexing;
  positions.clear();
}

//----------------------------------------------------------------------------------------------------

// myfseek function, eliminating bug with shift library
int TTPFile::Myfseek(FILE* stream, long int offset, int origin) {
#ifdef USE_CASTOR
  if (rfio_feof(stream)) 
    Reopen();  
  else rfio_fseek(stream, offset, origin);
#else
  return fseek(stream, offset, origin);
#endif
  return 1;
}

//----------------------------------------------------------------------------------------------------

unsigned int TTPFile::Reopen()
{
#ifdef USE_CASTOR
  rfio_fclose(file);
  file = rfio_fopen(filename, const_cast<char*>("rb"));
#else
  fclose(file);
  file = fopen(filename, "rb");
#endif
  if (!file) {
      ERROR("TPPFile::Reopen") << "Could not reopen file `" << filename << "'." << c_endl;
    return 1;
  }

  return 0;
}

} // namespace


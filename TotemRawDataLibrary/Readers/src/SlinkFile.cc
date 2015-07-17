/****************************************************************************
*
* This is a part of the TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "TotemRawDataLibrary/DataFormats/interface/CommonDef.h"
#include "TotemRawDataLibrary/Readers/interface/SlinkFile.h"

#include <cstdlib>

#ifdef USE_CASTOR
  #include "shift.h"
  // uncomment these below to make it work with c++ streams
  //#undef min
  //#undef log
#endif

using namespace std;

namespace Totem {

//----------------------------------------------------------------------------------------------------

  SlinkFile::SlinkFile() : filename(NULL)
{
  file = NULL;

  // default values for BOF, EOF
  seqBOF.clear(); seqEOF.clear();
  seqBOF.push_back(sequence(0x50c800c700c600c5ULL)); seqEOF.push_back(sequence(0xa0c800c700c600c5ULL));
  seqBOF.push_back(sequence(0x5000000000000b13ULL)); seqEOF.push_back(sequence(0xa000000000000b13ULL));
  seqBOF.push_back(sequence(0x5000000000000b14ULL)); seqEOF.push_back(sequence(0xa000000000000b14ULL));
  seqBOF.push_back(sequence(0x5000000000000000ULL)); seqEOF.push_back(sequence(0xa000000000000000ULL, 0xffffffff00000000ULL));
}

//----------------------------------------------------------------------------------------------------

SlinkFile::~SlinkFile()
{
  if (file) Close();
}

//----------------------------------------------------------------------------------------------------

DataFile::OpenStatus SlinkFile::Open(const std::string &filename)
{
  // look for .txt extension and set type
  size_t dotPos = filename.rfind('.');
  string extension = (dotPos == string::npos) ? "" : filename.substr(dotPos);
  if (extension.compare(".txt") == 0)
    type = tAscii;
  else
    type = tBinary;

  // open the file
  if (type == tAscii)
    return osWrongFormat;
  else
#ifdef USE_CASTOR
    file = rfio_fopen64(const_cast<char*>(filename.c_str()), const_cast<char*>("rb"));
#else
    file = fopen64(filename.c_str(), "rb");
#endif

  if (!file)
    return osCannotOpen;

  // reset counters
  indexStatus = isNotIndexed;
  corrNum = 0;
  
  this->filename = (char*)malloc(sizeof(char)*filename.size()+1);
  strcpy(this->filename, filename.c_str());

  return osOK;
}

//----------------------------------------------------------------------------------------------------

void SlinkFile::Close()
{
  if (file) 
#ifdef USE_CASTOR
    rfio_fclose(file);
#else
  fclose(file);
#endif

  file = NULL;

  if (filename) {
    delete filename;
    filename = NULL;
  }
}

//----------------------------------------------------------------------------------------------------

char SlinkFile::Read(unsigned long long &n)
{

#ifdef USE_CASTOR
  int ret = rfio_fread(&n, 1, 8, file);
#else
  int ret = fread(&n, 1, 8, file);
#endif

  if (ret != 8) {     
#ifdef USE_CASTOR
    int eofFlag = rfio_feof(file);
#else
    int eofFlag = feof(file);
#endif
    if (!eofFlag)
      ERROR("SlinkFile::Read") << "Problems with reading binary file. " << ret << " bytes read instead of 8." << c_endl;
    return 2;
  }

  return 0;  
}

//----------------------------------------------------------------------------------------------------

signed int SlinkFile::CheckSequence(unsigned long long n, const std::vector<sequence> &seqs)
{
  for (unsigned int i = 0; i < seqs.size(); i++)
    if (seqs[i].compare(n)) return i;
    
  return -1;
}

//----------------------------------------------------------------------------------------------------

char SlinkFile::SeekNextFrame()
{
  // check if file is open
  if (!file)
    return 1;

  bool found = false;
  unsigned long long n = 0;  

#ifdef USE_CASTOR
  while (!rfio_feof(file)) {
#else
  while (!feof(file)) {
#endif
    Read(n);
    if (CheckSequence(n, seqBOF) >= 0) {
      found = true;
      break;
    }
  }

  // remember position if indexing
#ifdef USE_CASTOR
  if (indexStatus == isIndexing && found) positions.push_back(rfio_ftell(file));
#else
  if (indexStatus == isIndexing && found) positions.push_back(ftell(file));
#endif

#ifdef USE_CASTOR
      int eofFlag = rfio_feof(file);
#else
      int eofFlag = feof(file);
#endif

  // if eof is found and was indexing, set indexed
  if (indexStatus == isIndexing && eofFlag)
    indexStatus = isIndexed;

  // return code
  if (!found)
    return 2;
  else
    return 0;
}

//----------------------------------------------------------------------------------------------------

unsigned char SlinkFile::LoadFrameCheckFEOF(const char *msg)
{
  // check for EOF
#ifdef USE_CASTOR
  int eofFlag = rfio_feof(file);
#else
  int eofFlag = feof(file);
#endif
  if (eofFlag) {
    // print message
    ERROR("SlinkFile::LoadFrame") << "Corrupted Slink frame (event " << positions.size() << "). " << msg << c_endl;

    // if indexing, set index and remove the last position - it's not complete
    if (indexStatus == isIndexing) {
      indexStatus = isIndexed;
      positions.pop_back();
    }

    return 1;
  }

  return 0;
}

//----------------------------------------------------------------------------------------------------

char SlinkFile::LoadFrame()
{
  if (LoadFrameCheckFEOF("File ends after BOF sequence. Frame ignored.")) return 10;

  unsigned long long n;

  // header
  // the line with flags (f7f7) whether fiber is connected
  unsigned long long mask = 0;
  Read(n);
  if ((n & 0x000000000000FFFFULL) == 0x000000000000F7F7ULL) mask |= 0x000000000000FFFFULL;
  if ((n & 0x00000000FFFF0000ULL) == 0x00000000F7F70000ULL) mask |= 0x00000000FFFF0000ULL;
  if ((n & 0x0000FFFF00000000ULL) == 0x0000F7F700000000ULL) mask |= 0x0000FFFF00000000ULL;
  if ((n & 0xFFFF000000000000ULL) == 0xF7F7000000000000ULL) mask |= 0xFFFF000000000000ULL;
  if (LoadFrameCheckFEOF("File ends after 1st line of header. Frame ignored.")) return 11;
  if (!mask)
    WARN("SlinkFile::LoadFrame") << "Frame mask is zero. No data will be read from this Slink frame. (row = 0x" << hex << n << dec << ")." << c_endl;

  // the line with tag
  Read(n);
  frame.SetTag(n);
  if (LoadFrameCheckFEOF("File ends after 2nd line of header. Frame ignored.")) return 12;
  if (n) WARN("SlinkFile::LoadFrame") << "Frame tag is " << hex << n << dec << " instead of 0 as it should be for this version." << c_endl;

  // data block
  unsigned long long* addr = (unsigned long long*) frame.GetSlinkFrame();
  for (int i = 0; i < VFAT_FRAME_BIT_SIZE; i++) {  
    // read data and remember previous position
    if (LoadFrameCheckFEOF("File ends while reading data block. Frame ignored.")) return 13;
#ifdef USE_CASTOR
    long int prevPos = rfio_ftell(file);
#else
    long int prevPos = ftell(file);
#endif
    Read(n);

    if (CheckSequence(n, seqBOF) >= 0) {
      ERROR("SlinkFile::LoadFrame") << "Corrupted Slink frame (event " << positions.size() << "). A BOF sequence found in data block (at position " 
        << i << "). Frame dropped." << c_endl;
      Myfseek(file, prevPos, SEEK_SET); // go one record back so as SeekNextFrame can find it
      if (indexStatus == isIndexing)
        positions.pop_back();
      return 1;
    }

    if (CheckSequence(n, seqEOF) >= 0) {
      ERROR("SlinkFile::LoadFrame") << "Corrupted Slink frame (event " << positions.size() << "). An EOF sequence found in data block (at position "
        << i << "). Frame dropped.)" << c_endl;
      if (indexStatus == isIndexing)
        positions.pop_back();
      return 1;
    }

    // save value
    addr[i] = n & mask;
  }

  // footer
#ifdef USE_CASTOR
  long int prevPos = rfio_ftell(file);
#else
  long int prevPos = ftell(file);
#endif

  Read(n);
  if (CheckSequence(n, seqBOF) >= 0) {
    ERROR("SlinkFile::LoadFrame") << "Corrupted Slink frame (event " << positions.size()
      << "). A BOF sequence found in place of EOF sequence. Frame dropped.)" << c_endl;
    Myfseek(file, prevPos, SEEK_SET); // go one record back so as SeekNextFrame can find it
    if (indexStatus == isIndexing)
      positions.pop_back();
    return 1;
  }

  if (CheckSequence(n, seqEOF) < 0) {
    ERROR("SlinkFile::LoadFrame") << "Corrupted Slink frame (event " << positions.size() << "). EOF sequence not found after data block. Frame dropped.)" << c_endl;
    if (indexStatus == isIndexing)
      positions.pop_back();
    return 1;
  }

  // everything OK
  return 0;  
}

//----------------------------------------------------------------------------------------------------

unsigned char SlinkFile::GetNextEvent(RawEvent *event)
{
  while (true) {
    // try to move to the next frame beggining
     if (SeekNextFrame()) return 1;    // eof without finding BOF or file not open

    // read data block to frame
     int resLoadFrame = LoadFrame();
    if (resLoadFrame == 1) corrNum++;  // corrupted frame
      if (resLoadFrame == 0) {       // success
      LoadEvent(event);
      return 0;  
    }
      if (resLoadFrame > 1)  return 2;     // eof
    
    }
}

//----------------------------------------------------------------------------------------------------

unsigned char SlinkFile::GetEvent(unsigned long n,RawEvent *event)
{
  // check if indexed and n <= number of frames
  if (indexStatus != 2) return 1;
  if (n >= positions.size()) return 2;

  // load the frame
  Myfseek(file, positions[n], SEEK_SET);
  int resLoadFrame = LoadFrame();
  LoadEvent(event);
  return resLoadFrame;
}

// myfseek function, eliminating bug with shift library
int SlinkFile::Myfseek(FILE* stream, long int offset, int origin) {
#ifdef USE_CASTOR
  if (rfio_feof(stream)) 
    Reopen();  
  else rfio_fseek(stream, offset, origin);
#else
  return fseek(stream, offset, origin);
#endif
  return 1;
}

unsigned int SlinkFile::Reopen()
{
#ifdef USE_CASTOR
  rfio_fclose(file);
  file = rfio_fopen(filename, const_cast<char*>("rb"));
#else
  fclose(file);
  file = fopen(filename, "rb");
#endif
  if (!file) {
      ERROR("SlinkFile::Reopen") << "Could not reopen file `" << filename << "'." << c_endl;
    return 1;
  }

  return 0;
}

} // namespace

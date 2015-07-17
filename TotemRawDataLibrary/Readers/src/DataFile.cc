/****************************************************************************
*
* This is a part of the TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/


#include "TotemRawDataLibrary/DataFormats/interface/RawEvent.h"
#include "TotemRawDataLibrary/DataFormats/interface/CommonDef.h"
#include "TotemRawDataLibrary/Readers/interface/DataFile.h"
#include "TotemRawDataLibrary/Readers/interface/SlinkFile.h"
#include "TotemRawDataLibrary/Readers/interface/MultiSlinkFile.h"
#include "TotemRawDataLibrary/Readers/interface/VMEFile.h"
#include "TotemRawDataLibrary/Readers/interface/VME2File.h"
#include "TotemRawDataLibrary/Readers/interface/VMEAFile.h"
#include "TotemRawDataLibrary/Readers/interface/VMEAStream.h"
#include "TotemRawDataLibrary/Readers/interface/VMEBFile.h"
#include "TotemRawDataLibrary/Readers/interface/TTPFile.h"



using namespace std;

//#define DEBUG 0

namespace Totem {

#if 1
  // fake implementation (vtable problem ...)
  DataFile::OpenStatus DataFile::Open(const string&) { return osOK; }
  void DataFile::Close() {}
  unsigned char DataFile::GetNextEvent(RawEvent *) { return 0; }
  unsigned char DataFile::GetEvent(unsigned long, RawEvent *) { return 0; }
  void DataFile::Rewind() {}
  void DataFile::StartIndexing() {}
  unsigned long DataFile::GetEventsCount() const { return 0; }
  unsigned long DataFile::GetEventsCorrupted() const { return 0; }
  bool DataFile::RandomAccessSupported() const { return true; }
#endif

//----------------------------------------------------------------------------------------------------


DataFile* DataFile::OpenStandard(const string &fileName)
{
#ifdef DEBUG
  printf(">> size test\n");
  printf("\tsize of int : %lu\n", sizeof(int));
  printf("\tsize of long : %lu\n", sizeof(long));
  printf("\tsize of long long : %lu\n", sizeof(long));
  printf("\tsize of std::streampos : %lu\n", sizeof(std::streampos));
  printf("\tsize of CircularBuffer::position_type : %lu\n", sizeof(CircularBuffer<unsigned int>::position_type));

  CircularBuffer<unsigned int>::position_type pt = 0xFF;
  streamoff one(1);
  for (unsigned int i = 8; i < 65; i += 8)
  {
    cout << dec << i << ": " << hex << pt << ", ";
    pt = pt + one;
    cout << pt << endl;
    pt = (pt - one) * 256;
    one = one * 256;
  }

  printf(">> Trying to open file %s.\n", fileName.c_str());
#endif

  // try to open the file with all known file formats
  for (unsigned char i = 0; ; i++)
  {
    DataFile *f = NULL;
    switch (i)
    {
      case 0: f = new VMEAStream(); break;
      case 1: f = new VMEBFile(); break;
      case 2: f = new VMEAFile(); break;
      case 3: f = new VME2File(); break;
      case 4: f = new VMEFile(); break;
      case 5: f = new MultiSlinkFile(); break;
      case 6: f = new TTPFile(); break;
      case 7: f = new SlinkFile(); break;
    }

    // this shall never happen as SlinkFile accepts everything
    if (!f)
    {
      printf("ERROR in DataFile::OpenStandard: error creating DataFile instance (i = %u).\n", i);
      break;
    }

    OpenStatus r = f->Open(fileName);

#ifdef DEBUG
    printf("\tclass `%s' gives result %u\n", f->GetClassName().c_str(), r);
#endif

    if (r == osOK)
    {
      printf("DataFile::OpenStandard: File `%s' has been opened by class `%s'.\n", fileName.c_str(), f->GetClassName().c_str());
      return f;
    }

    if (r == osCannotOpen)
    {
      printf("Error in DataFile::OpenStandard: File `%s' cannot be opened.\n", fileName.c_str());
      delete f;
      return NULL;
    }

    if (r == osWrongFormat)
    {
      // do not print anything: this situation is normal as one tries to open the file with all known file readers
      //printf("DataFile::OpenStandard: osWrongFormat\n");

      delete f;
    }
  }

  return NULL;
}

} // namespace


/****************************************************************************
*
* This is a part of the TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "TotemRawDataLibrary/DataFormats/interface/CommonDef.h"
#include "TotemRawDataLibrary/Readers/interface/DataFile.h"
#include "TotemRawDataLibrary/DataFormats/interface/VFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/VFATFrameCollection.h"
#include "TotemRawDataLibrary/Utilities/interface/CablingAnalyzer.h"

#ifdef _MONITOR_
	#ifdef _MONITOR_QT3_
		#include "VFAT2RegistersXML.h"
		#define REG_CLASS VFAT2RegistersXML
	#else
		#include "RegistersInXmlAttributes.h"
		#define REG_CLASS RegistersInXmlAttributes
	#endif
#else
  // TODO: solve the problem with missing VFAT2RegistersXML
#endif

#include <vector>


using namespace Totem;


int main(int argc, char *argv[])
{
#ifdef _MONITOR_
  setvbuf(stdout, NULL, _IONBF, 0);
  setvbuf(stderr, NULL, _IONBF, 0);

  if (argc != 3) {
    printf("Usage: cabling <data file> <I2C file>\n");
    return 1;
  }

  // try to open files
  DataFile* input = DataFile::OpenStandard(argv[1]);
  if (!input) {
    printf("Error in opening data file.\n");
    return 2;
  }

  REG_CLASS conf;
  if (conf.OpenFile(argv[2]))
  {
    printf("Error in opening I2C file.\n");
    return 3;
  }

  // run test
  CablingAnalyzer test;
  unsigned short result = test.run(input, &conf);
  if (result) return result;

  delete input;
#else
  printf("SORRY, for the time being, this program can only work in the Monitor framework.\n");
#endif

  return 0;
}

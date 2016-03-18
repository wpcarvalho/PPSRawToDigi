/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#include "EventFilter/TotemRawToDigi/interface/RawDataUnpacker.h"

//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

RawDataUnpacker::RawDataUnpacker(const edm::ParameterSet &conf)
{
}

//----------------------------------------------------------------------------------------------------

int RawDataUnpacker::Run(int fedId, const FEDRawData &data, VFATFrameCollection &coll, TotemRawEvent &rawEvent)
{
  /*
  unsigned short int length = fedData.size();
  const ScalersEventRecordRaw_v6 *raw = (struct ScalersEventRecordRaw_v6 *)fedData.data();
  */

  // TODO
  return 0;
}

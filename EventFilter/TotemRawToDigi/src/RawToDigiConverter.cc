/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#include "EventFilter/TotemRawToDigi/interface/RawToDigiConverter.h"

//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

RawToDigiConverter::RawToDigiConverter(const edm::ParameterSet &conf)
{
}

//----------------------------------------------------------------------------------------------------

int RawToDigiConverter::Run(const VFATFrameCollection &coll,
  const TotemDAQMapping &mapping, const TotemAnalysisMask &mask,
  edm::DetSetVector<RPStripDigi> &rpData, std::vector <RPCCBits> &rpCC, TotemRawToDigiStatus &status)
{
  // TODO
  return 0;
}

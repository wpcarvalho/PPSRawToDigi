/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#include "EventFilter/TotemRawToDigi/RawDataUnpacker.h"

//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

RawDataUnpacker::RawDataUnpacker(const edm::ParameterSet &conf)
{
}

//----------------------------------------------------------------------------------------------------

int RawDataUnpacker::Run(int fedId, const edm::FEDRawData &data, VFATFrameCollection &coll)
{
}

/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#ifndef _RawDataUnpacker_h_
#define _RawDataUnpacker_h_

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/FEDRawData/interface/FEDRawData.h"

#include "DataFormats/TotemRawData/interface/TotemRawEvent.h"

#include "EventFilter/TotemRawToDigi/interface/VFATFrameCollection.h"

//----------------------------------------------------------------------------------------------------

/// \brief Collection of code for unpacking of TOTEM raw-data.
class RawDataUnpacker
{
  public:
    RawDataUnpacker(const edm::ParameterSet &conf);

    /// unpack data from FED with fedId into `coll' collection
    int Run(int fedId, const FEDRawData &data, VFATFrameCollection &coll, TotemRawEvent &rawEvent);
};

#endif

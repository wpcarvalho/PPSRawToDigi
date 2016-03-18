/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#ifndef _RawToDigiConverter_h_
#define _RawToDigiConverter_h_

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "EventFilter/TotemRawToDigi/interface/VFATFrameCollection.h"

//----------------------------------------------------------------------------------------------------

/// \brief Collection of code to convert TOTEM raw data into digi.
class RawToDigiConverter
{
  public:
    RawToDigiConverter(const edm::ParameterSet &conf);

    /// Converts vfat data in `coll'' into digi.
    // TODO: add outputs
    int Run(const VFATFrameCollection &coll);
};

#endif

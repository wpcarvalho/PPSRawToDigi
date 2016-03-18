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
#include "DataFormats/Common/interface/DetSetVector.h"

#include "EventFilter/TotemRawToDigi/interface/VFATFrameCollection.h"

#include "CondFormats/TotemReadoutObjects/interface/TotemDAQMapping.h"
#include "CondFormats/TotemReadoutObjects/interface/TotemAnalysisMask.h"

#include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"
#include "DataFormats/TotemL1Trigger/interface/RPCCBits.h"
#include "DataFormats/TotemRawData/interface/TotemRawEvent.h"
#include "DataFormats/TotemRawData/interface/TotemRawToDigiStatus.h"

//----------------------------------------------------------------------------------------------------

/// \brief Collection of code to convert TOTEM raw data into digi.
class RawToDigiConverter
{
  public:
    RawToDigiConverter(const edm::ParameterSet &conf);

    /// Converts vfat data in `coll'' into digi.
    // TODO: add outputs
    int Run(const VFATFrameCollection &coll,
      const TotemDAQMapping &mapping, const TotemAnalysisMask &mask,
      edm::DetSetVector<RPStripDigi> &rpData, std::vector<RPCCBits> &rpCC, TotemRawToDigiStatus &status);
};

#endif

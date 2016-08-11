/****************************************************************************
 *  Seyed Mohsen Etesami
 ****************************************************************************/

#ifndef EventFilter_CTPPSRawToDigi_DiamondRawDataUnpacker
#define EventFilter_CTPPSRawToDigi_DiamondRawDataUnpacker

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/CTPPSDigi/interface/DiamondFEDInfo.h"

#include "EventFilter/CTPPSRawToDigi/interface/DiamondVFATInterface.h"
#include "EventFilter/CTPPSRawToDigi/interface/DiamondVFATFrameCollection.h"

//----------------------------------------------------------------------------------------------------

/// \brief Collection of code for unpacking of Diamond timing detector raw-data.
class DiamondRawDataUnpacker
{
 public:
  typedef uint64_t word;

  /// VFAT transmission modes
  enum { vmCluster = 0x80, vmRaw = 0x90 };

  DiamondRawDataUnpacker() {}
    
  DiamondRawDataUnpacker(const edm::ParameterSet &conf);

  /// Unpack data from FED with fedId into `coll' collection.
  int Run(int fedId, const FEDRawData &data, std::vector<DiamondFEDInfo> &fedInfoColl, DiamondVFATFrameCollection &coll) const;

  /// Process one Opto-Rx (or LoneG) frame.
  int ProcessOptoRxFrame(const word *buf, unsigned int frameSize, DiamondFEDInfo &fedInfo, DiamondVFATFrameCollection *fc) const;


  /// Process one Opto-Rx frame in parallel (new) format
  int ProcessOptoRxFrameParallel(const word *buffer, unsigned int frameSize, DiamondFEDInfo &fedInfo, DiamondVFATFrameCollection *fc) const;

  /// Process data from one VFAT in parallel (new) format
  int ProcessVFATDataParallel(const uint16_t *buf, unsigned int OptoRxId, DiamondVFATFrameCollection *fc) const;


  /// Process data from one VFAT in parallel Diamond Format
  int ProcessVFATDataFED(const uint16_t *buf, unsigned int FEDId, DiamondVFATFrameCollection *fc) const;
};

#endif

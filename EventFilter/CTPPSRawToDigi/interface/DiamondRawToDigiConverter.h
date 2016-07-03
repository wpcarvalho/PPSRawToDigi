/****************************************************************************
*  Seyed Mohsen Etesami
****************************************************************************/

#ifndef EventFilter_CTPPSRawToDigi_DimondRawToDigiConverter
#define EventFilter_CTPPSRawToDigi_DiamondRawToDigiConverter

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "EventFilter/CTPPSRawToDigi/interface/DiamondVFATInterface.h"

#include "CondFormats/CTPPSReadoutObjects/interface/TotemDAQMappingDiamond.h"
#include "CondFormats/CTPPSReadoutObjects/interface/DiamondAnalysisMask.h"

#include "DataFormats/CTPPSDigi/interface/DiamondDigi.h"
#include "DataFormats/CTPPSDigi/interface/DiamondVFATStatus.h"

//----------------------------------------------------------------------------------------------------

/// \brief Collection of code to convert PPS Timing raw data into digi.
class DiamondRawToDigiConverter
{
  public:
    DiamondRawToDigiConverter(const edm::ParameterSet &conf);

    /// Creates Diamond digi.
    void Run(const DiamondVFATInterface &coll, const TotemDAQMappingDiamond &mapping, const DiamondAnalysisMask &mask,
      edm::DetSetVector<DiamondDigi> &digi, edm::DetSetVector<DiamondVFATStatus> &status);

    /// Print error summaries.
    void PrintSummaries();

  private:
    struct Record
    {
      const DiamondVFATInfo *info;
      const DiamondVFATFrame *frame;
      DiamondVFATStatus status;
    };

    unsigned char verbosity;
    
    unsigned int printErrorSummary;
    unsigned int printUnknownFrameSummary;

    enum TestFlag { tfNoTest, tfWarn, tfErr };

    /// flags for which tests to run
    unsigned int testFootprint;
    unsigned int testCRC;
    unsigned int testID;
    unsigned int testECRaw;
    unsigned int testECDAQ;
    unsigned int testECMostFrequent;
    unsigned int testBCMostFrequent;

    /// the minimal required number of frames to determine the most frequent counter value
    unsigned int EC_min, BC_min;

    /// the minimal required (relative) occupancy of the most frequent counter value to be accepted
    double EC_fraction, BC_fraction;

    /// error summaries
    std::map<DiamondFramePosition, std::map<DiamondVFATStatus, unsigned int> > errorSummary;
    std::map<DiamondFramePosition, unsigned int> unknownSummary;

    /// Common processing for all VFAT based sub-systems.
    void RunCommon(const DiamondVFATInterface &input, const TotemDAQMappingDiamond &mapping,
      std::map<DiamondFramePosition, Record> &records);
};

#endif

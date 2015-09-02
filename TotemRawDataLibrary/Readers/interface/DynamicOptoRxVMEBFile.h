/****************************************************************************
*
* This is a part of the TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors: 
*
****************************************************************************/

#ifndef _Totem_DynamicOptoRxVMEBFile_h_
#define _Totem_DynamicOptoRxVMEBFile_h_

#include "TotemRawDataLibrary/Readers/interface/BaseVMEBFile.h"
#include "TotemRawDataLibrary/Readers/interface/CircularBuffer.h"
#include "TotemRawDataLibrary/DataFormats/interface/OptoRxSupplementalData.h"
#include "TotemRawDataLibrary/DataFormats/interface/PositionedVFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/RawVFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/ClusterizationVFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/PositionedVFATFrameCollection.h"
#include "TotemRawDataLibrary/DataFormats/interface/SecondLevelTriggerVFATFrame.h"

#include <vector>
#include <map>

namespace Totem {

/**
 * \ingroup TotemRawDataLibrary
 * TODO: describe
 **/ 
class DynamicOptoRxVMEBFile : public BaseVMEBFile
{
  public:
    DynamicOptoRxVMEBFile();
    virtual ~DynamicOptoRxVMEBFile();

    virtual VFATFrameCollection* CreateCollection() const{
        return new PositionedVFATFrameCollection;
    }
    virtual bool IsCollectionCompatible(VFATFrameCollection *frames) const{
        return (frames->GetClassName().compare("PositionedVFATFrameCollection") == 0);
    }
    virtual std::string GetClassName() const {
        return "DynamicOptoRxVMEBFile";
    }

    ///\brief processes a VMEB Event
    /// returns the number of GOH blocks that failed consistency checks
    static unsigned int ProcessVMEBEvent(char *ptr,  RawEvent *);
    virtual unsigned char GetNextEvent(RawEvent*);

  protected:
    /// process one LDC event
    /// returns the number of GOH blocks that failed consistency checks
    static unsigned int ProcessSubEvent(char *ptr, RawEvent *);
    virtual std::string getExtension() const {
        return ".UNKNOWN"; //todo backpatch it once extension is known
    }

    static unsigned int ProcessOptoRxFrame(unsigned long long* buffer, unsigned int frameSize, RawEvent *event);
    static bool IsCrcForOptoRxFrameValid(unsigned long long* buffer, unsigned int frameSize);
    static unsigned int ProcessRawVFATFrames(std::vector<PositionedVFATFrame*> framesRead, RawEvent *event, OptoRxSupplementalData* optoRx);
    static unsigned int ProcessClusterizationVFATFrames(std::vector<PositionedVFATFrame*> framesRead, RawEvent *event, OptoRxSupplementalData* optoRx);
    static unsigned int ProcessMixedRawAndClusterizationVFATFrames(std::vector<PositionedVFATFrame*> framesRead, RawEvent *event, OptoRxSupplementalData* optoRx);

  private:
    static int32_t crcLookupTable[256]; /* we need a sign for initialization; we'll use only positive values  */
    static void initializeCrcLookupTable();

    static unsigned int FilterOutGohWithoutAllFibers(std::vector<PositionedVFATFrame*>* framesRead, OptoRxSupplementalData* optoRx);
    static unsigned int FilterOutInconsitientGoh(std::vector<PositionedVFATFrame*>* framesRead, OptoRxSupplementalData* optoRx);
    static unsigned int FilterOutInvalidFrames(std::vector<PositionedVFATFrame*>* framesRead, OptoRxSupplementalData* optoRx);
    static void AddFrameToCollection(PositionedVFATFrame* frame, OptoRxSupplementalData* optoRx, PositionedVFATFrameCollection* frames);
};

}
#endif

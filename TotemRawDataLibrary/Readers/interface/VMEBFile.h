/****************************************************************************
*
* This is a part of the TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors:
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#ifndef _Totem_VMEBFile_h_
#define _Totem_VMEBFile_h_

#include "TotemRawDataLibrary/Readers/interface/BaseVMEBFile.h"
#include "TotemRawDataLibrary/Readers/interface/CircularBuffer.h"
#include "TotemRawDataLibrary/DataFormats/interface/OptoRxVFATFrameCollection.h"

#include <vector>
//#include <fstream>

namespace Totem {

/**
 * \ingroup TotemRawDataLibrary
 * Reads a file in VMEB format, without using the DAQA monitoring library.
**/

class VMEBFile : public BaseVMEBFile
{
  public:
    VMEBFile();
    virtual ~VMEBFile();

    virtual VFATFrameCollection* CreateCollection() const
	{
        return new OptoRxVFATFrameCollection;
    }

    virtual bool IsCollectionCompatible(VFATFrameCollection *c) const
	{
        return (c->GetClassName().compare("OptoRxVFATFrameCollection") == 0);
    }

    virtual std::string GetClassName() const
	{
        return "VMEBFile";
    }

    /// Processes a VMEB Event.
    /// returns the number of GOH blocks that failed consistency checks
    static unsigned int ProcessVMEBEvent(char *ptr, OptoRxVFATFrameCollection *, RawEvent *);

    virtual unsigned char GetNextEvent(RawEvent*);

  protected:
    /// Process one LDC event.
    /// returns the number of GOH blocks that failed consistency checks
    static unsigned int ProcessSubEvent(char *ptr, OptoRxVFATFrameCollection *, RawEvent *);

    virtual std::string getExtension() const
	{
        return ".vmeb";
    }
};

}
#endif

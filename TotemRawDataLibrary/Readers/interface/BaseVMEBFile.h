#ifndef _Totem_BaseVMEBFile_h_
#define _Totem_BaseVMEBFile_h_

#include "TotemRawDataLibrary/Readers/interface/DataFile.h"
#include "TotemRawDataLibrary/Readers/interface/CircularBuffer.h"
#include "TotemRawDataLibrary/DataFormats/interface/OptoRxVFATFrameCollection.h"

#include <vector>

namespace Totem {

/**
 * \ingroup TotemRawDataLibrary
 * TODO: describe
 **/ 
class BaseVMEBFile : public DataFile
{
public:
    /// standard equipment types
    enum { etOptoRxOld = 100, etLoneG = 110, etOptoRx = 120, etDynamicOptoRx = 999 /*todo UNKNOWN */ };
    typedef unsigned long long word;

    BaseVMEBFile();
    virtual ~BaseVMEBFile();

    virtual OpenStatus Open(const std::string &);
    virtual void Close();
    virtual void Rewind();
    virtual unsigned char GetEvent(unsigned long, RawEvent*);
    virtual unsigned char GetNextEvent(RawEvent*);

    virtual VFATFrameCollection* CreateCollection() const = 0;
    virtual bool IsCollectionCompatible(VFATFrameCollection *c) const = 0;
    virtual std::string GetClassName() const = 0;

    virtual unsigned long GetEventsCount() const {
        return positions.size();
    }
    virtual unsigned long GetEventsCorrupted() const {
        return corruptedEventCounter;
    }
    virtual void StartIndexing(){
        positions.clear();
        indexStatus = isIndexing;
    }
    virtual void StopIndexing(){
        indexStatus = isIndexed;
    }
    virtual bool RandomAccessSupported() const{
        return true;
    }

protected:
    static const unsigned int eventHeaderSize;

    virtual std::string getExtension() const = 0;
    /// reads 'bytesToRead' bytes from the file to buffer, starting at the given offset
    virtual unsigned char ReadToBuffer(unsigned int bytesToRead, unsigned int offset);
    /// processes a LoneG frame (trigger data)
    static unsigned int ProcessLoneGFrame(word *ptr, signed long size, RawEvent *);

    /// data pointer, to be allocated one time only
    char *dataPtr;
    /// data buffer size
    unsigned int dataPtrSize;
    /// stores event positions in file
    std::vector<std::streampos> positions;
    /// input file pointer
    FILE* infile;
    /// number of corrupted frames
    unsigned int corruptedEventCounter;
    /// name of opened file
    std::string filename;

#ifdef USE_CASTOR
    /// used to reset eof flag when using rfio_api
    /// returns 0 if OK
    virtual unsigned int Reopen();
#endif
};

}
#endif

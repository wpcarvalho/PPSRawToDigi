#ifndef _Totem_SRSFile_h_
#define _Totem_SRSFile_h_

#include "TotemRawDataLibrary/Readers/interface/DataFile.h"
#include "TotemRawDataLibrary/DataFormats/interface/OptoRxVFATFrameCollection.h"
#include "TotemRawDataLibrary/Readers/interface/StorageFile.h"

#include <vector>

namespace Totem {

/**
 * Reads a file in SRS format.
 **/ 
class SRSFile : public DataFile
{
public:
    /// standard equipment types
    enum { etOptoRx = 22 };

    /// 1 word = 64 bits = 8 bytes
    typedef unsigned long long word;

    SRSFile();
    virtual ~SRSFile();

    virtual OpenStatus Open(const std::string &);
    virtual OpenStatus Open(StorageFile* storageFile);

    virtual void Close();
    virtual void Rewind();
    virtual unsigned char GetEvent(unsigned long, RawEvent*);
    virtual unsigned char GetNextEvent(RawEvent*);

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
        return "SRSFile";
    }

    virtual unsigned long GetEventsCount() const
    {
        return positions.size();
    }

    virtual unsigned long GetEventsCorrupted() const
    {
        return corruptedEventCounter;
    }

    virtual void StartIndexing()
    {
        positions.clear();
        indexStatus = isIndexing;
    }

    virtual void StopIndexing()
    {
        indexStatus = isIndexed;
    }

    virtual bool RandomAccessSupported() const
    {
        return true;
    }

    /// Processes an SRS Event.
    /// returns the number of GOH blocks that failed consistency checks
    static unsigned int ProcessSRSEvent(char *ptr, OptoRxVFATFrameCollection *, RawEvent *);

    /// TODO: update
    /// Extract VFAT data from an OptoRx frame and saves them to the collection.
    /// BufferType can be any type that has [] operator (array, CircularBuffer, ...)
    /// returns the number of cells that failed consistency checks
    /// \param buf is the OptoRx frame buffer, including the header and footer
    /// \param frameSize is the size (in words) of the OptoRx frame including the header and footer
    static unsigned int ProcessOptoRxFrame(word *buf, unsigned int frameSize,
                                           OptoRxVFATFrameCollection *, RawEvent *);

protected:
    static const unsigned int eventHeaderSize;

    /// Process one LDC event.
    /// returns the number of GOH blocks that failed consistency checks
    static unsigned int ProcessSubEvent(char *ptr, OptoRxVFATFrameCollection *, RawEvent *);

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
 	StorageFile *infile;

	/// number of corrupted frames
    unsigned int corruptedEventCounter;
   
   	/// name of opened file
    std::string filename;
};


}
#endif

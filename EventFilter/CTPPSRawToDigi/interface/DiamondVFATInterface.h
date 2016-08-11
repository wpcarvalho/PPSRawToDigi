/****************************************************************************
 *  Seyed Mohsen Etesami
 ****************************************************************************/


#ifndef EventFilter_CTPPSRawToDigi_DiamondVFATInterface
#define EventFilter_CTPPSRawToDigi_DiamondVFATInterface

#include "CondFormats/CTPPSReadoutObjects/interface/DiamondFramePosition.h"

#include "EventFilter/CTPPSRawToDigi/interface/DiamondVFATFrame.h"

#include <string>

/**
 * Interface class for VFAT-frame collections.
 **/
class DiamondVFATInterface 
{
 public:
  DiamondVFATInterface() {}
  virtual ~DiamondVFATInterface() {}

  /// returns pointer to frame with ID, performs NO duplicity check (if there is precisely one frame with this 12bit ID)
  virtual const DiamondVFATFrame* GetFrameByID(unsigned int ID) const = 0;

  /// returns frame at given position in Slink frame
  virtual const DiamondVFATFrame* GetFrameByIndex(DiamondFramePosition index) const = 0;

  /// returns frame at given position in Slink frame and checks 12bit ID
  virtual const DiamondVFATFrame* GetFrameByIndexID(DiamondFramePosition index, unsigned int ID);

  /// return the number of VFAT frames in the collection
  virtual unsigned int Size() const = 0;

  /// returns whether the collection is empty
  virtual bool Empty() const = 0;

  /// pair: frame DAQ position, frame data
  typedef std::pair<DiamondFramePosition, const DiamondVFATFrame*> value_type;

  /// the DiamondVFATInterface interator
  class Iterator {
  protected:
    /// interator value
    value_type value;

    /// the pointer to the collection
    const DiamondVFATInterface* collection;

  public:
    /// constructor, automatically sets the iterator to the beginning
  Iterator(const DiamondVFATInterface* c = NULL) : collection(c)
    { if (collection) value = collection->BeginIterator(); }

    /// returns the DAQ position of the current element
    DiamondFramePosition Position()
    { return value.first; }

    /// returns the frame data of the current element
    const DiamondVFATFrame* Data()
    { return value.second; }

    /// shifts the iterator
    void Next()
    { value = collection->NextIterator(value); }

    /// returns whether the iterator points over the end of the collection
    bool IsEnd()
    { return collection->IsEndIterator(value); }
  };
    
 protected:
  /// returns the beginning of the collection
  virtual value_type BeginIterator() const = 0;

  /// shifts the iterator
  virtual value_type NextIterator(const value_type&) const = 0;

  /// checks whether the iterator points over the end of the collection
  virtual bool IsEndIterator(const value_type&) const = 0;
};

#endif

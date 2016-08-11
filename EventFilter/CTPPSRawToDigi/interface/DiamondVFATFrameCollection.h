/****************************************************************************
 *   Seyed Mohsen Etesami    
 ****************************************************************************/


#ifndef EventFilter_CTPPSRawToDigi_DiamondVFATFrameCollection
#define EventFilter_CTPPSRawToDigi_DiamondVFATFrameCollection

#include "EventFilter/CTPPSRawToDigi/interface/DiamondVFATInterface.h"

#include <map>

/**
 * A basic implementation of VFAT frame collection, as map: DiamondFramePosition --> VFATFrame.
 **/
class DiamondVFATFrameCollection : public DiamondVFATInterface
{
 protected:
  typedef std::map<DiamondFramePosition, DiamondVFATFrame> MapType;

  MapType data;

  virtual value_type BeginIterator() const;
  virtual value_type NextIterator(const value_type&) const;
  virtual bool IsEndIterator(const value_type&) const;

 public:
  DiamondVFATFrameCollection();
  ~DiamondVFATFrameCollection();

  const DiamondVFATFrame* GetFrameByID(unsigned int ID) const;
  const DiamondVFATFrame* GetFrameByIndex(DiamondFramePosition index) const;

  virtual unsigned int Size() const
  {
    return data.size();
  }

  virtual bool Empty() const
  {
    return (data.size() == 0);
  }

  void Insert(const DiamondFramePosition &index, const DiamondVFATFrame &frame)
  {
    data.insert({index, frame});
  }

  /// inserts an empty (default) frame to the given position and returns pointer to the frame
  DiamondVFATFrame* InsertEmptyFrame(DiamondFramePosition index)
  {
    return &data.insert({index, DiamondVFATFrame()}).first->second;
  }

  /// cleans completely the collection
  void Clear()
  {
    data.clear();
  }
};

#endif

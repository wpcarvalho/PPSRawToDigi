/****************************************************************************
*
* This is a part of the TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors: 
*   Mate Csanad (mate.csanad@cern.ch)
*   Jan Kašpar (jan.kaspar@gmail.com) 
*    
****************************************************************************/


#ifndef _Totem_SimpleVFATFrameCollection_h_
#define _Totem_SimpleVFATFrameCollection_h_

#include "TotemRawDataLibrary/DataFormats/interface/VFATFrameCollection.h"
#include "TotemRawDataLibrary/DataFormats/interface/OldVFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/VFATFrame.h"

#include <map>

  #include <cstdio>

namespace Totem {

/**
\ingroup TotemRawDataLibrary
\brief A simple implementation of VFAT frame collection.

Based on std::map<FramePosition, VFATFrame>
**/
class SimpleVFATFrameCollection : public VFATFrameCollection
{
  protected:
    typedef std::map<FramePosition, OldVFATFrame> MapType;

    MapType data;

    virtual value_type BeginIterator() const;
    virtual value_type NextIterator(const value_type&) const;
    virtual bool IsEndIterator(const value_type&) const;

  public:
    SimpleVFATFrameCollection();
    ~SimpleVFATFrameCollection();

    virtual std::string GetClassName() const
      { return "SimpleVFATFrameCollection"; }

    const OldVFATFrame* GetFrameByID(unsigned int ID) const;
    const OldVFATFrame* GetFrameByIndex(FramePosition index) const;

    virtual unsigned int Size() const
      { return data.size(); }

    virtual bool Empty() const
      { return (data.size() == 0); }

    /// inserts an empty (default) frame to the given position and returns pointer to the frame
    OldVFATFrame* InsertEmptyFrame(FramePosition index)
      { return &data.insert(std::pair<FramePosition, OldVFATFrame>(index, OldVFATFrame())).first->second; }

    /// cleans completely the collection
    void Clear()
      { data.clear(); }
};

} // namespace
#endif


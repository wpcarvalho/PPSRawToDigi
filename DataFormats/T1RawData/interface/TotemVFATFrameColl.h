/**************************************************************
 * $Id: TotemVFATFrameColl.h,v 1.1 2008/09/12 16:17:05 lgrzanka Exp $
 * $Revision: 1.1 $
 * $Date: 2008/09/12 16:17:05 $
 ***************************************************************/

#ifndef _TotemVFATFrameColl_h_
#define _TotemVFATFrameColl_h_

#include "DataFormats/T1RawData/interface/TotemRawVFATFrame.h"
#include <vector>
#include <ext/hash_map>

#define NUMBER_OF_VFATS 64

class TotemVFATFrameColl
{
 public:
  TotemVFATFrameColl(char PlugVFATobjectsManualy = 0);	/// constructor, by default puts 'NUMBER_OF_VFATS' VFATframes into frames array
  ~TotemVFATFrameColl();									/// destructor

  std::vector<TotemRawVFATFrame>& GetFrames() 			/// getter for frame array
    { return frames; };
  const	std::vector<TotemRawVFATFrame>& GetFrames() const 			/// getter for frame array
    { return frames; };
  TotemRawVFATFrame* operator[] (int i)					/// convenience operator to get pointer to 'i'-th frame
    { return &(frames[i]); };

  TotemRawVFATFrame* GetVFATFrameByID(int ID);			/// returns pointer to frame with ID 'ID' or NULL if there's no such a frame
  TotemRawVFATFrame* operator() (int ID)					/// convenience operator doing the same as method above
    { return GetVFATFrameByID(ID); };

  void CreateBufferList();								/// if you use TotemVFATFrameColl(1), you shoud you this method after inserting your VFAT frames into 'frames' array
  std::vector<unsigned short*> GetBufferList()			/// returns array of buffers ('data' fields) of all VFAT frames
    { return buffers; };

 private:
  std::vector<TotemRawVFATFrame> frames;
  std::vector<unsigned short*> buffers;
  __gnu_cxx::hash_map<int, int> map;
};

#endif

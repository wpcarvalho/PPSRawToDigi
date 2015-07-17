/**************************************************************
 * $Id: TotemVFATFrameColl.cc,v 1.1 2008/09/12 16:17:05 lgrzanka Exp $
 * $Revision: 1.1 $
 * $Date: 2008/09/12 16:17:05 $
 ***************************************************************/


#include "DataFormats/T1RawData/interface/TotemVFATFrameColl.h"

///----------------------------------------------------------------------------------------------------

TotemVFATFrameColl::TotemVFATFrameColl(char PlugVFATobjectsManualy)
{
  if (PlugVFATobjectsManualy) return;

  /// create frame and buffer arrays
  frames.resize(NUMBER_OF_VFATS);
  CreateBufferList();
}

///----------------------------------------------------------------------------------------------------

TotemVFATFrameColl::~TotemVFATFrameColl()
{
}

///----------------------------------------------------------------------------------------------------

void TotemVFATFrameColl::CreateBufferList()
{
  /// reset
  buffers.clear();

  /// add data of existing frames
  unsigned int i;
  for (i = 0; i < frames.size(); i++) {
    buffers.push_back(frames[i].getData());
  }

  ///  fill up by NULLs
  for (; i < NUMBER_OF_VFATS; i++) {
    buffers.push_back(NULL);
  }
}

///----------------------------------------------------------------------------------------------------

TotemRawVFATFrame* TotemVFATFrameColl::GetVFATFrameByID(int ID)
{
  /// if map is empty, create it
  /// in the map, channels are numbered from 1, it is because map[undefined ID] gives 0
  if (map.empty()) {
    for (unsigned int i = 0; i < frames.size(); i++) {
      map[frames[i].getChipID()] = i + 1;
    }
  }

  /// return the right frame
  if (map[ID]) return &(frames[map[ID] - 1]);
  else return NULL;
}

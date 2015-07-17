/*********************************************************
 * $Id: TotemRawVFATFrame.cc,v 1.1 2008/09/12 16:17:05 lgrzanka Exp $
 * $Revision: 1.1 $
 * $Date: 2008/09/12 16:17:05 $
 **********************************************************/

// Corresponds to data from one frame (data from one VFAT)
// Gives access to output from VFAT (flags,parameters,data section)

#include "DataFormats/T1RawData/interface/TotemRawVFATFrame.h"


//----------------------------------------------------------------------------------------------------

TotemRawVFATFrame::TotemRawVFATFrame(unsigned short *_data)
{
  if (_data) setData(_data);
}

//----------------------------------------------------------------------------------------------------

TotemRawVFATFrame::~TotemRawVFATFrame()
{
}


//----------------------------------------------------------------------------------------------------

bool TotemRawVFATFrame::channelActive(unsigned char channel)
{
  return ( data[1 + (channel / 16)] & (1 << (15 - channel % 16)) ) ? 1 : 0;
}

//----------------------------------------------------------------------------------------------------

std::vector<unsigned char> TotemRawVFATFrame::getActiveChannels()
{
  std::vector<unsigned char> channels;

  for (int i = 0; i < 8; i++) {
    /// quick check
    if (!data[1 + i]) continue;

    /// go throug bits
    unsigned short mask;
    char offset;
    for (mask = 1 << 15, offset = 15; mask; mask >>= 1, offset--) {
      //		for (mask = 1 << 15, offset = 0; mask; mask >>= 1, offset++) {
      if (data[1 + i] & mask) channels.push_back( i * 16 + offset );
    }
  }

  return channels;
}

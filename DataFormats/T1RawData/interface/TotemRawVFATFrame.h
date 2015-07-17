/*********************************************************
 * $Id: TotemRawVFATFrame.h,v 1.1 2008/09/12 16:17:05 lgrzanka Exp $
 * $Revision: 1.1 $
 * $Date: 2008/09/12 16:17:05 $
 **********************************************************/

#ifndef _TotemRawVFATFrame_h_
#define _TotemRawVFATFrame_h_

#include <vector>
#include <cstring>

/**
   Corresponds to data from one frame (data from one VFAT).
   The buffer is organized as following (reversed Figure 8 at page 23 of VFAT2 manual):

   buffer index	content			size
   ---------------------------------------------------------------
   0			CRC				16 bits
   1->8		Channel data	128 bits, channel 0 first
   9			ChipID			4 dummy bits (1110) + 12 bits
   10			EC, Flags		4 dummy bits (1100) + 8, 4 bits 
   11			BC				4 dummy bits (1010) + 12 bits
**/

class TotemRawVFATFrame
{
 public:
  TotemRawVFATFrame(unsigned short* _data = 0);		/// Constructor - copy data to private buffer
  TotemRawVFATFrame(const TotemRawVFATFrame& copy)	/// Copy constructor
    { this->setData(copy.data); };

  ~TotemRawVFATFrame();								/// Destructor - does nothing, actually

  void setData(const unsigned short *data);			/// copy data
  unsigned short* getData();							/// Return pointer to Channel Data (16 bytes)

  unsigned short getBC();								/// Return Bunch Crossing number (BC<11:0>)
  unsigned short getEC();								/// Return Event Counter (EV<7:0>)
  unsigned short getFlags();							/// Return flags
  unsigned short getChipID();							/// Return ChipID (ChipID<11:0>)
  unsigned short getCRC();							/// Return Check sum

  bool channelActive(unsigned char channel);			/// Check if channel number 'channel' was active
  /// 	returns positive number if it was active, 0 otherwise
  std::vector<unsigned char> getActiveChannels();		/// Returns array of active channels
  /// 	it's MUCH MORE EFFICIENT than the previous method for events with low channel occupancy

 private:
  unsigned short data[12];							/// data buffer
};

///----------------------------------------------------------------------------------------------------

inline void TotemRawVFATFrame::setData(const unsigned short *_data)
{
  memcpy(data, _data, 24);
}

///----------------------------------------------------------------------------------------------------

inline unsigned short* TotemRawVFATFrame::getData()
{
  return data;
}

///----------------------------------------------------------------------------------------------------

inline unsigned short TotemRawVFATFrame::getBC()
{
  return data[11] & 0x0FFF;
}

///----------------------------------------------------------------------------------------------------

inline unsigned short TotemRawVFATFrame::getEC()
{
  return (data[10] & 0x0FF0) >> 4;
}

///----------------------------------------------------------------------------------------------------

inline unsigned short TotemRawVFATFrame::getFlags()
{
  return data[10] & 0x000F;
}

///----------------------------------------------------------------------------------------------------

inline unsigned short TotemRawVFATFrame::getChipID()
{
  return data[9] & 0x0FFF;
}

///----------------------------------------------------------------------------------------------------

inline unsigned short TotemRawVFATFrame::getCRC()
{
  return data[0];
}

#endif

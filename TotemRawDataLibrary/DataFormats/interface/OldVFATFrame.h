/****************************************************************************
*
* This is a part of the TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com)
*   Leszek Grzanka 
*
****************************************************************************/

#ifndef _Totem_OldVFATFrame_h_
#define _Totem_OldVFATFrame_h_

#include <vector>
#include "TotemRawDataLibrary/DataFormats/interface/VFATFrame.h"
#include <cstring>

namespace Totem {

/**
 * \brief VFAT frame class.
 * 
 * Corresponds to data from one frame (data from one VFAT). The buffer is organized as following (reversed
 * Figure 8 at page 23 of VFAT2 manual):
 * \verbatim
 * buffer index  content      size
 * ---------------------------------------------------------------
 *   0      CRC        16 bits
 *   1->8    Channel data  128 bits, channel 0 first
 *   9      ChipID      4 dummy bits (1110) + 12 bits
 *   10      EC, Flags    4 dummy bits (1100) + 8, 4 bits 
 *   11      BC        4 dummy bits (1010) + 12 bits
 * \endverbatim
 *
 * \ingroup TotemRawDataLibrary
**/
class OldVFATFrame: public VFATFrame
{
  public:

    OldVFATFrame(word* _data = NULL);                        ///< Constructor - copy data to private buffer
    OldVFATFrame(const OldVFATFrame& copy)                      ///  Copy constructor
      { setData(copy.data); }

    ~OldVFATFrame();                                         ///< Destructor - does nothing, actually

    void setData(const word *data);                       ///< copy data
    const word* getData() const;                          ///< Returns pointer to Channel Data (16 bytes)

    word getBC() const;
    word getEC() const;
    word getFlags() const;
    word getChipID() const;
    word getCRC() const;                                  ///< Returns Check sum
    bool checkFootprint() const;
    bool checkCRC() const;
    bool checkBitShift() const;                           ///< Checks whether the frame is missing the first 1 bit.
    bool channelActive(unsigned char channel) const;
    std::vector<unsigned char> getActiveChannels() const;


    void Print(bool binary = false) const;

  protected:
    word data[12];                                        ///< data buffer

    /// internaly used to check CRC
    static word crc_calc(word crc_in, word dato);
};

//----------------------------------------------------------------------------------------------------

inline void OldVFATFrame::setData(const word *_data)
{
  memcpy(data, _data, 24);
}

//----------------------------------------------------------------------------------------------------

inline const VFATFrame::word* OldVFATFrame::getData() const
{
  return data;
}

//----------------------------------------------------------------------------------------------------

inline VFATFrame::word OldVFATFrame::getBC() const
{
  return data[11] & 0x0FFF;
}

//----------------------------------------------------------------------------------------------------

inline VFATFrame::word OldVFATFrame::getEC() const
{
  return (data[10] & 0x0FF0) >> 4;
}

//----------------------------------------------------------------------------------------------------

inline VFATFrame::word OldVFATFrame::getFlags() const
{
  return data[10] & 0x000F;
}

//----------------------------------------------------------------------------------------------------

inline VFATFrame::word OldVFATFrame::getChipID() const
{
  return data[9] & 0x0FFF;
}

//----------------------------------------------------------------------------------------------------

inline bool OldVFATFrame::checkFootprint() const
{
  return ((data[9] & 0xF000) == 0xE000) && ((data[10] & 0xF000) == 0xC000) && ((data[11] & 0xF000) == 0xA000);
}

//----------------------------------------------------------------------------------------------------

inline VFATFrame::word OldVFATFrame::getCRC() const
{
  return data[0];
}

} // namespace
#endif


/*********************************************************
* Seyed Mohsen Etesami
**********************************************************/

#ifndef EventFilter_CTPPSRawToDigi_DiamondVFATFrame
#define EventFilter_CTPPSRawToDigi_DiamondVFATFrame

#include <vector>
#include <cstddef>
#include <stdint.h>

/**
 * Representation of VFAT frame plus extra info added by DAQ.
**/
class DiamondVFATFrame
{
  public:
    typedef uint16_t word;
    typedef uint32_t timeinfo;

  public:
    DiamondVFATFrame(const word* _data = NULL);

    DiamondVFATFrame(const DiamondVFATFrame& copy)
    {
      setDataAndCast(copy.data);
      presenceFlags = copy.presenceFlags;
    }

    virtual ~DiamondVFATFrame() {}

     DiamondVFATFrame::timeinfo getLeadingEtime() const;
     DiamondVFATFrame::timeinfo getTrailingEtime()          const;
     DiamondVFATFrame::timeinfo getThresholdVolt() const;
//     DiamondVFATFrame::timeinfo getMultihit() const;
//     DiamondVFATFrame::timeinfo getHptdcerrorflag() const;
   
    /// Copies a memory block to data buffer.
    void setDataAndCast(const word *_data);

    DiamondVFATFrame::word* getData()
    {
      return data;
    }

    /// Returns Bunch Crossing number (BC<11:0>).
    DiamondVFATFrame::word getBC() const
    {
      return data[11] & 0x0FFF;
    }

    /// Returns Event Counter (EV<7:0>).
    DiamondVFATFrame::word getEC() const
    {
      return (data[10] & 0x0FF0) >> 4;
    }

    /// Returns flags.
    DiamondVFATFrame::word getFlags() const
    {
      return data[10] & 0x000F;
    }

    /// Returns vfatID .
    DiamondVFATFrame::word getChannelID() const
    {
      return data[9] & 0x0FFF;
    }

    DiamondVFATFrame::word getMultihit() const
    {
      return data[2] & 0x01;
    }

    DiamondVFATFrame::word getHptdcerrorflag() const
    {
      return data[1] & 0xFFFF;
    }

    /// Returns the CRC.
    DiamondVFATFrame::word getCRC() const
    { 
      return data[0];
    }

    /// Sets presence flags.
    void setPresenceFlags(uint8_t v)
    {
      presenceFlags = v;
    }

    /// Returns true if the BC word is present in the frame.
    bool isBCPresent() const
    {
      return presenceFlags & 0x1;
    }

    /// Returns true if the EC word is present in the frame.
    bool isECPresent() const
    {
      return presenceFlags & 0x2;
    }

    /// Returns true if the ID word is present in the frame.
    bool isIDPresent() const
    {
      return presenceFlags & 0x4;
    }

   /// Returns true if the leading edge time  word is present in the frame.
    bool isLEDTimePresent() const
    {
      return presenceFlags & 0x8;
    }

      /// Returns true if the trainling edge time  word is present in the frame.
    bool isTEDTimePresent() const
    {
      return presenceFlags & 0x10;
    }
     /// Returns true if the threshold voltage  word is present in the frame.
    bool isThVolPresent() const
    {
      return presenceFlags & 0x20;
    }
 
   /// Returns true if the multi hit  word is present in the frame.
    bool isMuHitPresent() const
    {
      return presenceFlags & 0x40;
    }

    /// Checks the fixed bits in the frame.
    /// Returns false if any of the groups (in BC, EC and ID words) is present but wrong.
    bool checkFootprint() const;

    /// Returns false if any of the groups (in LEDTime, TEDTime, Threshold Voltage and Multi hit  words) is present but wrong.
    bool checkTimeinfo() const;

    /// Prints the frame.
    /// If binary is true, binary format is used.
    void Print(bool binary = false) const;

  private:
    /** Raw data frame as sent by electronics.
    * The container is organized as follows:
    * \verbatim
    * buffer index   content       size
    * ---------------------------------------------------------------
    *   0            CRC                    16 bits
    *   1            HPTDC error            16 bits
    *   2            Multi hit              5 constant bits (01111) 10 empty bits + 1 bit
    *   3,4          Threshold Voltage      5 constant bits (01110) +27 bits
    *   5,6          Trailing Etime         5 constant bits (01101) +27 bits
    *   7,8          Leading Etime          5 constant bits (01100) +27 bits
    *   9            ChipID                 4 constant bits (1110) + 12 bits
    *   10           EC, Flags              4 constant bits (1100) + 8, 4 bits
    *   11           BC                     4 constant bits (1010) + 12 bits
    * \endverbatim
    **/
    word data[12];
    word semidata[12];
    timeinfo timedata[6];

    /// Flag indicating the presence of various components.
    ///   bit 1: "BC word" (buffer index 11)
    ///   bit 2: "EC word" (buffer index 10)
    ///   bit 3: "ID word" (buffer index 9)
    ///   bit 4: "LEDTime" (buffer index 7)
    ///   bit 5: "TEDTime" (bufer index 5)
    ///   bit 6: "Thresholdvoltage" (bufer index 3)
    ///   bit 7: "Multi hit" (bufer index 2)

    uint8_t presenceFlags;

};                                                                     

#endif

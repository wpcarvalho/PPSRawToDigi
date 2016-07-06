/****************************************************************************
*  Seyed Mohsen Etesami   
****************************************************************************/

#ifndef DataFormats_CTPPSDigi_DiamondVFATStatus
#define DataFormats_CTPPSDigi_DiamondVFATStatus

#include <bitset>
#include <map>

#include "DataFormats/CTPPSDigi/interface/HPTDCErrorFlags.h"

//----------------------------------------------------------------------------------------------------

/**
 * Class which contains information about conversion from RAW to DIGI for a single read-out (VFAT).
 */
class DiamondVFATStatus
{
  public:
    DiamondVFATStatus() : status(0), hptdc_error_status(0) {}

    typedef std::bitset<8> DiamondVFATStatusWord;

    /// VFAT is present in mapping but no data is present int raw event
    inline bool isMissing() const { return status[0]; }

    /// 12-bit hw id from the header of the vfat frame is diffrent from the 16-bit one from hw mapping.
    inline bool isIDMismatch() const { return status[1]; }
    
    /// Footprint error
    inline bool isFootprintError() const { return status[2]; }

    /// CRC error
    inline bool isCRCError() const { return status[3]; }

    /// VFATFrame event number doesn't follow the number derived from DAQ
    inline bool isECProgressError() const { return status[4]; }

    /// BC number is incorrect
    inline bool isBCProgressError() const { return status[5]; }

    /// All channels from that VFAT are not taken into account
    inline bool isFullyMaskedOut() const { return status[6]; }

    /// Some channels from VFAT ale masked out, but not all
    inline bool isPartiallyMaskedOut() const { return status[7]; }

    /// No channels are masked out
    inline bool isNotMasked() const { return !(status[6] || status[7]); }

    
    inline void setMissing(bool val = true) { status[0] = val; }
    inline void setIDMismatch(bool val = true) { status[1] = val; }
    inline void setFootprintError(bool val = true) { status[2] = val; }
    inline void setCRCError(bool val = true) { status[3] = val; }
    inline void setECProgressError(bool val = true) { status[4] = val; }
    inline void setBCProgressError(bool val = true) { status[5] = val; }

    inline void setFullyMaskedOut() { status[6]=true; }
    inline void setPartiallyMaskedOut() { status[7]=true; }
    inline void setNotMasked() { status[6]=status[7]=false; }


    /**  HPTDC Error flags inside each VFAT Raw data frame.
     * \verbatim
     *   15 0      
     *	14 Hit lost in group 0 from read-out fifo overflow.													
     *	13 Hit lost in group 0 from L1 buffer overflow													
     *	12 Hit error have been detected in group 0.													
     *	11 Hit lost in group 1 from read-out fifo overflow.													
     *	10 Hit lost in group 1 from L1 buffer overflow													
     *	9  Hit error have been detected in group 1.													
     *	8  Hit data lost in group 2 from read-out fifo overflow.													
     *	7  Hit lost in group 2 from L1 buffer overflow													
     *	6  Hit error have been detected in group 2.													
     *	5  Hit lost in group 3 from read-out fifo overflow.													
     *	4  Hit lost in group 3 from L1 buffer overflow													
     *	3  Hit error have been detected in group 3.													
     *	2  Hits rejected because of programmed event size limit													
     *	1  Event lost (trigger fifo overflow).													
     *   0  Internal fatal chip error has been detected.													
     * \endverbatim
     **/

    inline void setHPTDCErrors(uint16_t hpdcer) {
      hptdc_error_status = HPTDCErrorFlags(hpdcer);
    }

    inline bool  getHPTDCErrors( int errornum_)
    {
      return hptdc_error_status.ErrorId(errornum_);
    }

    bool isOK() const
    {
      return !(status[0] || status[1] || status[2] || status[3] || status[4] || status[5]);
    }

    bool operator < (const DiamondVFATStatus &cmp) const
    {
      return (status.to_ulong() < cmp.status.to_ulong());
    }
  
    friend std::ostream& operator << (std::ostream& s, const DiamondVFATStatus &st);

  private:

    /// the status bits
    DiamondVFATStatusWord status;
    HPTDCErrorFlags hptdc_error_status;

};

#endif

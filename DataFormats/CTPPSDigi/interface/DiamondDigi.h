#ifndef CTPPSDigi_DiamondDigi_h
#define CTPPSDigi_DiamondDigi_h

/** \class DiamondDigi
 *
 * Digi for PPSTiming.
 *  
 *
 * \author Seyed Mohsen Etesami
 * March 2016
 */

#include <boost/cstdint.hpp>
#include "DataFormats/CTPPSDigi/interface/HPTDCErrorFlags.h"

class DiamondDigi{

 public:
  
  DiamondDigi(unsigned short chid_, unsigned int ledgt_, unsigned int tedgt_, unsigned int threvolt, bool mhit_, unsigned short hptdcerror_);
  DiamondDigi();
  ~DiamondDigi() {};
  
  /// Digis are equal if they are have same chid, ledt and tedt, threshold voltage, multihit flag, hptdcerror flags
  bool operator==(const DiamondDigi& digi) const;

  /// Return digi values number
  unsigned short getChannelId() const { return chid; }
  unsigned int getLeadingEdge() const { return ledgt; }
  unsigned int getTrailingEdge() const { return tedgt; }
  unsigned int getThresholdVoltage() const { return threvolt; }
  bool getMultipleHit() const { return mhit; }
  HPTDCErrorFlags getHPTDCErrorFlags() const { return hptdcerror; }

  /// Set digi values
  inline void setChannelId(unsigned char chid_) { chid = chid_; }
  inline void setLeadingEdge(unsigned int ledgt_) { ledgt = ledgt_; }
  inline void setTrailingEdge(unsigned int tedgt_) { tedgt = tedgt_; }
  inline void setThresholdVoltage(unsigned int threvolt_) { threvolt = threvolt_; }
  inline void setMultipleHit(bool mhit_) { mhit = mhit_; }
  inline void setHPTDCErrorFlags(const HPTDCErrorFlags& hptdcerror_) { hptdcerror = hptdcerror_; }

  /// Print content of digi
  //void print() const;

 private:
  unsigned char chid;
  unsigned int ledgt;
  unsigned int tedgt;
  unsigned int threvolt;
  bool mhit;
  HPTDCErrorFlags hptdcerror;
};

inline bool operator< (const DiamondDigi& one, const DiamondDigi& other)
{
  return one.getChannelId() < other.getChannelId();
}

#include <iostream>

inline std::ostream & operator<<(std::ostream & o, const DiamondDigi& digi)
{
  return o << " " << digi.getChannelId()
           << " " << digi.getLeadingEdge()
	   << " " << digi.getTrailingEdge()
           << " " << digi.getThresholdVoltage()
           << " " << digi.getMultipleHit()
           << " " << digi.getHPTDCErrorFlags().error_flags;
}

#endif


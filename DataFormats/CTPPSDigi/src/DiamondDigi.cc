/** \file
 * 
 *
 * \author Seyed Mohsen Etesami
 */

#include <DataFormats/CTPPSDigi/interface/DiamondDigi.h>

using namespace std;

DiamondDigi::DiamondDigi(unsigned short chid_, unsigned int ledgt_, unsigned int tedgt_, unsigned int threvolt_, bool mhit_, unsigned short hptdcerror_) :
  chid(chid_), ledgt(ledgt_), tedgt(tedgt_), threvolt(threvolt_), mhit(mhit_), hptdcerror(HPTDCErrorFlags(hptdcerror_))
{}

DiamondDigi::DiamondDigi() :
  chid(0), ledgt(0), tedgt(0), threvolt(0), mhit(false), hptdcerror(HPTDCErrorFlags(0))
{}

// Comparison
bool
DiamondDigi::operator==(const DiamondDigi& digi) const
{
  if (chid      != digi.getChannelId() 
   || ledgt     != digi.getLeadingEdge()
   || tedgt     != digi.getTrailingEdge()
   || threvolt  != digi.getThresholdVoltage()
   || mhit      != digi.getMultipleHit()
   || hptdcerror.error_flags != digi.getHPTDCErrorFlags().error_flags) return false;
  return true; 
}


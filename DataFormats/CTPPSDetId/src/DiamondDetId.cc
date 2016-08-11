/****************************************************************************
 * Author: Seyed Mohsen Etesami
 ****************************************************************************/


#include "DataFormats/CTPPSDetId/interface/DiamondDetId.h"
#include "FWCore/Utilities/interface/Exception.h"

using namespace std;


DiamondDetId::DiamondDetId():DetId(DetId::VeryForward,pps_timing_diamond_subdet_id)  //first constructor and defult
{}


DiamondDetId::DiamondDetId(uint32_t id):DetId(id)   //second constructor
{
  if (! Check(id))
    {
      throw cms::Exception("InvalidDetId") << "PPSTimingDetId ctor:"
					   << " det: " << det()
					   << " subdet: " << subdetId()
					   << " is not a valid PPS Timing id";  
    }
}


void DiamondDetId::init(unsigned int Arm, unsigned int Station,
			unsigned int RomanPot,unsigned int Plane, unsigned int Detector)
{
  if( Arm>=2 || Station>=3 || RomanPot>=1 || Plane>3 || Detector>16)   //in timing station=1 we have one RP for timing therefore RomanPot>=1
    {
      throw cms::Exception("InvalidDetId") << "PPSTimingDetId ctor:" 
					   << " Invalid parameters: " 
					   << " Arm "<<Arm
					   << " Station "<<Station
					   << " RomanPot "<<RomanPot
                                           << " PLane "<<Plane
					   << " Detector "<<Detector
					   << std::endl;
    }

  uint32_t ok=0xfe000000;
  id_ &= ok;
  //0x3f= 3X8+15=39
  id_ |= ((Arm&0x1) << startArmBit);
  id_ |= ((Station&0x3) << startStationBit);
  id_ |= ((RomanPot&0x7) << startRPBit);
  id_ |= ((Plane&0x3) << startPlaneBit);
  id_ |= ((Detector&0xF) << startDetBit);
}


DiamondDetId::DiamondDetId(unsigned int Arm, unsigned int Station,
			   unsigned int RomanPot,unsigned int Plane, unsigned int Detector):       //third constructor with full id of all parts
  DetId(DetId::VeryForward, pps_timing_diamond_subdet_id)
{
  this->init(Arm,Station,RomanPot,Plane,Detector);
}


std::ostream& operator<<( std::ostream& os, const DiamondDetId& id )
{
  os << " Arm "<<id.Arm()
     << " Station "<<id.Station()
     << " RomanPot "<<id.RomanPot()
     << " Palne "<<id.Plane()
     << " Detector "<<id.Detector();

  return os;
}




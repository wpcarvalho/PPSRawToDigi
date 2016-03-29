/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors: 
 *	Hubert Niewiadomski
 *	Jan Kašpar (jan.kaspar@gmail.com) 
 *
 ****************************************************************************/


#include "DataFormats/TotemRPDetId/interface/TotemRPDetId.h"
#include "FWCore/Utilities/interface/Exception.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

TotemRPDetId::TotemRPDetId():DetId(DetId::VeryForward, totem_rp_subdet_id)
{
}

//----------------------------------------------------------------------------------------------------

TotemRPDetId::TotemRPDetId(uint32_t id):DetId(id)
{
  if (! Check(id))
    {
      throw cms::Exception("InvalidDetId") << "TotemRPDetId ctor:"
					   << " det: " << det()
					   << " subdet: " << subdetId()
					   << " is not a valid Totem RP id";  
    }
}

//----------------------------------------------------------------------------------------------------

void TotemRPDetId::init(unsigned int Arm, unsigned int Station, unsigned int RomanPot, unsigned int Detector)
{
  if( Arm>=2 || Station>=3 || RomanPot>=6 || Detector>=10)
    {
      throw cms::Exception("InvalidDetId") << "TotemRPDetId ctor:" 
					   << " Invalid parameters: " 
					   << " Arm "<<Arm
					   << " Station "<<Station
					   << " RomanPot "<<RomanPot
					   << " Detector "<<Detector
					   << std::endl;
    }

  uint32_t ok=0xfe000000;
  id_ &= ok;

  id_ |= ((Arm&0x1) << startArmBit);
  id_ |= ((Station&0x3) << startStationBit);
  id_ |= ((RomanPot&0x7) << startRPBit);
  id_ |= ((Detector&0xf) << startDetBit);
}

//----------------------------------------------------------------------------------------------------

TotemRPDetId::TotemRPDetId(unsigned int Arm, unsigned int Station, unsigned int RomanPot, unsigned int Detector):       
  DetId(DetId::VeryForward,totem_rp_subdet_id)
{
  this->init(Arm,Station,RomanPot,Detector);
}

//----------------------------------------------------------------------------------------------------

std::ostream& operator<<( std::ostream& os, const TotemRPDetId& id )
{
  os << " Arm "<<id.Arm()
     << " Station "<<id.Station()
     << " RomanPot "<<id.RomanPot()
     << " Detector "<<id.Detector();

  return os;
}

//----------------------------------------------------------------------------------------------------

string TotemRPDetId::SystemName(NameFlag flag)
{
  string name;
  if (flag == nFull) name = "rp";
  if (flag == nPath) name = "RP";

  return name;
}

//----------------------------------------------------------------------------------------------------

string TotemRPDetId::ArmName(unsigned int id, NameFlag flag)
{
  string name;
  if (flag == nFull) name = SystemName(flag) + "_";
  if (flag == nPath) name = SystemName(flag) + "/sector ";

  if (id == 0) name += "45";
  if (id == 1) name += "56";

  return name;
}

//----------------------------------------------------------------------------------------------------

string TotemRPDetId::StationName(unsigned int id, NameFlag flag)
{
  string name;
  if (flag == nFull) name = ArmName(id / 10, flag) + "_";
  if (flag == nPath) name = ArmName(id / 10, flag) + "/station ";

  if ((id % 10) == 0) name += "210";
  if ((id % 10) == 2) name += "220";

  return name;
}

//----------------------------------------------------------------------------------------------------

string TotemRPDetId::RPName(unsigned int id, NameFlag flag)
{
  string name; 
  if (flag == nFull) name = StationName(id / 10, flag) + "_";
  if (flag == nPath) name = StationName(id / 10, flag) + "/";

  if ((id % 10) == 0) name += "nr_tp";
  if ((id % 10) == 1) name += "nr_bt";
  if ((id % 10) == 2) name += "nr_hr";
  if ((id % 10) == 3) name += "fr_hr";
  if ((id % 10) == 4) name += "fr_tp";
  if ((id % 10) == 5) name += "fr_bt";
  return name;
}

//----------------------------------------------------------------------------------------------------

string TotemRPDetId::PlaneName(unsigned int id, NameFlag flag)
{
  string name;
  if (flag == nFull) name = RPName(id / 10, flag) + "_";
  if (flag == nPath) name = RPName(id / 10, flag) + "/plane ";

  char buf[10];
  sprintf(buf, "%02u", (id % 10) + 1);

  return name + buf;
}

//----------------------------------------------------------------------------------------------------

string TotemRPDetId::ChipName(unsigned int id, NameFlag flag)
{
  string name;
  if (flag == nFull) name = PlaneName(id / 10, flag) + "_";
  if (flag == nPath) name = PlaneName(id / 10, flag) + "/chip ";

  char buf[10];
  sprintf(buf, "%u", (id % 10) + 1);

  return name + buf;
}

//----------------------------------------------------------------------------------------------------

string TotemRPDetId::StripName(unsigned int id, unsigned char strip, NameFlag flag)
{
  string name;
  if (flag == nFull) name = ChipName(id, flag) + "_";
  if (flag == nPath) name = ChipName(id, flag) + "/strip";

  char buf[10];
  sprintf(buf, "%u", strip);

  return name + buf;
}

//----------------------------------------------------------------------------------------------------

string TotemRPDetId::OfficialName(ElementLevel level, unsigned int id, NameFlag flag, unsigned char strip)
{
  switch (level) {
    case lSystem: return SystemName(flag);
    case lArm: return ArmName(id, flag);
    case lStation: return StationName(id, flag);
    case lRP: return RPName(id, flag);
    case lPlane: return PlaneName(id, flag);
    case lChip: return ChipName(id, flag);
    case lStrip: return StripName(id, flag);
    default: return "";
  }
}



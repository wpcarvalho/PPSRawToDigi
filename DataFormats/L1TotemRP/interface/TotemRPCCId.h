/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *  Leszek Grzanka (braciszek@gmail.com)
 *
 ****************************************************************************/

#ifndef DataFormats_L1TotemRP_TotemRPCCId
#define DataFormats_L1TotemRP_TotemRPCCId

#include "DataFormats/DetId/interface/DetId.h"

#include <iosfwd>
#include <iostream>

/**
 * TODO: describe
 **/
class TotemRPCCId: public DetId
{
  public:
    TotemRPCCId();

    /// Construct from a packed id. It is required that the Detector part of
    /// id is Totem and the SubDet part is RP, otherwise an exception is thrown.
    explicit TotemRPCCId(uint32_t id);

    /// Construct from fully qualified identifier.
    TotemRPCCId(unsigned int Arm, unsigned int Station, unsigned int RomanPot, unsigned int Direction);

    /// Bit 24 = Arm: 1=z>0 0=z<0
    /// Bits [22:23] Station
    /// Bits [19:21] Roman Pot number
    /// Bits [15:18] CC direction: 0 - V, 1 - U
    inline int arm() const
    {
      return int((id_>>startArmBit) & 0x1);
    }

    inline int station() const
    {
      return int((id_>>startStationBit) & 0x3);
    }

    inline int romanPot() const
    {
      return int((id_>>startRPBit) & 0x7);
    }

    inline int direction() const
    {
      return int((id_>>startDirBit) & 0xf);
    }

    inline bool isStripsCoordinateUDirection() const
    {
      return direction()%2;
    }

    inline bool IsStripsCoordinateVDirection() const
    {
      return !isStripsCoordinateUDirection();
    }

    static const int startArmBit = 24;
    static const int startStationBit = 22;
    static const int startRPBit = 19;
    static const int startDirBit = 15;
    static const int totem_rp_subdet_id = 3;

  private:
    inline void init(unsigned int Arm, unsigned int Station, unsigned int RomanPot, unsigned int Direction);
};

std::ostream& operator<<(std::ostream& os, const TotemRPCCId& id);

#endif

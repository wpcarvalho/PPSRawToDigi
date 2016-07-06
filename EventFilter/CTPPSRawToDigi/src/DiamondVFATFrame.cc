/****************************************************************************
*  Seyed Mohsen Etesami
****************************************************************************/

#include "EventFilter/CTPPSRawToDigi/interface/DiamondVFATFrame.h"

#include <stdio.h>
#include <cstring>
#include <iostream>

using namespace std;

//----------------------------------------------------------------------------------------------------

DiamondVFATFrame::DiamondVFATFrame(const DiamondVFATFrame::word *_data) :
  presenceFlags(127)    // by default BC, EC, ID, LEDTime, TEDtime,ThreVolt and MuHit are present
{
  if (_data) setDataAndCast(_data);
  else memset(data, 0, 12 * sizeof(word));
}

//----------------------------------------------------------------------------------------------------

void
DiamondVFATFrame::setDataAndCast(const DiamondVFATFrame::word *_data)
{
  memcpy(data, _data, 24);
  memcpy(semidata, _data, 24);
  memmove (semidata,semidata+3,18);
  DiamondVFATFrame::timeinfo* timedata_ = (DiamondVFATFrame::timeinfo*) semidata;
  memcpy(timedata, timedata_, 24);
}

//----------------------------------------------------------------------------------------------------

DiamondVFATFrame::timeinfo
DiamondVFATFrame::getLeadingEtime() const
{
  //return (timedata[2] >> 11) & 0x1FFFFF;
  return ((data[8]&0x1f)<<16)+data[7]; //FW bug
}

DiamondVFATFrame::timeinfo
DiamondVFATFrame::getTrailingEtime() const
{
  //return (timedata[1] >> 11) & 0x1FFFFF;
  return ((data[6]&0x1f)<<16)+data[5]; //FW bug
}

DiamondVFATFrame::timeinfo
DiamondVFATFrame::getThresholdVolt() const
{
  //return (timedata[0] >>5 ) & 0x7FFFFFF;
  return ((data[4]&0x3ff)<<16)+data[3]; //FW bug
}


//----------------------------------------------------------------------------------------------------

bool
DiamondVFATFrame::checkTimeinfo() const
{
  if (isLEDTimePresent() && (data[7] & 0xF800) != 0x6000)
    return false;

  if (isTEDTimePresent() && (data[5] & 0xF800) != 0x6800)
    return false;

  if (isThVolPresent() && (data[3] & 0xF800) != 0x7000)
    return false;

  if (isMuHitPresent() && (data[2] & 0xF800) != 0x7800)
    return false;

  return true;
}

//----------------------------------------------------------------------------------------------------

bool
DiamondVFATFrame::checkFootprint() const
{
  if (isIDPresent() && (data[9] & 0xF000) != 0xE000)
    return false;

  if (isECPresent() && (data[10] & 0xF000) != 0xC000)
    return false;

  if (isBCPresent() && (data[11] & 0xF000) != 0xA000)
    return false;

  return true;
}

//----------------------------------------------------------------------------------------------------

void
DiamondVFATFrame::Print(bool binary) const
{
  if (binary)
  {
    for (int i = 0; i < 12; i++)
    {
      const word &w = data[11 - i];
      word mask = (1 << 15);
      for (int j = 0; j < 16; j++)
      {
        if (w & mask)
          printf("1");
        else
          printf("0");
        mask = (mask >> 1);
        if ((j + 1) % 4 == 0)
          printf("|");
      }
      printf("\n");
    }
  } else {
    printf("ID = %03x, BC = %04u, EC = %03u, flags = %2u, CRC = %04x ", getChannelID(), getBC(), getEC(), getFlags(), getCRC());


    if (checkFootprint())
      printf("  OK");
    else
      printf("FAIL");

    printf(", frame = %04x|%04x|%04x|", data[11], data[10], data[9]);
    for (int i = 8; i > 0; i--)
      printf("%04x", data[i]);
    printf("|%04x", data[0]);

    printf(", presFl=%x", presenceFlags);

    printf("\n");
  }
}

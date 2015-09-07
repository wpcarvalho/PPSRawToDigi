/**********************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors:
*   Michal Zmuda (m.zmuda@aol.com)
*
**********************************************************/


#include <stdio.h>
#include <algorithm>

#include "TotemRawDataLibrary/DataFormats/interface/CommonDef.h"
#include "TotemRawDataLibrary/DataFormats/interface/RawVFATFrame.h"

using namespace std;

// improves readability when having an array improves memory allocation
//      single short has 2B but much more is allocated - even 8B for 64bit OS
//      an array of 5 shorts would then occupy 16B (5x2b ceil) instead of  
//      40B for five single short fields - this is 60% of gain
//      in context of whole frame it depends (there is an short (8B allocated) per active channel)
#define header      infoHolder[0]
#define tailer      infoHolder[1]
#define wordWithEC  infoHolder[2]
#define wordWithBC  infoHolder[3]
#define wordWithID  infoHolder[4]

namespace Totem {


RawVFATFrame::RawVFATFrame(word data[14])
{
    for (unsigned int i = 0; i < 5; i++)
      infoHolder[i] = 0;

    header = data[0];
    if((data[1] & 0xFF00) == 0xF000){
        tailer = data[1];
    } else {
        wordWithBC = data[1];
        wordWithEC = data[2];
        wordWithID = data[3];
        tailer = data[12];

        //active channels
        for (int channel = 0; channel < 128; channel++) {
            if ((data[4 + (channel / 16)] & (1 << (channel % 16)))) {
                activeChannels.push_back(channel);
            }
        }
    }
}

RawVFATFrame::~RawVFATFrame(){}

bool RawVFATFrame::channelActive(unsigned char channel) const{
    return find(activeChannels.begin(), activeChannels.end(), channel) != activeChannels.end();
}

std::vector<unsigned char> RawVFATFrame::getActiveChannels() const{
    //todo consider changes as it exposes internal structure of class (maybe const for vector?)
    return activeChannels;
}

VFATFrame::word RawVFATFrame::getBC() const
{
  return wordWithBC & 0x0FFF;
}

VFATFrame::word RawVFATFrame::getEC() const
{
  return (wordWithEC & 0x0FF0) >> 4;
}
VFATFrame::word RawVFATFrame::getFlags() const
{
  return wordWithEC & 0x000F;
}

VFATFrame::word RawVFATFrame::getChipID() const
{
  return wordWithID & 0x0FFF;
}

VFATFrame::word RawVFATFrame::getSize() const
{
  return tailer & 0x00FF;
}

bool RawVFATFrame::checkCRC() const
{
  return true;//no crc at all
}

bool RawVFATFrame::checkFootprint() const
{
  return ((tailer & 0xFF00) == 0xF000)
    && ((wordWithID & 0xF000) == 0xE000)
    && ((wordWithEC & 0xF000) == 0xC000)
    && ((wordWithBC & 0xF000) == 0xA000)
    && ((header & 0xF000) == 0x9000);
}

VFATFrame::word RawVFATFrame::getRxFlags() const
{
  return (header & 0x0F00) >> 8;
}

VFATFrame::word RawVFATFrame::getFiberIdx() const
{
  return (header & 0x00F0) >> 4;
}

VFATFrame::word RawVFATFrame::getGohIdx() const{
  return header & 0x000F;
}

bool RawVFATFrame::passedHardwareSynchronisationChecks() const{
  return  !(getRxFlags() & 1);
}

void printWord(unsigned short word){
    unsigned short mask = (1 << 15);
    for (int j = 0; j < 16; j++) {
        printf((word & mask) ? "1" : "0");
        mask = (mask >> 1);
        if ((j + 1) % 4 == 0) printf("|");
    }
    printf("\n");
}

void RawVFATFrame::Print(bool binary) const{
    if (binary) {
        //header
        printWord(header);
        printWord(wordWithBC);
        printWord(wordWithEC);
        printWord(wordWithID);
        //data
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 16; j++) {
                printf(channelActive(i*16 + (15-j)) ? "1" : "0");
                if ((j + 1) % 4 == 0) printf("|");
            }
            printf("\n");
        }
        //tailer
        printWord(tailer);
    } else {
        printf("ID = %03x, BC = %04u, EC = %03u, flags = %2u, ", getChipID(), getBC(), getEC(), getFlags());
        printf("RxFlags = %2u, FiberIdx = %2u, GohIdx = %2u ", getRxFlags(), getFiberIdx(), getGohIdx());
        printf(checkFootprint() ? "footprint   OK" : "footprint FAIL");
        //header
        printf(", frame = %04x|%04x|%04x|%04x|", header, wordWithBC, wordWithEC, wordWithID);
        //data
        for (int i = 8; i > 0; i--){
            unsigned short word = 0;
            for (int j = 0; j < 16; j++) {
                if(channelActive(i*16 + j)){
                    word = word | (1<<j);
                }
            }
            printf("%04x", word);
        }
        //tailer
        printf("|%04x\n", tailer);
    }
}


}

#undef header
#undef tailer
#undef wordWithEC
#undef wordWithBC
#undef wordWithID

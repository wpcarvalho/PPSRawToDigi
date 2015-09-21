/**********************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors:
*   Michal Zmuda (m.zmuda@aol.com)
*
**********************************************************/


#include <cstdio>
#include <cstring>
#include <algorithm>

#include "TotemRawDataLibrary/DataFormats/interface/CommonDef.h"
#include "TotemRawDataLibrary/DataFormats/interface/ClusterizationVFATFrame.h"

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

const int BC_FOOTPRINT = 0xA000;
const int EC_FOOTPRINT = 0xC000;
const int ID_FOOTPRINT = 0xE000;

namespace Totem {

ClusterizationVFATFrame::ClusterizationVFATFrame(word* data)
{
    for (unsigned int i = 0; i < 5; i++)
      infoHolder[i] = 0;

    dataHolder = NULL;
    dataSize = 0;

    unsigned cnt = 0;

    header = data[cnt];
    cnt++;

    if((data[cnt] & 0xF000) == BC_FOOTPRINT){
        wordWithBC = data[cnt];
        cnt++;
    }
    if((data[cnt] & 0xF000) == EC_FOOTPRINT){
        wordWithEC = data[cnt];
        cnt++;
    }
    if((data[cnt] & 0xF000) == ID_FOOTPRINT){
        wordWithID = data[cnt];
        cnt++;
    }

    int dataStart = cnt;
    //all cluster words start with 0
    while((data[cnt] & (1<<15)) == 0){
        cnt++;
    }

    //copy payload
    dataSize = cnt - dataStart;
    if(dataSize > 0) {
        dataHolder = new word[dataSize];
        memcpy(dataHolder, &data[dataStart], sizeof(word) * dataSize);
    }

    tailer = data[cnt];
}

ClusterizationVFATFrame::~ClusterizationVFATFrame(){}

bool ClusterizationVFATFrame::channelActive(unsigned char channel) const{
    vector<unsigned char> activeChannels = getActiveChannels();
    return find(activeChannels.begin(), activeChannels.end(), channel) != activeChannels.end();
}

pair<unsigned char, unsigned char> ClusterizationVFATFrame::getClustersAt(unsigned i) const{
    return i < dataSize ? make_pair((dataHolder[i] & 0x7F00)>>8, dataHolder[i] & 0x00FF) : make_pair(0,0);
}

vector<unsigned char> ClusterizationVFATFrame::getActiveChannels() const{
    vector<unsigned char> activeChannels;
    for(unsigned cnt = 0; cnt < dataSize; cnt++){
        pair<unsigned char, unsigned char> cluster = getClustersAt(cnt);
        unsigned short size = cluster.first;
        unsigned short start = cluster.second;
        for(int i = 0; i < size; i++){
            activeChannels.push_back(start + i);
        }
    }
    return activeChannels;
}

vector< pair<unsigned char, unsigned char> > ClusterizationVFATFrame::getActiveClusters() const
{
    vector< pair<unsigned char, unsigned char> > clusters;
    for(unsigned cnt = 0; cnt < dataSize; cnt++){
        clusters.push_back(getClustersAt(cnt));
    }
    return clusters;
}

VFATFrame::word ClusterizationVFATFrame::hasBC() const
{
    return (wordWithBC & 0xF000) == BC_FOOTPRINT;
}

VFATFrame::word ClusterizationVFATFrame::hasEC() const
{
    return (wordWithEC & 0xF000) == EC_FOOTPRINT;
}

VFATFrame::word ClusterizationVFATFrame::hasChipID() const
{
    return (wordWithID & 0xF000) == ID_FOOTPRINT;
}

VFATFrame::word ClusterizationVFATFrame::getBC() const
{
    return hasBC() ? wordWithBC & 0x0FFF : 0;
}

VFATFrame::word ClusterizationVFATFrame::getEC() const
{
    return hasEC() ? (wordWithEC & 0x0FF0) >> 4 : 0;
}

VFATFrame::word ClusterizationVFATFrame::getChipID() const
{
    return hasChipID() ? wordWithID & 0x0FFF : 0;
}

VFATFrame::word ClusterizationVFATFrame::getFlags() const
{
  return hasEC() ? wordWithEC & 0x000F : 0;
}

VFATFrame::word ClusterizationVFATFrame::getSize() const
{
  return tailer & 0x00FF;
}

bool ClusterizationVFATFrame::checkCRC() const
{
  return true;//no crc at all
}

bool ClusterizationVFATFrame::checkFootprint() const
{
  return ((tailer & 0xFF00) == 0xF000)
    && (!hasChipID() || ((wordWithID & 0xF000) == ID_FOOTPRINT))
    && (!hasEC() || ((wordWithEC & 0xF000) == EC_FOOTPRINT))
    && (!hasBC() || ((wordWithBC & 0xF000) == BC_FOOTPRINT))
    && ((header & 0xF000) == 0x8000);
}

VFATFrame::word ClusterizationVFATFrame::getRxFlags() const
{
  return (header & 0x0F00) >> 8;
}

VFATFrame::word ClusterizationVFATFrame::getFiberIdx() const
{
  return (header & 0x00F0) >> 4;
}

VFATFrame::word ClusterizationVFATFrame::getGohIdx() const{
  return header & 0x000F;
}

bool ClusterizationVFATFrame::passedHardwareSynchronisationChecks() const{
  return !(getRxFlags() & 1)                    //if passed crc etc
         && (hasEC() || !(getRxFlags() & 3))    //if BC consistient
         && (hasBC() || !(getRxFlags() & 4));   //if EC consistient
}

void printVfatWord(unsigned short word){
    unsigned short mask = (1 << 15);
    for (int j = 0; j < 16; j++) {
        printf((word & mask) ? "1" : "0");
        mask = (mask >> 1);
        if ((j + 1) % 4 == 0) printf("|");
    }
    printf("\n");
}

void ClusterizationVFATFrame::Print(bool binary) const{
    if (binary) {
        printVfatWord(header);
        if(hasBC()) printVfatWord(wordWithBC);
        if(hasEC()) printVfatWord(wordWithEC);
        if(hasChipID()) printVfatWord(wordWithID);
        for (unsigned i = 0; i < dataSize; i++) {
            printVfatWord(dataHolder[i]);
            printf("\n");
        }
        printVfatWord(tailer);
    } else {
        if(hasChipID()) printf("ID = %03x, ", getChipID());
        if(hasBC()) printf("BC = %04u, ", getBC());
        if(hasEC()) printf("EC = %03u, ", getEC());
        printf("RxFlags = %2u, FiberIdx = %2u, GohIdx = %2u ", getRxFlags(), getFiberIdx(), getGohIdx());
        printf(checkFootprint() ? "footprint   OK" : "footprint FAIL");
        //header
        printf(", frame = %04x|", header);
        if(hasBC()) printf("%04x|", wordWithBC);
        if(hasEC()) printf("%04x|", wordWithEC);
        if(hasChipID()) printf("%04x|", wordWithID);
        //data
        for (unsigned i = 0; i < dataSize; i++) {
            printf("%04x|", dataHolder[i]);
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

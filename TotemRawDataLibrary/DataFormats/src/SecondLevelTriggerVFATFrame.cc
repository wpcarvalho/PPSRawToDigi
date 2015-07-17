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
#include "TotemRawDataLibrary/DataFormats/interface/SecondLevelTriggerVFATFrame.h"

using namespace std;



namespace Totem {

SecondLevelTriggerVFATFrame::SecondLevelTriggerVFATFrame(unsigned short *rawData) {
    unsigned cnt = 0;
    header = rawData[cnt];
    cnt++;
    int dataStart = 1;

    //all words start with 0
    while ((rawData[cnt] & (1 << 15)) == 0) {
        cnt++;
    }

    //copy payload
    dataSize = cnt - dataStart;
    if (dataSize > 0) {
        data = new unsigned short[dataSize];
        memcpy(data, &rawData[dataStart], sizeof(unsigned short) * dataSize);
    }

    tailer = rawData[cnt];
}

SecondLevelTriggerVFATFrame::~SecondLevelTriggerVFATFrame() {
}

unsigned short SecondLevelTriggerVFATFrame::SecondLevelTriggerVFATFrame::getSize() const {
    return tailer & 0x00FF;
}

bool SecondLevelTriggerVFATFrame::checkFootprint() const {
    return ((header & 0xF000) == 0xB000)
        && ((tailer & 0xFF00) == 0xF000);
}

vector< pair<unsigned char, unsigned char> > *SecondLevelTriggerVFATFrame::getData() const
{
    vector< pair<unsigned char, unsigned char> >* result = new vector< pair<unsigned char, unsigned char> >();
    for (unsigned i = 0; i < dataSize; i++) {
        result->push_back(make_pair((data[i] & 0xFF)<<8,data[i] & 0xFF));
    }
    return result;
}

void print2ndLevWordBinary(unsigned short word){
    unsigned short mask = (1 << 15);
    for (unsigned i = 0; i < 16; i++) {
        printf((word & mask) ? "1" : "0");
        mask = (mask >> 1);
        if ((i + 1) % 4 == 0) printf("|");
    }
    printf("\n");
}

void SecondLevelTriggerVFATFrame::Print() const {
    print2ndLevWordBinary(header);
    for (unsigned i = 0; i < dataSize; i++) {
        print2ndLevWordBinary(data[i]);
    }
    print2ndLevWordBinary(tailer);
}

void SecondLevelTriggerVFATFrame::PrintBinary() const {
    printf("header = %04x|", tailer);
    printf("payload = ");
    for (unsigned i = 0; i < dataSize; i++) {
        printf("%04x|", data[i]);
    }
    printf("tailer = %04x|", tailer);
}

}

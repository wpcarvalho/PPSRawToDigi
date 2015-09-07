/**********************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* Authors:
*   Michal Zmuda (m.zmuda@aol.com)
*
**********************************************************/

#pragma once

namespace Totem {

/**
 * TODO: describe.
 *
 * \ingroup TotemRawDataLibrary
 **/
class OptoRxSupplementalData
{
private:
    unsigned long long rawHeader;
    unsigned long long rawTailer;
    unsigned orbitCounter;

public:
    OptoRxSupplementalData(unsigned long long _rawHeader = 0, unsigned long long _rawTailer = 0, unsigned _orbitCounter = 0):
        rawHeader(_rawHeader), rawTailer(_rawTailer), orbitCounter(_orbitCounter) {}
    ~OptoRxSupplementalData() {}

    unsigned short getFOV() const {
        return (rawHeader >> 4) & 0xF;
    }
    unsigned short getOptoRxID() const {
        return (rawHeader >> 8) & 0xFFF;
    }
    unsigned short getBC() const {
        return (rawHeader >> 20) & 0xFFF;
    }
    unsigned getEC() const {
        return (rawHeader >> 32) & 0xFFFFFF;
    }
    unsigned short getEvtTy() const {
        return (rawHeader >> 56) & 0xF;
    }
    unsigned short getBOE() const {
        return (rawHeader >> 60) & 0xF;
    }

    unsigned short getTTCstat() const {
        return (rawTailer >> 4) & 0xF;
    }
    unsigned short getCRC() const {
        return (rawTailer >> 16) & 0xFFFF;
    }
    unsigned short getEventSize() const {
        return (rawTailer >> 32) & 0x3FF;
    }
    unsigned short getEOE() const {
        return (rawTailer >> 60) & 0xF;
    }
    unsigned getOrbitCounter() const {
        return orbitCounter;
    }

    void setHeader(unsigned long long _rawHeader){
        rawHeader = _rawHeader;
    }
    void setTailer(unsigned long long _rawTailer){
        rawTailer = _rawTailer;
    }
};

}


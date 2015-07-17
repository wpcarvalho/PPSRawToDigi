// -*- C++ -*-
//
// Package:    CoincidenceChip
// Class:      CoincidenceChip
//
// Author: Leszek Grzanka

#ifndef _L1TriggerTotemCoincidenceChipCoincidenceChip_H_
#define _L1TriggerTotemCoincidenceChipCoincidenceChip_H_

#include "L1TriggerTotem/CoincidenceChip/interface/CoincidenceChipConfiguration.h"
#include <bitset>

class CoincidenceChip: public CoincidenceChipConfiguration{
  public:
    static const unsigned int OutputBits_size = 16;
    typedef std::bitset<80> InputBits;
    typedef std::bitset<OutputBits_size> OutputBits;
    CoincidenceChip():CoincidenceChipConfiguration(){};
    ~CoincidenceChip(){};
    OutputBits process(InputBits inputBits) const;
    unsigned short NumberOfDetectorPlanes() const;
};


#endif

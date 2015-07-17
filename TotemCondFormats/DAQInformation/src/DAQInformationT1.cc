

#include "FWCore/Utilities/interface/typelookup.h"

#include "TotemCondFormats/DAQInformation/interface/DAQInformationT1.h"

void DAQInformationT1::Reset()
{
  readoutPositionToId.clear();
  readoutIdToRegisters.clear();
/*
  coincidencePositionToId.clear();
  coincidenceIdToRegisters.clear();
*/
}


TYPELOOKUP_DATA_REG(DAQInformationT1);

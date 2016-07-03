
#include "FWCore/Utilities/interface/typelookup.h"

#include "CondFormats/CTPPSReadoutObjects/interface/DiamondAnalysisMask.h"

//----------------------------------------------------------------------------------------------------

void DiamondAnalysisMask::insert(const CTPPSTimingSymbID &sid, const DiamondVFATAnalysisMask &vam)
{
  analysisMask[sid] = vam;
}

//----------------------------------------------------------------------------------------------------

TYPELOOKUP_DATA_REG(DiamondAnalysisMask);

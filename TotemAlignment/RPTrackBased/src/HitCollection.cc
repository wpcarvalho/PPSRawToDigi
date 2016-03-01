/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $Id: StraightTrackAlignment.cc 3363 2010-09-14 19:12:36Z jkaspar $
* $Revision: 2529 $
* $Date: 2010-04-22 17:30:01 +0200 (Thu, 22 Apr 2010) $
*
****************************************************************************/

#include "DataFormats/TotemRPDataTypes/interface/RPRecoHit.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "TotemAlignment/RPTrackBased/interface/HitCollection.h"

#include <cmath>

Hit::Hit(const RPRecoHit &rh) : 
  id(TotRPDetId::RawToDecId(rh.DetId())),
  position(rh.Position()),
  sigma(rh.Sigma())
{
  // prevents null uncertainty, this would make certain algorithms crash
  if (sigma < 1E-6)
    sigma = 66E-3/sqrt(12.);
}


/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $Id: HitCollection.h 9977 2015-01-12 14:00:26Z tsodzawi $
* $Revision: 9977 $
* $Date: 2015-01-12 16:00:26 +0200 (pon, 12 sty 2015) $
*
****************************************************************************/

#ifndef _HitCollection_h_
#define _HitCollection_h_

#include <vector>

class RPRecoHit;

struct Hit {
  unsigned int id;    ///< detector decimal id
  double position;    ///< position in mm
  double sigma;       ///< uncertainty in mm

  Hit(unsigned int i=0, double p=0., double s=0.) : id(i), position(p), sigma(s) {}
  Hit(const RPRecoHit&);
};

typedef std::vector<Hit> HitCollection;

#endif


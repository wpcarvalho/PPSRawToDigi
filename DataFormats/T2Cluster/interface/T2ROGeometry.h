/****************************************************************************
 *
 * This is a part of TOTEM analysis/reconstruction software.
 * Authors:
 *       Giuseppe Latino (Siena Univ. & Pisa INFN), created on 2007/05/10
 *
 ****************************************************************************/
#ifndef T2ROGEOMETRY__
#define T2ROGEOMETRY__

#include "DataFormats/T2DetId/interface/T2DetId.h"
//#include "T2DetId.h"

using namespace std;

class T2ROGeometry
{
 public:
  T2ROGeometry();
  T2ROGeometry(uint32_t det_id);
  virtual ~T2ROGeometry();
  void SetT2DetType(uint32_t det_id);
  inline char GetT2DetType() const {return T2DetType_;}
  double pi;
  //
  // *NOTE* 
  //  when using the following methods:
  //   det_id = T2 detector plane id in cms standard format (32 bit unsigned integer) 
  //   i_raw = strip/pad raw, in range 0-255(0-23) for strip(pads)
  //   j_col = strip/pad column, in range 0-1(0-64) for strip(pads)
  //
  double GetStripPhiMin(int i_raw, int j_col);
  double GetStripPhiMax(int i_raw, int j_col);
  double GetStripRMin(int i_raw, int j_col);
  double GetStripRMax(int i_raw, int j_col);
  double GetPadPhiMin(int i_raw, int j_col);
  double GetPadPhiMax(int i_raw, int j_col);
  double GetPadRMin(int i_raw, int j_col);
  double GetPadRMax(int i_raw, int j_col);
  
  // some methods to retrieve local azimuthal angle in rad
  // starting from 0 close at pad-column 0.
  double GetPadPhiMinLocal(int i_raw, int j_col);
  double GetPadPhiMaxLocal(int i_raw, int j_col);


  unsigned int GetDetId() const {return det_id_;} // Messa dal BERI
  //
 private:  
  uint32_t det_id_;   // detector id
  unsigned int T2arm_, T2halfTelescope_, T2plane_, T2planeSide_;
  char T2DetType_;
  double StripWidth, StripPitch, StripMinR; // strip readout board parameters
  double StripMinPhiDeg, StripMaxPhiDeg, StripCutSchiftY, StripCut;
  double PadRCut, PadMinR, PadRDPhiCut; // pad readout board parameters
  double PadMaxR[24], PadMaxPhiDeg[24], PadDPhiDeg[24];
  double StripPhiMin[256][2];      // Strip min. phi (deg)
  double StripPhiMax[256][2];      // Strip max. phi (deg)
  double StripDPhiMaxMin[256][2];  // Strip max.-min. phi (deg)
  double StripRMin[256][2];        // Strip min. R (mm)
  double StripRMax[256][2];        // Strip max. R (mm)
  double PadRMin[24][65];          // Pad min. R (mm)
  double PadRMax[24][65];          // Pad max. R (mm)
  double PadPhiMin[24][65];        // Pad min. phi (deg)
  double PadPhiMax[24][65];        // Pad max. phi (deg)
};
#endif


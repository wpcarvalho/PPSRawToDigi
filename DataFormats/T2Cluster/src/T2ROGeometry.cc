/****************************************************************************
 *
 * This is a part of TOTEM analysis/reconstruction software.
 * Task: define T2 readout plane (strip/pad) geometry
 *  
 * Authors: Giuseppe Latino (Siena Univ. & Pisa INFN), created on 2007/05/10 
 *
 *
 ****************************************************************************/
//Mirko Berretti: class copied from DataFormat/T2Cluster; The missing methods found 
//in the  Geometry/TotemGeometry/../T2ROGeometry Erik class are added. Geom::pi replaced by pi everywhere


#include "TMath.h"
#include "DataFormats/T2Cluster/interface/T2ROGeometry.h"
//#include "T2ROGeometry.h"
#include "FWCore/Utilities/interface/Exception.h"

//ClassImp(T2ROGeometry)
T2ROGeometry::T2ROGeometry()
{
 pi = 3.1415927;
  //
  //   ** define Strip Parameters **
  //
  StripWidth = 0.080, StripPitch = 0.400;
  StripMinR = 42.460, StripMinPhiDeg = 84.0, StripMaxPhiDeg = 276.0;
  StripCutSchiftY = 1., StripCut = 0.2;
  //
  //   ** define Pad Parameters **
  //
  PadRCut = 0.100, PadMinR = 42.550, PadRDPhiCut = 0.100;
  
  double PadMaxR_vec[24] = {44.677,47.019,49.483,52.076,54.806,57.677,60.700,63.880,
                            67.227,70.726,74.455,78.356,82.461,86.781,91.327,96.111,
			    101.145,106.443,112.011,117.878,124.050,130.547,137.384,144.580};
  double PadMaxPhiDeg_vec[24] = {275.933,275.936,275.939,275.942,275.945,275.948,275.950,275.953,
                                 275.955,275.957,275.960,275.962,275.963,275.965,275.967,275.969,
                                 275.970,275.972,275.973,275.974,275.976,275.977,275.978,275.979};
  double PadDPhiDeg_vec[24] = {2.820,2.826,2.832,2.838,2.844,2.850,2.854,2.860,
                               2.864,2.868,2.873,2.877,2.880,2.884,2.888,2.892,
                               2.894,2.898,2.900,2.902,2.906,2.908,2.910,2.912};
  
  for(int i = 0; i < 24; i++){
    PadMaxR[i] = PadMaxR_vec[i];
    PadMaxPhiDeg[i] = PadMaxPhiDeg_vec[i];
    PadDPhiDeg[i] = PadDPhiDeg_vec[i];
  }
  //
  // *** Define Strip Geometry for design reference T2 readout plane ***
  //
  for(int j = 0; j < 256; j++) {
    for(int k = 0; k < 2; k++) {
      StripRMin[j][k] = StripMinR + j*StripPitch;
      StripRMax[j][k] = StripMinR + StripWidth + j*StripPitch;
    }
    StripPhiMin[j][0] = (acos(StripCutSchiftY/StripRMin[j][0]) + pi/2.)*180./pi;
    StripPhiMax[j][0] = StripMaxPhiDeg;
    StripPhiMin[j][1] = StripMinPhiDeg;
    StripPhiMax[j][1] = (acos((StripCutSchiftY+StripCut)/StripRMin[j][1]) + pi/2.)*180./pi;
    StripDPhiMaxMin[j][0] = StripPhiMax[j][0] - StripPhiMin[j][0];
    StripDPhiMaxMin[j][1] = StripPhiMax[j][1] - StripPhiMin[j][1];
  }
  //
  // *** Define Pad Geometry for design reference T2 readout plane ***
  //
  for(int j = 0; j < 24; j++) {
    for(int k = 0; k < 65; k++) {
      PadRMax[j][k] = PadMaxR[j];
      if(j==0){
	PadRMin[j][k] = PadMinR;
      }else{
	PadRMin[j][k] = PadRMax[j-1][k] + PadRCut;
      }
      if(k==0){
	PadPhiMax[j][k] = PadMaxPhiDeg[j];
      }else{
	PadPhiMax[j][k] = PadPhiMin[j][k-1] - PadRDPhiCut/PadRMin[j][k]*180./pi;
      }
      PadPhiMin[j][k] = PadPhiMax[j][k] - PadDPhiDeg[j];
    }
  }

}

T2ROGeometry::T2ROGeometry(uint32_t det_id)
{
  pi = 3.1415927;
  det_id_ = det_id;
  //
  //   ** define T2 detector numbers (arm, halfTelescope, plane, planeSide) from detector ID number **
  //
  //      det_id_ = T2 detector plane ID (uint32_t) in CMS standard format
  //      Bit 24 = arm (ar): 0==z>0 1==z<0
  //      Bit 23 = halfTelescope (ht): 0==Inner 1==Outer (w.r.t. the beampipe)
  //      Bits [20:22] = plane (pl): [0:4] (starting from IP5)
  //      Bit 19 = planeSide (ps): 0==Front 1==Back (Front = "active area looking at IP5")
  //      - Each T2 detector will be uniquely characterized by these 4 numbers (ar,ht,pl,ps)
  //      - For RO geometry convention purpouses they can be classifies into 4 types:
  //        type "A": (0,1,j,1) (T2-NW) and (1,1,j,0) (T2-SE) j=0->4, as from reference readout board design
  //        type "B": (0,1,j,0) (T2-NW) and (1,1,j,1) (T2-SE) j=0->4, from "A" + 180 deg. rotation around X-axis
  //        type "C": (0,0,j,0) (T2-NW) and (1,0,j,1) (T2-SE) j=0->4, from "A" + 180 deg. rotation around Y-axis
  //        type "D": (0,0,j,1) (T2-NW) and (1,0,j,0) (T2-SE) j=0->4, from "A" + 180 deg. rotation around Z-axis
  //
  T2DetId giveNum;
  T2arm_ = giveNum.arm(det_id_); 
  T2halfTelescope_ = giveNum.halfTelescope(det_id_); 
  T2plane_ = giveNum.plane(det_id_); 
  T2planeSide_ = giveNum.planeSide(det_id_);
  //
  //   ** define T2 detector type **
  //
  T2DetType_ = 'N';
  if( (T2arm_ == 0 && T2halfTelescope_ == 1 && T2plane_ <=4 && T2planeSide_ == 1) ||
      (T2arm_ == 1 && T2halfTelescope_ == 1 && T2plane_ <=4 && T2planeSide_ == 0) )
    T2DetType_ = 'A';

  if( (T2arm_ == 0 && T2halfTelescope_ == 1 && T2plane_ <=4 && T2planeSide_ == 0) ||
      (T2arm_ == 1 && T2halfTelescope_ == 1 && T2plane_ <=4 && T2planeSide_ == 1) )
    T2DetType_ = 'B';

  if( (T2arm_ == 0 && T2halfTelescope_ == 0 && T2plane_ <=4 && T2planeSide_ == 0) ||
      (T2arm_ == 1 && T2halfTelescope_ == 0 && T2plane_ <=4 && T2planeSide_ == 1) )
    T2DetType_ = 'C';

  if( (T2arm_ == 0 && T2halfTelescope_ == 0 && T2plane_ <=4 && T2planeSide_ == 1) ||
      (T2arm_ == 1 && T2halfTelescope_ == 0 && T2plane_ <=4 && T2planeSide_ == 0) )
    T2DetType_ = 'D';
  //
  if(T2DetType_ == 'N') cout<<"T2ROGeometry: warning, wrong T2DetID as input ! => DetID = "<<det_id_<<endl;

  //  cout<<"From T2ROGeometry Sbe: pl:" << T2plane_ <<", pls: "<< T2planeSide_  <<", ht: " <<T2halfTelescope_<<", ID = "<<det_id_<<"\t"<<"Type = "<<T2DetType_<<endl;


  //
  //   ** define Strip Parameters **
  //
  StripWidth = 0.080, StripPitch = 0.400;
  StripMinR = 42.460, StripMinPhiDeg = 84.0, StripMaxPhiDeg = 276.0;
  StripCutSchiftY = 1., StripCut = 0.2;
  //
  //   ** define Pad Parameters **
  //
  PadRCut = 0.100, PadMinR = 42.550, PadRDPhiCut = 0.100;

  double PadMaxR_vec[24] = {44.677,47.019,49.483,52.076,54.806,57.677,60.700,63.880,
                            67.227,70.726,74.455,78.356,82.461,86.781,91.327,96.111,
			    101.145,106.443,112.011,117.878,124.050,130.547,137.384,144.580};
  double PadMaxPhiDeg_vec[24] = {275.933,275.936,275.939,275.942,275.945,275.948,275.950,275.953,
                                 275.955,275.957,275.960,275.962,275.963,275.965,275.967,275.969,
                                 275.970,275.972,275.973,275.974,275.976,275.977,275.978,275.979};
  double PadDPhiDeg_vec[24] = {2.820,2.826,2.832,2.838,2.844,2.850,2.854,2.860,
                               2.864,2.868,2.873,2.877,2.880,2.884,2.888,2.892,
                               2.894,2.898,2.900,2.902,2.906,2.908,2.910,2.912};

  for(int i = 0; i < 24; i++){
    PadMaxR[i] = PadMaxR_vec[i];
    PadMaxPhiDeg[i] = PadMaxPhiDeg_vec[i];
    PadDPhiDeg[i] = PadDPhiDeg_vec[i];
  }
  //
  //
  // *** Define Strip Geometry for design reference T2 readout plane ***
  //
  for(int j = 0; j < 256; j++) {
    for(int k = 0; k < 2; k++) {
      StripRMin[j][k] = StripMinR + j*StripPitch;
      StripRMax[j][k] = StripMinR + StripWidth + j*StripPitch;
    }
    StripPhiMin[j][0] = (acos(StripCutSchiftY/StripRMin[j][0]) + pi/2.)*180./pi;
    StripPhiMax[j][0] = StripMaxPhiDeg;
    StripPhiMin[j][1] = StripMinPhiDeg;
    StripPhiMax[j][1] = (acos((StripCutSchiftY+StripCut)/StripRMin[j][1]) + pi/2.)*180./pi;
    StripDPhiMaxMin[j][0] = StripPhiMax[j][0] - StripPhiMin[j][0];
    StripDPhiMaxMin[j][1] = StripPhiMax[j][1] - StripPhiMin[j][1];
  }
  //
  // *** Define Pad Geometry for design reference T2 readout plane ***
  //
  for(int j = 0; j < 24; j++) {
    for(int k = 0; k < 65; k++) {
      PadRMax[j][k] = PadMaxR[j];
      if(j==0){
	PadRMin[j][k] = PadMinR;
      }else{
	PadRMin[j][k] = PadRMax[j-1][k] + PadRCut;
      }
      if(k==0){
	PadPhiMax[j][k] = PadMaxPhiDeg[j];
      }else{
	PadPhiMax[j][k] = PadPhiMin[j][k-1] - PadRDPhiCut/PadRMin[j][k]*180./pi;
      }
      PadPhiMin[j][k] = PadPhiMax[j][k] - PadDPhiDeg[j];
    }
  }
  //
} // T2ROGeometry Constructor
//
//


void T2ROGeometry::SetT2DetType(uint32_t det_id) {

  det_id_ = det_id;
  //
  //   ** define T2 detector numbers (arm, halfTelescope, plane, planeSide) from detector ID number **
  //
  //      det_id_ = T2 detector plane ID (uint32_t) in CMS standard format
  //      Bit 24 = arm (ar): 0==z>0 1==z<0
  //      Bit 23 = halfTelescope (ht): 0==Inner 1==Outer (w.r.t. the beampipe)
  //      Bits [20:22] = plane (pl): [0:4] (starting from IP5)
  //      Bit 19 = planeSide (ps): 0==Front 1==Back (Front = "active area looking at IP5")
  //      - Each T2 detector will be uniquely characterized by these 4 numbers (ar,ht,pl,ps)
  //      - For RO geometry convention purpouses they can be classifies into 4 types:
  //        type "A": (0,1,j,1) (T2-NW) and (1,1,j,0) (T2-SE) j=0->4, as from reference readout board design
  //        type "B": (0,1,j,0) (T2-NW) and (1,1,j,1) (T2-SE) j=0->4, from "A" + 180 deg. rotation around X-axis
  //        type "C": (0,0,j,0) (T2-NW) and (1,0,j,1) (T2-SE) j=0->4, from "A" + 180 deg. rotation around Y-axis
  //        type "D": (0,0,j,1) (T2-NW) and (1,0,j,0) (T2-SE) j=0->4, from "A" + 180 deg. rotation around Z-axis
  //
  T2DetId giveNum;
  T2arm_ = giveNum.arm(det_id_); 
  T2halfTelescope_ = giveNum.halfTelescope(det_id_); 
  T2plane_ = giveNum.plane(det_id_); 
  T2planeSide_ = giveNum.planeSide(det_id_);
  //
  //   ** define T2 detector type **
  //
  T2DetType_ = 'N';
  if( (T2arm_ == 0 && T2halfTelescope_ == 1 && T2plane_ <=4 && T2planeSide_ == 1) ||
      (T2arm_ == 1 && T2halfTelescope_ == 1 && T2plane_ <=4 && T2planeSide_ == 0) )
    T2DetType_ = 'A';
  
  if( (T2arm_ == 0 && T2halfTelescope_ == 1 && T2plane_ <=4 && T2planeSide_ == 0) ||
      (T2arm_ == 1 && T2halfTelescope_ == 1 && T2plane_ <=4 && T2planeSide_ == 1) )
    T2DetType_ = 'B';
  
  if( (T2arm_ == 0 && T2halfTelescope_ == 0 && T2plane_ <=4 && T2planeSide_ == 0) ||
      (T2arm_ == 1 && T2halfTelescope_ == 0 && T2plane_ <=4 && T2planeSide_ == 1) )
    T2DetType_ = 'C';
  
  if( (T2arm_ == 0 && T2halfTelescope_ == 0 && T2plane_ <=4 && T2planeSide_ == 1) ||
      (T2arm_ == 1 && T2halfTelescope_ == 0 && T2plane_ <=4 && T2planeSide_ == 0) )
    T2DetType_ = 'D';
  //
  if(T2DetType_ == 'N') cout<<"T2ROGeometry: warning, wrong T2DetID as input ! => DetID = "<<det_id_<<endl;
  
  //  cout<<"From T2ROGeometry Er: pl:" << T2plane_ <<", pls: "<< T2planeSide_  <<", ht: " 
  //     <<T2halfTelescope_<<", ID = "<<det_id_<<"\t"<<"Type = "<<T2DetType_<<endl;

} // SetT2DetType



T2ROGeometry::~T2ROGeometry()
{
} // T2ROGeometry Destructor


//
//
// Methods
// *NOTE* =>  when using the following methods you need 2 inputs:
//
// 1)  i_raw = strip/pad raw, in range 0-255(0-23) for strip(pads)
// 2)  j_col = strip/pad column, in range 0-1(0-64) for strip(pads)
//

double T2ROGeometry::GetStripPhiMin(int i_raw, int j_col)
{
  if(T2DetType_ == 'A') {
    return StripPhiMin[i_raw][j_col];
  }
  //
  if(T2DetType_ == 'B') {
    if(j_col == 0) return StripMinPhiDeg;
    if(j_col == 1) return (StripMaxPhiDeg - StripDPhiMaxMin[i_raw][1]); 
  }
  //
  if(T2DetType_ == 'C') {
    if(j_col == 0) return (StripMinPhiDeg + 180.);
    if(j_col == 1) return (StripMaxPhiDeg -180. - StripDPhiMaxMin[i_raw][1]);
  }
  //
  if(T2DetType_ == 'D') { 
    if(StripPhiMin[i_raw][j_col] + 180. <= 360.) {
      return StripPhiMin[i_raw][j_col] + 180.;
    }else{
      return StripPhiMin[i_raw][j_col] - 180.;
    }
  }
  return 0.0;
  //
}
//
//
double T2ROGeometry::GetStripPhiMax(int i_raw, int j_col)
{
  if(T2DetType_ == 'A') {
    return StripPhiMax[i_raw][j_col];
  }
  //
  if(T2DetType_ == 'B') {
    if(j_col == 0) return (StripMinPhiDeg + StripDPhiMaxMin[i_raw][0]);
    if(j_col == 1) return StripMaxPhiDeg;
  }
  //
  if(T2DetType_ == 'C') {
    if(j_col == 0)  return (StripMinPhiDeg + StripDPhiMaxMin[i_raw][0] - 180.);
    if(j_col == 1) return (StripMaxPhiDeg - 180.);
  }
  //
  if(T2DetType_ == 'D') {
    if(StripPhiMax[i_raw][j_col] + 180. <= 360.) {
      return StripPhiMax[i_raw][j_col] + 180.;
    }else{
      return StripPhiMax[i_raw][j_col] - 180.;
    }
  }
  throw cms::Exception("ERROR") << "This point should not be reached!";
  //
}
//
//
double T2ROGeometry::GetStripRMin(int i_raw, int j_col)
{
  return StripRMin[i_raw][j_col]; // No R dependence on detector ID
}
//
//
double T2ROGeometry::GetStripRMax(int i_raw, int j_col)
{
  return StripRMax[i_raw][j_col]; // No R dependence on detector ID
}
//
//
double T2ROGeometry::GetPadPhiMin(int i_raw, int j_col)
{
  //cout<<" dentro GetPadPhiMin"<<endl;
  if(T2DetType_ == 'A') {
    // cout<<"GetPadPhiMin typeA: return  PadPhiMin["<<i_raw<<"]["<<j_col<<"] ="<<  PadPhiMin[i_raw][j_col]<<endl;
    return PadPhiMin[i_raw][j_col];
  }
  //
  if(T2DetType_ == 'B') {
    // cout<<"GetPadPhiMin typeB: return  PadPhiMin["<<i_raw<<"]["<<64 - j_col<<"] ="<<  PadPhiMin[i_raw][64 - j_col]<<endl;
    return PadPhiMin[i_raw][64 - j_col];
  
  }
  //
  if(T2DetType_ == 'C') {

    // if(j_col <= 32)
    //cout<<"GetPadPhiMin typeC: jcol<=32("<<j_col<<"): "<<(3.*180. - PadPhiMax[i_raw][j_col])<<endl; 
    //if(j_col > 32)
    //cout<<"GetPadPhiMin typeC: jcol>32("<<j_col<<"): "<<(180. - PadPhiMax[i_raw][j_col])<<endl; 

    if(j_col <= 32) return (3.*180. - PadPhiMax[i_raw][j_col]);

    if(j_col > 32) return (180. - PadPhiMax[i_raw][j_col]);
  }
  //
  if(T2DetType_ == 'D') {
    // cout<<"GetPadPhiMin typeD:"<<endl; 
    if(PadPhiMin[i_raw][j_col] + 180. <= 360.) {      
      return PadPhiMin[i_raw][j_col] + 180.;
    }else{
      
      return PadPhiMin[i_raw][j_col] - 180.;
    }
  }
 throw cms::Exception("ERROR") << "This point should not be reached!";
}
//
//
double T2ROGeometry::GetPadPhiMax(int i_raw, int j_col)
{
  if(T2DetType_ == 'A') {
    return PadPhiMax[i_raw][j_col];
  }
  //
  if(T2DetType_ == 'B') {
    return PadPhiMax[i_raw][64 - j_col];
  }
  //
  if(T2DetType_ == 'C') {
    //  if(j_col < 32)
      //cout<<"GetPadPhiMax typeC: jcol<32("<<j_col<<"): "<<(3.*180. - PadPhiMin[i_raw][j_col])<<endl; 
    //if(j_col >= 32)
    //cout<<"GetPadPhiMax typeC: jcol>=32("<<j_col<<"): "<<(180. - PadPhiMin[i_raw][j_col])<<endl; 
    
    if(j_col < 32) return (3.*180. - PadPhiMin[i_raw][j_col]);
    if(j_col >= 32) return (180. - PadPhiMin[i_raw][j_col]);
  }
  //
  if(T2DetType_ == 'D') {
    if(PadPhiMax[i_raw][j_col] + 180. <= 360.) {
      return PadPhiMax[i_raw][j_col] + 180.;
    }else{
      return PadPhiMax[i_raw][j_col] - 180.;
    }
  }
  throw cms::Exception("ERROR") << "This point should not be reached!";
}
//
//
double T2ROGeometry::GetPadRMin(int i_raw, int j_col)
{
  return PadRMin[i_raw][j_col]; // No R dependence on detector ID
}
//
//
double T2ROGeometry::GetPadRMax(int i_raw, int j_col)
{
  return PadRMax[i_raw][j_col]; // No R dependence on detector ID
}


double T2ROGeometry::GetPadPhiMinLocal(int i_raw, int j_col) {
  
  return (PadPhiMin[i_raw][j_col] - 180) /180 *pi /*Geom::pi()*/;
  
} // GetPadPhiMinLocal

/**
 *
 */

double T2ROGeometry::GetPadPhiMaxLocal(int i_raw, int j_col) {
  
  return (PadPhiMax[i_raw][j_col] - 180) / 180 * pi/*Geom::pi()*/;
  
} // GetPadPhiMaxLocal

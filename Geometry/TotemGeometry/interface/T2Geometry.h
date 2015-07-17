#ifndef TotemGeometry_T2Geometry_h
#define TotemGeometry_T2Geometry_h

/** 
 * Class T2Digitizer for digitization of T2 GEM detector
 * Here CMSConventions are used especially their 
 * coordinate system
 * 
 * Author: Erik Brücken / University of Helsinki
 * Email: brucken@cc.helsinki.fi
 *
 */

#include "DataFormats/T2DetId/interface/T2DetId.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"
#include "DataFormats/T2Cluster/interface/T2ROGeometry.h"
#include <iostream>

//class T2Pad;
//class T2Strip; 

class T2Geometry {

 private:

  int sRows, sColumns, pRows, pColumns; // total number of rows and columns
  int planeNr;
  /*unsigned int*/uint32_t unitId;
  double rMinPlane, rMaxPlane; 
  double phiInitPlane, phiTotPlane;
  double zPosPlane;
  double cAlpha; // factor to get dimension of pad-rows

  //  T2Strip* strip;


 public:

  T2Geometry();
  ~T2Geometry(){};

  void setPlane(/*unsigned int*/uint32_t unitId); // set the correct dimensions
  int getNearestStripRow(Local3DPoint* hitPos);
  int getNearestStripCol(Local3DPoint* hitPos);
  int getNearestPadRow(Local3DPoint* hitPos);
  int getNearestPadCol(Local3DPoint* hitPos);
  int getNearestPadRow(double radius, double phi);
  // phi in radians
  int getNearestPadCol(double radius, double phi);
  inline /*unsigned int*/uint32_t getUnitId() const { return unitId; };

  
  //Mirko Berretti

  int getNearestStripRow0_256(double radius,unsigned int theid); 
  //Strip row from 0-127
  int getNearestStripRow_(double radius,unsigned int theid);
  //Strip col from 0-1; phi in radians
  int getNearestStripCol_(double phi,unsigned int theid);
  int getNearestPadRow_(double radius, double phi);
  // phi in radians
  int getNearestPadCol_(double radius,double phi,unsigned int theid);
  int getNearestPadRow2_(Local3DPoint* hitPos,unsigned int theid);
};

#endif // TotemGeometry_T2Geometry_h
/*
int T2Geometry::getNearestPadCol_(double Y,double phi,unsigned int theid)
int T2Geometry::getNearestStripRow_(double radius) 
int T2Geometry::getNearestPadRow_(double radius, double phi) 
int T2Geometry::getNearestStripCol_(double Y,double phi,unsigned int theid) getNearestPadRow2_(Local3DPoint* hitPos,unsigned int theid)
*/

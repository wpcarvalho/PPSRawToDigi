#ifndef TotemGeometry_T2GeometryUtil_h
#define TotemGeometry_T2GeometryUtil_h

/** 
 * Class T2GeometryUtil. Utilities for T2 detector position / identification
 * Here CMSConventions are used especially their 
 * coordinate system
 * 
 * Author: Mirko Berretti
 * Email: mirko.berretti@gmail.com
 *
 */

#include "DataFormats/T2DetId/interface/T2DetId.h"
#include "Geometry/TotemGeometry/interface/T2Geometry.h"
#include "DataFormats/T2Hit/interface/T2Hit.h"


class T2GeometryUtil {

 private:

  T2DetId t2detid;

 

 public:

  T2GeometryUtil();
  virtual ~T2GeometryUtil();

    struct T2DetInfo{
      unsigned int arm;
      unsigned int ht; 
      unsigned int pl; 
      unsigned int plside;  
      unsigned int pl_0to9;
      unsigned int symb;
      uint32_t cmsswid;
      double Zdet;
    };

    typedef struct{
      int vfatiid;
      unsigned int channel;
  } vfatid_channel;

    //T2DetInfo thet2plane;

    std::string DetInfoFromRawId(uint32_t);

    int HV_Ysign(uint32_t thedetID);
    
    //if vector size=0, means there are convertion problems.
    vfatid_channel PadVfatsIdFromRecoHit(T2Hit hit);
    std::vector<vfatid_channel> PadVfatsIdsFromPadVect(T2Hit hit);  
    vfatid_channel PadVfatsIdFromRowCol(int row, int col, uint32_t thedet); 
    

    vfatid_channel StripVfatsIdFromRecoHit(T2Hit hit);   
    std::vector<vfatid_channel> StripVfatsIdsFromStripVect(T2Hit hit);   
    vfatid_channel StripVfatsIdFromRowCol(int row, int col, uint32_t thedet);   

    static unsigned int RawtoSymb(uint32_t thedet);

    static double Zposition(unsigned int arm, unsigned int ht, unsigned int pl, unsigned int plside);

    static double Zposition(uint32_t cmsswid);
  
    static T2DetInfo GetT2Info(uint32_t  cmsswid);
    //T2DetInfo GetT2Info(unsigned int symbid);
    // unsigned int GetSymbId(uint32_ cmsswid);


};

#endif





/****************************************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* Authors: 
*	Erik Brücken, University of Helsinki, email: brucken@cc.helsinki.fi
*	Jan Kaspar (jan.kaspar@gmail.com) 
*
* Copied by Fredrik Oljemark 2009/08/20
****************************************************************************/

/*   
 *   Program to convert vfatIID and channel to geometrical strips and pads 
 *   numbers (Giuseppe's layout).
 *   Hardware layout and geometrical layout are shown in the documentation 
 *   that can be found in the same directory.
 *
 *   Author: Erik Brücken, University of Helsinki
 *    email: brucken@cc.helsinki.fi
 */


//include "Rtypes.h"
//include <vector>

//namespace Totem {
//	class VFATFrameColl;
//}

using namespace std;

class T2TrigMap
{
 public:

  int getMetaSector(int TrigVFat,int VFatChannel) {
    if ((TrigVFat==0||TrigVFat==1)&&(VFatChannel%2==0)&&(VFatChannel<128)) {
      iMetaS=VFatChannel/16 + 8*TrigVFat;
      return iMetaS;
    }
    else {return -1;}
  }
  
  int getMetaTrigSector(int TrigVFat,int VFatChannel) {
    if ((TrigVFat==0||TrigVFat==1)&&(VFatChannel%2==0)&&(VFatChannel<128)) {
      iMetaTS=(VFatChannel/2)%8;
      return iMetaTS;
    }
    else {return -1;}
  }
  
  void setMap(const std::map<int,int>& theMap) {
    tmpChaVsMetapad=theMap;
  }
 
 private:
  int iMetaTS,iMetaS,iMapLoaded;
  std::map<int,int> tmpChaVsMetapad;
		
};


/*

class elementContainer {
  
private:
  
  Int_t elementFlavor;
  Int_t elementNR;
  
public:
  
  elementContainer(Bool_t elementFlavor, Int_t elementNR);
  
  inline Int_t getElementFlavor() { return elementFlavor; };
  inline Int_t getElementNR() { return elementNR; };
  
}; // elementContainer

elementContainer::elementContainer(Bool_t elementFlavor, Int_t elementNR) {

  this->elementFlavor = elementFlavor;
  this->elementNR = elementNR;

} // elementContainer::elementContainer

*/

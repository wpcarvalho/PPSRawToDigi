/****************************************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* Authors: 
*	Erik Brücken, University of Helsinki, email: brucken@cc.helsinki.fi
*	Jan Kaspar (jan.kaspar@gmail.com) 
*    
* $Id: T2Mapping.h 2135 2010-02-16 13:26:06Z jkaspar $
* $Revision: 2135 $
* $Date: 2010-02-16 14:26:06 +0100 (Tue, 16 Feb 2010) $
*
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


#include "Rtypes.h"
#include <vector>

namespace Totem {
	class VFATFrameCollection;
}

using namespace std;

class T2Mapping
{
	public:
		
		/// Flavour flags
		enum {fStrip, fPad};

		/// Function to convert vfatIID and channelnumber to geometrical strip and padnumber
		/// The IID is Internal ID, i.e. it is NOT that 16 or 12 bit number, it is describing possition
		/// of VFAT within one plane; the range is 0..16
		Bool_t convertToGeo(Int_t vfatIID, Int_t channelNR, Int_t &elementFlag, Int_t &elementNR);

		/// An interface between VFAT frame collection and cluster reconsctruction
		/// input: collection of VFAT frames 'frames' and vector of IDs indexed by IID 0..16
		/// output: vectors of active strips and pads
		/// returns number of conversion errors
		int convertToStripsPads(Totem::VFATFrameCollection *frames, vector<int> IDs, vector<int>& strips, vector<int>& pads);
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

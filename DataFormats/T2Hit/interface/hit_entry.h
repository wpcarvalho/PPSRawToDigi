#ifndef hit_entry_h
#define hit_entry_h

#include <vector>
//#include "Rtypes.h"

class hit_entry
{
 public:
 
  hit_entry() 
    {rCl_=0; phiCl_=0; dphiCl_=0; drCl_=0; zCl_=0; dzCl_=0;}

  hit_entry(float rCl, float phiCl, float drCl, float dphiCl, float zCl, float dzCl, 
	    bool isPad, unsigned int numentry) 
    {rCl_=rCl; phiCl_=phiCl; dphiCl_=dphiCl; drCl_=drCl; zCl_=zCl; 
    dzCl_=dzCl;  isPad_=isPad; nument_=numentry;}
    
  float rCl_; 
  float phiCl_;
  float dphiCl_;
  float drCl_;
  float zCl_;
  float dzCl_;	 
  bool isPad_;
  unsigned int nument_;	
   
};  

#endif

/*
  Class T2Road
  Author: Mirko Berretti
  mirko.berretti@gmail.com
*/
#include "DataFormats/T2Road/interface/T2Road.h"

void T2Road::CalculateRoadExtreme()
{
  double rmin=0.;
  double rmax=0.;
  double phimin=0.;
  double phimax=0.;
  std::vector<T2Hit> thisRd=this->thisRoad;
  if(thisRd.size()==0)
    cout<<"Empty road, impoossible to calculate the extreme"<<endl;
  else
    {
      rmin=thisRd[0].GetHitR();
      rmax=thisRd[0].GetHitR();
      phimin=thisRd[0].GetHitPhi();
      phimax=thisRd[0].GetHitPhi();
      for(unsigned int k=0; k<thisRd.size();k++ )
	{
	  if(thisRd[k].GetHitR()>rmax)
	    rmax=thisRd[k].GetHitR();
	  if(thisRd[k].GetHitR()<rmin)
	    rmin=thisRd[k].GetHitR();
	  if(thisRd[k].GetHitPhi()<phimin)
	    phimin=thisRd[k].GetHitPhi();
	  if(thisRd[k].GetHitPhi()>phimax)
	    phimax=thisRd[k].GetHitPhi();	       
	}
      this->SetRoadPhimin(phimin);
      this->SetRoadPhimax(phimax);
      this->SetRoadRmin(rmin);
      this->SetRoadRmax(rmax);	 

    }
  
}

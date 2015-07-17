/*
  Class T2Hit
  Author: Mirko Berretti
  mirko.berretti@gmail.com
*/
/*
  Author: Mirko Berretti
  mirko.berretti@gmail.com
  University of Siena
*/
#include "DataFormats/T2Hit/interface/T2Hit.h"
#include <TMath.h>



T2Hit::T2Hit()
{}

T2Hit::~T2Hit()
{}

//Used for vtx inclusions.


T2Hit::T2Hit(double X, double Y, double Z, double EX, double EY, double EZ, uint32_t cmsswid)
{
  double R=sqrt(X*X+Y*Y); 

  double Phi=0;

  if(X==0)
    Phi=0;
  else
    Phi=atan(fabs(Y)/fabs(X));

  //std::cout<<"A"<<std::endl;
  Phi=Phi*180.0/3.14159;	
			    
  if((Y<0)&&(X>0))
    Phi=360.0-Phi;
      
  if((Y>0)&&(X<0))
    Phi=180.0-Phi;

      
  if((Y<0)&&(X<0))
    Phi=Phi+180.;

  hitR_=R;
  hitPhi_=Phi;
  hitZ_=Z;
  hitX_=X;
  hitY_=Y;

  hitDX_=EX;
  hitDY_=EY;
 
  //Pag 70 Statistics for nuclear and particle physics.
  if(R!=0)
    {
      hitDR_=sqrt((1.0/R)*(1.0/R)*(X*X*EX*EX+Y*Y*EY*EY));
      hitDPhi_=sqrt((1.0/R)*(1.0/R)*(1.0/R)*(1.0/R)*(Y*Y*EX*EX+X*X*EY*EY))*180/3.14159;
    }
  else
    {
      hitDR_= std::max(hitDX_,hitDY_);					
      hitDPhi_=1.5;
    }

  hitDZ_=EZ;
  VtxHits=true;
  HitInTrk=true;
  hitClass_=9;//Only used for this purpose
  detrawid_=cmsswid;
  
}


T2Hit::T2Hit(double expX_2, double expY_2, double Z_, uint32_t cmsswid)
{

  double expR_2=sqrt(expX_2*expX_2+expY_2*expY_2);  
  
  double expPhi_2=atan(fabs(expY_2)/fabs(expX_2));
  // std::cout<<"B"<<std::endl;
  expPhi_2=expPhi_2*180.0/3.14159;	
			    
  if((expY_2<0)&&(expX_2>0))
    expPhi_2=360.0-expPhi_2;
      
  if((expY_2>0)&&(expX_2<0))
    expPhi_2=180.0-expPhi_2;

      
  if((expY_2<0)&&(expX_2<0))
    expPhi_2=expPhi_2+180.;

  hitR_=expR_2;
  hitPhi_=expPhi_2;
  hitZ_=Z_;
  hitX_=expX_2;
  hitY_=expY_2;
  detrawid_=cmsswid;
 
}

bool operator<(const T2Hit &h1,const T2Hit &h2)
{
return h1.hitZ_ < h2.hitZ_;
}
bool operator>(const T2Hit &h1,const T2Hit &h2)
{
return h1.hitZ_ > h2.hitZ_;
}
bool operator==(const T2Hit &h1,const T2Hit &h2)
{
return h1.hitZ_ == h2.hitZ_;
}
bool operator!=(const T2Hit &h1,const T2Hit &h2)
{
return h1.hitZ_ != h2.hitZ_;
}


void T2Hit::ComputeHit()
{

  float BestdR, BestR, BestdPhi, BestPhi;
  bool Allcheck=true;
  unsigned int cltype;
  unsigned int numpad=0;
  unsigned int numstrip=0;


  if((thehitentries.size()==0)||(thehitentries.size()>2))
    {
      if(thehitentries.size()==0)
	cout<<" No Cluster present to make a Hit"<<endl;
      else
	cout<<" ()()()() Cannot produce a hit with three clusters ()()()() "<<endl;
      Allcheck=false;
     
      cout << "Warning T2Hit: number of matched cluster to find hit = "<<thehitentries.size();
    }
  else
    {
      BestdR= thehitentries[0].drCl_;
      BestR= thehitentries[0].rCl_;
      BestdPhi= thehitentries[0].dphiCl_;
      BestPhi=thehitentries[0].phiCl_;
      unsigned int i=0;
      while((i<thehitentries.size())&&(Allcheck))
	{
	  if (thehitentries[i].zCl_ != thehitentries[0].zCl_)
	    {
	      Allcheck=false;
	      cout<<" Cannot produce a hit with clusters at different z "<<endl;
	    }  
	  if(Allcheck)
	    {
              if(thehitentries[i].isPad_)
		numpad=numpad+thehitentries[i].nument_;
	      else
		numstrip=numstrip+thehitentries[i].nument_;

    
	      if(thehitentries[i].drCl_< BestdR)
		{
		  BestdR=thehitentries[i].drCl_;
		  BestR=thehitentries[i].rCl_;
		}
	      if(thehitentries[i].dphiCl_< BestdPhi)
		{
		  BestdPhi=thehitentries[i].dphiCl_;
		  BestPhi=thehitentries[i].phiCl_;
		}
	    }
	  i++;
	}
      
     

      if(Allcheck)
	{
	  /*
	  bool aflag=false;
	  this->SetVtxHit(aflag); 
	  this->SetInTrk(aflag);
	  */

	  this->SetHitNumStrip(numstrip);
	  this->SetHitNumPad(numpad);
	   
	  if((numpad>0)&&(numstrip>0)) 
	    cltype=1;
	  else
	    cltype=2;
	  
	  if((numpad==1)&&(cltype==1))   //has to be changed after T2ROGeometry modification
	    {
	      BestdPhi=1.5;
	    }
	  if((numpad>1)&&(cltype==1))   //has to be changed after T2ROGeometry modification
	    {
	      BestdPhi=3.0;
	    }
	   
      
	  
	  this->SetHitClass(cltype);
	  this->SetHitR(BestR);
	  this->SetHitDR(BestdR);
	  this->SetHitPhi(BestPhi);
	  this->SetHitDPhi(BestdPhi);
	  this->SetHitZ(thehitentries[0].zCl_);
	  this->SetHitDZ(thehitentries[0].dzCl_);
	  
	 
	  double thephirad=BestPhi*3.14159265/180.0;
	  double theDphirad=BestdPhi*3.14159265/180.0;
	 

	  this->SetHitX(BestR*cos(thephirad));
	  this->SetHitY(BestR*sin(thephirad));

	  //The error on X,Y are from propagation formula
	  
	  
	  this->SetHitDX(sqrt(BestdR*BestdR*cos(thephirad)*cos(thephirad)+theDphirad*theDphirad*BestR*BestR*sin(thephirad)*sin(thephirad)));
          
	  
	  this->SetHitDY(sqrt(BestdR*BestdR*sin(thephirad)*sin(thephirad)+theDphirad*theDphirad*BestR*BestR*cos(thephirad)*cos(thephirad)));
	  

	}  
    }


}




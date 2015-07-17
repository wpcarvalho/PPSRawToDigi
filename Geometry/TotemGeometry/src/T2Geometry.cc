/** 
 * Class T2Digitizer for digitization of T2 GEM detector
 * Here CMSConventions are used especially their 
 * coordinate system
 * 
 * Author: Erik Brücken / University of Helsinki
 * Email: brucken@cc.helsinki.fi
 * Mirko Berretti: University of Siena and Pisa INFN
 * Email: mirko.berretti@gmail.com
 *
 * At the moment the dimensions are in "mm" not in "cm".
 * That will change in the furure to be compilant to the 
 * CMS conventions.
 */

#include "Geometry/TotemGeometry/interface/T2Geometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iostream>
#include <cmath>


/**
 * Constructor to initalize geometry related dimensions 
 */

T2Geometry::T2Geometry() {
  
  rMinPlane =  42.46; // all dimensions should be in cm (CMS convention)
  // for pads 42.55 and for strips 42.46 TO CHECK!!!!
  rMaxPlane = 144.58; // for pads 144.58 and for strips 144.54 TO CHECK!!!!
  phiTotPlane = 3.351032165; // rad (192°)
  phiInitPlane = 0; // default
  zPosPlane = 0; // default
  planeNr = -1;
  unitId = 0;
  sRows = 256;
  sColumns = 2;
  pRows = 24;
  pColumns = 65;
  cAlpha = log(pow(rMaxPlane/rMinPlane, 1.0/pRows));

}

/**
 * Set plane using the detector unitId. 
 */

void T2Geometry::setPlane(/*unsigned int*/uint32_t unitId) {

  // std::cout<<"T2GeometryA called with"<<unitId<<std::endl;
  this->unitId = unitId;
  
  int arm = T2DetId::arm(unitId);
  int halfTelescope = T2DetId::halfTelescope(unitId);
  int plane = T2DetId::plane(unitId);
  int planeSide = T2DetId::planeSide(unitId);
  LogDebug("T2Geometry") << "Arm " << arm << ", halfTelescope " << halfTelescope 
			 << ", plane " << plane << ", planeSide " << planeSide;

  //  phiInitPlane = (halfTelescope == 0) ? -1.675516082 : 1.466076572;
  phiInitPlane = (halfTelescope == 0) ? -1.675516082 : -1.675516082;
  //std::cout<<"T2GeometryC"<< "Arm " << arm << ", halfTelescope " << halfTelescope << ", plane " << plane << ", planeSide " << planeSide<<" AbsId:"<<this->getUnitId()<<std::endl;
}
 
/**
 *
 */

int T2Geometry::getNearestStripRow(Local3DPoint* hitPos) {
  
  double x = hitPos->x();
  double y = hitPos->y();
  double r = sqrt(x*x + y*y);
  
  // Calculate row of strips:
  int sRow = (int) floor((r-rMinPlane) / (rMaxPlane-rMinPlane-0.04) * sRows);
  
  // Small corrections if we are on the "border" or hit by a rounding mistake...
  if (sRow>=sRows)      sRow = sRows-1;
  if (sRow<0)           sRow = 0;

  LogDebug("T2Geometry") << "Nearest strip radius: " << r << ", sRow: " << sRow;    

  return sRow;


} // getNearestStripRow

/**
 *
 */

int T2Geometry::getNearestStripCol(Local3DPoint* hitPos) {
  
  double x = hitPos->x();
  double y = hitPos->y();
  double phi = (x==0 && y==0) ? 0 : atan2(y,x); // -pi < phi <= +pi
  double phiRel;

  /*
    if(phi < 0 && phiInitPlane > 0) {
    
    phiRel = phi - phiInitPlane + 2 * M_PI; 

    } else phiRel = phi - phiInitPlane;
  */

  phiRel = phi - phiInitPlane;
  // Calculate col of strips
  int sCol = (int) floor(phiRel / phiTotPlane * sColumns);
 
  // Small corrections if we are on the "border" or hit by a rounding mistake...
  if (sCol>=sColumns)   sCol = sColumns-1;
  if (sCol<0)           sCol = 0;
  
  LogDebug("T2Geometry") << "Strip: phi " << phi <<", phiRel "<< phiRel << ", phiinit " << phiInitPlane 
			 << ", phitot " << phiTotPlane <<", Column "<< (sColumns-1) - sCol;
  
  return (sColumns-1) - sCol;

} // getNearestStripCol

/**
 *
 */

int T2Geometry::getNearestPadRow(Local3DPoint* hitPos) {
  
  double x = hitPos->x();
  double y = hitPos->y();
  double r = sqrt(x*x + y*y);
  
  // Calculate row of pads
  int pRow = (int) floor(log(r/(rMinPlane+0.09))/cAlpha);
  
  // Small corrections if we are on the "border" or hit by a rounding mistake...
  if (pRow>=pRows)      pRow = pRows-1;
  if (pRow<0)           pRow = 0;

  LogDebug("T2Geometry") << "Nearest pad radius: " << r << ", pRow: " << pRow;    

  return pRow;

} // getNearestPadRow

/**
 *
 */

int T2Geometry::getNearestPadCol(Local3DPoint* hitPos) {
  
  double x = hitPos->x();
  double y = hitPos->y();
  double phi = (x==0 && y==0) ? 0 : atan2(y,x); // -pi < phi <= +pi
  double phiRel;
 
  /* 
     if(phi < 0 && phiInitPlane > 0) {
    
     phiRel = phi - phiInitPlane + 2 * M_PI; 

     } else phiRel = phi - phiInitPlane;
  */

  phiRel = phi - phiInitPlane;

  // Calculate col of pads 
  int pCol = (int) floor(phiRel / phiTotPlane * pColumns);
  
  // Small corrections if we are on the "border" or hit by a rounding mistake...
  if (pCol>=pColumns)   pCol = pColumns-1;
  if (pCol<0)           pCol = 0;

  // std::cout<<"Pad: phi " << phi <<", phiRel "<< phiRel << ", phiinit " << phiInitPlane 
  //   << ", phitot " << phiTotPlane << ", Column "<< (pColumns-1) - pCol<<std::endl;

  LogDebug("T2Geometry") << "Pad: phi " << phi <<", phiRel "<< phiRel << ", phiinit " << phiInitPlane 
			 << ", phitot " << phiTotPlane << ", Column "<< (pColumns-1) - pCol;

  return (pColumns-1) - pCol;

} // getNearestPadCol

/**
 *
 */

int T2Geometry::getNearestPadRow(double radius, double phi) {
  
  // Calculate row of pads
  int pRow = (int) floor(log(radius/(rMinPlane+0.09))/cAlpha);
  
  // Small corrections if we are on the "border" or hit by a rounding mistake...
  if (pRow>=pRows)      pRow = pRows-1;
  if (pRow<0)           pRow = 0;

  LogDebug("T2Geometry") << "Nearest pad radius: " << radius << ", pRow: " << pRow;    

  return pRow;

} // getNearestPadRow

/**
 *
 */


int T2Geometry::getNearestPadCol(double radius, double phi) {
  
  double phiRel;
 
  /* 
     if(phi < 0 && phiInitPlane > 0) {
    
     phiRel = phi - phiInitPlane + 2 * M_PI; 

     } else phiRel = phi - phiInitPlane;
  */

  phiRel = phi - phiInitPlane;

  // Calculate col of pads 
  int pCol = (int) floor(phiRel / phiTotPlane * pColumns);
  
  // Small corrections if we are on the "border" or hit by a rounding mistake...
  if (pCol>=pColumns)   pCol = pColumns-1;
  if (pCol<0)           pCol = 0;

  LogDebug("T2Geometry") << "Pad: phi " << phi <<", phiRel "<< phiRel << ", phiinit " << phiInitPlane 
			 << ", phitot " << phiTotPlane << ", Column "<< (pColumns-1) - pCol;

  return (pColumns-1) - pCol;

} // getNearestPadCol







//Added in 29 March 2010 - Mirko Berretti


//Strip row from 0-127
int T2Geometry::getNearestStripRow_(double radius,unsigned int theid) {
  
  double StriprMin=  42.46; 
  double StriprMax = 144.54; 
  int striprow=-1;

  if(((radius-StriprMin)<=0)||((radius-StriprMax)>=0))
    {
       if((radius-StriprMin)<=0)
	 striprow=0;

       if((radius-StriprMax)>=0)
	 striprow=255;
    }
  else
    {
      //0.04= half strip width;
      //0.4= strip pitch
      striprow= (radius-(StriprMin-0.04));//(StriprMin+0.04-(StriprMax-0.04))/0.4;
    }
  
  int striprowInVfat=striprow%128;
  

  return striprowInVfat;
}



int T2Geometry::getNearestStripRow0_256(double radius,unsigned int theid) {
  
  double StriprMin=  42.46; 
  double StriprMax = 144.54; 
  int striprow=-1;

  if(((radius-StriprMin)<=0)||((radius-StriprMax)>=0))
    {
       if((radius-StriprMin)<=0)
	 striprow=0;

       if((radius-StriprMax)>=0)
	 striprow=255;
    }
  else
    {
      //0.04= half strip width;
      //0.4= strip pitch
      striprow= (radius -(StriprMin-0.04))/0.4;
    }
  
  return striprow;
}










//Added in 29 March 2010 - Mirko Berretti
int T2Geometry::getNearestPadRow_(double radius, double phi) {
  
  // Calculate row of pads
  int pRow = (int) floor(log(radius/(rMinPlane+0.09))/cAlpha);  
  // Small corrections if we are on the "border" or hit by a rounding mistake...
  if (pRow>=pRows)      pRow = pRows-1;
  if (pRow<0)           pRow = 0;
  return pRow;

} 

//Added in 29 March 2010 - Mirko Berretti
//Strip col from 0-1; phi in radians
int T2Geometry::getNearestStripCol_(double phi,unsigned int theid) {

  //phiTotPlane;
  
  int stripcol=-1;
  T2ROGeometry t2rogeo(theid);
  double phiminc0=t2rogeo.GetStripPhiMin(0,0);//riga-colonna
  double phimaxc0=t2rogeo.GetStripPhiMax(0,0);//riga-colonna

  double phiminc1=t2rogeo.GetStripPhiMin(0,1);//riga-colonna
  double phimaxc1=t2rogeo.GetStripPhiMax(0,1);//riga-colonna
  
  //T2R0Geo ritorna solo valori da 0 a 360 gradi ma puo' succedere che: 
  //ad ex phimin=358 phimax=96   ...  phimin=264  phimax=1.34953
  double dphic0=10.;
  double dphic1=10.;
  if(phimaxc0<phiminc0)
    {      
      //traslo tutti i phi di - phimax 
      //264 min 1.3 max  phi=1.4 
      double newphimax=360.0;//=phimaxc0-phimaxc0=0;
      double newphimin=phiminc0-phimaxc0;
      double newphi=phi-phimaxc0;
      
      if(newphi<0)
	newphi=newphi+360.;

      if((newphi<newphimax)&&(newphi>newphimin))
	{
	  //std::cout<<"Info getNearestStripCol_: stripcol=0 because col0 phistrip min|max= "<<phiminc0<<" | "<<phimaxc0<<" translated to "<<newphimin<<" | "<<newphimax<<".  Hit phi from "<<phi<<" to "<<newphi<<std::endl;
	  
	  stripcol=0;
	}
      else //Indeed can happen: Phi: 1.44689 But:   Min Max C0: 264-1.34953 |||   Min Max C1: 1.6195-96
	{
	  dphic0=fabs(phi-phimaxc0);
	  if(fabs(phi-phiminc0)<dphic0)
	    dphic0=fabs(phi-phiminc0);
	}
    }

  if(phimaxc1<phiminc1)
    {      
      //traslo tutti i phi di - phimax
      double newphimax=360.0;//=phimaxc0-phimaxc0;
      double newphimin=phiminc1-phimaxc1;
      double newphi=phi-phimaxc1;
      
      if(newphi<0)
	newphi=newphi+360.;
      
      if((newphi<newphimax)&&(newphi>newphimin))
	{
	  //std::cout<<"Info getNearestStripCol_: stripcol=1 because col1 phistrip min|max= "<<phiminc1<<" | "<<phimaxc1<<" translated to "<<newphimin<<" | "<<newphimax<<".  Hit phi from "<<phi<<" to "<<newphi<<std::endl;
	  stripcol=1;
	}
      else
	{
	  dphic1=fabs(phi-phimaxc1);
	  if(fabs(phi-phiminc1)<dphic1)
	    dphic1=fabs(phi-phiminc1);
	}
    }

  

  if((phi<phimaxc0)&&(phi>phiminc0))
    stripcol=0;
  
  if((phi<phimaxc1)&&(phi>phiminc1))
    stripcol=1;



  if(stripcol==-1)
    if((dphic0<0.5)||(dphic1<0.5))
      {
	if(dphic0<dphic1)
	  stripcol=0;
	else
	  stripcol=1;
      }
    else
      std::cout<<"ERROR in getNearestStripCol_!!  Phi: "<<phi<<" But:    "<< phiminc0<<"-"<< phimaxc0<<" |||   "<< phiminc1<<"-"<< phimaxc1<<" "<<dphic0<<" "<<dphic0<<std::endl; 

  
  //ERROR in getNearestStripCol_!!  Phi: 1.4016 But:    264-1.34953 |||   1.6195-96 262.598 262.598

  return stripcol;
}


//Added in 29 March 2010 - Mirko Berretti
//Pad col from 0-64; phi in radians

int T2Geometry::getNearestPadRow2_(Local3DPoint* hitPos,unsigned int theid)
{
  
  T2ROGeometry t2rogeo(theid);
  std::vector<double> ParRminV;std::vector<double> ParRmaxV;
  ParRminV.push_back(42.5500); ParRmaxV.push_back(44.6770);
  ParRminV.push_back(44.7770); ParRmaxV.push_back(47.0190);
  ParRminV.push_back(47.1190); ParRmaxV.push_back(49.4830);
  ParRminV.push_back(49.5830); ParRmaxV.push_back(52.0760);
  ParRminV.push_back(52.1760); ParRmaxV.push_back(54.8060);
  ParRminV.push_back(54.9060); ParRmaxV.push_back(57.6770);
  ParRminV.push_back(57.7770); ParRmaxV.push_back(60.7000);
  ParRminV.push_back(60.8000); ParRmaxV.push_back(63.8800);
  ParRminV.push_back(63.9800); ParRmaxV.push_back(67.2270);
  ParRminV.push_back(67.3270); ParRmaxV.push_back(70.7260);

  ParRminV.push_back(70.8260); ParRmaxV.push_back(74.4550);
  ParRminV.push_back(74.5550); ParRmaxV.push_back(78.3560);
  ParRminV.push_back(78.4560); ParRmaxV.push_back(82.4160);
  ParRminV.push_back(82.5160); ParRmaxV.push_back(86.7810);

  ParRminV.push_back(86.8810); ParRmaxV.push_back(91.3270);
  ParRminV.push_back(91.4270); ParRmaxV.push_back(96.1110);
  ParRminV.push_back(96.2110); ParRmaxV.push_back(101.1450);
  ParRminV.push_back(101.2450); ParRmaxV.push_back(106.4430);
  ParRminV.push_back(106.5430); ParRmaxV.push_back(112.0110);
  ParRminV.push_back(112.1110); ParRmaxV.push_back(117.8780);
  ParRminV.push_back(117.9780); ParRmaxV.push_back(124.0500);
  ParRminV.push_back(124.1500); ParRmaxV.push_back(130.5470);
  ParRminV.push_back(130.6470); ParRmaxV.push_back(137.3840);
  ParRminV.push_back(137.4840); ParRmaxV.push_back(144.5800);




  double x = hitPos->x();
  double y = hitPos->y();
  double r = sqrt(x*x + y*y);


    int rindex=-1;
    bool RIndexFound=false;
    for(unsigned int ii=0;ii<ParRmaxV.size();ii++){
      if((r<ParRmaxV.at(ii))&&(r>ParRminV.at(ii))){
	rindex=ii;
	RIndexFound=true;
	break;
      }
    }
    if(RIndexFound==false){// hitR in the spacer, check it.. could happen for cls>1 hits.
      for(unsigned int ii=0;ii<ParRmaxV.size();ii++){
	if((ii+1)<ParRminV.size())
	  if((r<ParRminV.at(ii+1))&&(r>ParRmaxV.at(ii))){
	    rindex=ii;
	    RIndexFound=true;
	    break;
	  }
      }
      
    }    

  
  return rindex;
}


int T2Geometry::getNearestPadCol_(double R,double phi,unsigned int theid)
{

  int pCol=-1;

  

  //int DetHVysign=geoutil.HV_Ysign(theid);
  
  T2ROGeometry t2rogeo(theid);

  int row=getNearestPadRow_(R, phi);
   
  bool foundcol=false;
  unsigned int j=0;
  
  while((foundcol==false)&&(j<65))
      {
	double phiminj=t2rogeo.GetPadPhiMin(row,j);//riga-colonna
	double phimaxj=t2rogeo.GetPadPhiMax(row,j);//riga-colonna
	//	std::cout<<"??? row: "<<row<<"  Phi: "<<phi<<" Phimin: "<<phiminj<<" Phimax: "<<phimaxj<<std::endl;
	//There cpuld be a gap of about 0.1 between pad. Try to look for this->enlarge pad azimuthal area.
	if((phi>=(phiminj-0.1))&&(phi<=(phimaxj+0.1)))
	  {
	    pCol=j;
	    foundcol=true;
	  }
	else
	  {
	    if(fabs(phiminj-phimaxj)>355.) // This can happen: Phimin: 358.534 Phimax: 1.41393
	      {
		if((phi<1.5)||(phi>358.5))//((phi<2.)||(phi>358.))
		  {
		     pCol=j;
		     foundcol=true;
		  }
	      }
	  }
	j++;
      }
  
  /*
    if(pCol==-1)
      {
	std::cout<<"Error in getNearestPadCol_ : col not found !!!!. Report of the problem:";
	T2GeometryUtil geoutil;
	T2GeometryUtil::T2DetInfo detinfo=geoutil.GetT2Info(theid);
	
	
	bool foundcol=false;
	unsigned int j=0;
	while((foundcol==false)&&(j<65))
	  {
	    double phiminj=t2rogeo.GetPadPhiMin(row,j);//riga-colonna
	    double phimaxj=t2rogeo.GetPadPhiMax(row,j);//riga-colonna
	    std::cout<<"??? row: "<<row<<"  Hit-Phi: "<<phi<<" Phimin: "<<phiminj<<" Phimax: "<<phimaxj<<std::endl;
	    if((phi>=phiminj)&&(phi<=phimaxj))	      	
	      foundcol=true;
	      


	    j++;
	  }
	
      }
  */

  return pCol;

}

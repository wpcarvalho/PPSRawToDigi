/** 
 * Class T2GeometryUtil. Utilities for T2 detector position / identification
 * Here CMSConventions are used especially their 
 * coordinate system
 * 
 * Author: Mirko Berretti
 * Email: mirko.berretti@gmail.com
 *
 */
#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"

T2GeometryUtil::T2GeometryUtil()
{}

T2GeometryUtil::~T2GeometryUtil()
{}

unsigned int T2GeometryUtil::RawtoSymb(uint32_t thedet)
{

  unsigned int pl=T2DetId::plane(thedet);
  unsigned int pls=T2DetId::planeSide(thedet);
  unsigned int ht=T2DetId::halfTelescope(thedet);
  unsigned int arm=T2DetId::arm(thedet);
  
  //code in the monitor
  //simbolicDetid=(*TrkCit).GetHitT2(j).GetHitPlane()*4+(*TrkCit).GetHitT2(j).GetHitHalftele()*2+(*TrkCit).GetHitT2(j).GetHitPlaneSide()+20*(*TrkCit).GetHitT2(j).GetHitArm();
  
  //hitv2.at(a).GetHitPlane()*2+hitv2.at(a).GetHitPlaneSide() -------------WAY OF THE RAWDATA
  //unsigned int symbolic=pl*4+ht*2+pls+20*arm;               -------------WAY OF THE ALIGMENT SW (OLD) 

  unsigned int symbolic=pl*2+pls+ht*10+20*arm;	  
  
  return symbolic;
}


double T2GeometryUtil::Zposition(unsigned int myarm, unsigned int myhalftele, unsigned int myplane, unsigned int myplaneside)
{
  double z1=13817;//13828.3;          //14035.605; //first Gem first drift gas zone (mm)
  double planedist=  91;//91;//86.0;//
  double btbdist=24.6;   //25.0;
  double ovdist=46;//43.0;
  double zinsidedet= 0;//1.5;  //PUT 0!!

double zdetshift=myplane*planedist+myplaneside*btbdist;

// ht==0 innner; ht==1 outer respect to the ring center.

 if((myhalftele==0)&&(myarm==0))
   zdetshift=zdetshift+ovdist;
 
 if((myhalftele==1)&&(myarm==1))
   zdetshift=zdetshift+ovdist;


 double Zdet=zdetshift+z1+zinsidedet;

 if(myarm==1)
   Zdet=Zdet*(-1.0);

 return Zdet;
}



double T2GeometryUtil::Zposition(uint32_t cmsswid)
{
  unsigned int myarm = T2DetId::arm(cmsswid);
  unsigned int myhalftele = T2DetId::halfTelescope(cmsswid);
  unsigned int myplane = T2DetId::plane(cmsswid);
  unsigned int myplaneside = T2DetId::planeSide(cmsswid);
  double Zdet=T2GeometryUtil::Zposition(myarm, myhalftele, myplane, myplaneside);
 
  return Zdet;
}


  

T2GeometryUtil::T2DetInfo T2GeometryUtil::GetT2Info(uint32_t detectid)
{
  T2GeometryUtil::T2DetInfo thet2plane;
  if(detectid>100)//argument detectid=cmsswid
    {
      thet2plane.Zdet= T2GeometryUtil::Zposition(detectid);
      thet2plane.arm=T2DetId::arm(detectid);
      thet2plane.ht=T2DetId::halfTelescope(detectid); 
      thet2plane.pl=T2DetId::plane(detectid); 
      thet2plane.plside=T2DetId::planeSide(detectid);
      thet2plane.pl_0to9=thet2plane.pl*2+thet2plane.plside;
      thet2plane.symb=RawtoSymb(detectid); 
      thet2plane.cmsswid=detectid;
    }
  
  if(detectid<100)//argument detectid=symbid
    {  
      //std::cout<<"I will assign"<<detectid/20<<std::endl;
       thet2plane.arm=detectid/20;
       thet2plane.ht=(detectid%20)/10;
       thet2plane.pl=(detectid%10)/2;  
       thet2plane.plside=(detectid%10)%2;
       thet2plane.pl_0to9=thet2plane.pl*2+thet2plane.plside;
       thet2plane.symb=detectid;
       thet2plane.cmsswid=T2DetId::calculateRawId(thet2plane.arm, thet2plane.ht,thet2plane.pl, thet2plane.plside);
       thet2plane.Zdet= T2GeometryUtil::Zposition(thet2plane.cmsswid);	 
    }
 return thet2plane;  
}








std::string T2GeometryUtil::DetInfoFromRawId(uint32_t detidr)
{
  std::string description;
  
  description="";

  unsigned int ss=RawtoSymb(detidr);
  // std::cout<<"DetInfoFromRawId Initialized with: "<<ss<<std::endl;
  if((ss/20)==0)
    {
      description+="Arm: PLUS     ";
    }
  else
    {
       description+="Arm: MINUS     ";
    }

  if(((ss/10)%2)==0)
    {
      description+="Side: NEAR     ";
    }
  else
    {
      description+="Side: FAR     ";
    }
  
 
  char str1[10];
//  int i = sprintf(str1, "%d", ss%10);

  string hs_str(str1);
  
  description+="      Horse Shoe: ";

  description+=hs_str;

  return description;
}












T2GeometryUtil::vfatid_channel T2GeometryUtil::PadVfatsIdFromRowCol(int row, int col, uint32_t detid)
{
  
  T2GeometryUtil::vfatid_channel toret;
  unsigned int vfchannel;
//  unsigned int pls=(RawtoSymb(detid)%2);

  //std::cout<<"PadVfatsIdFromRowCol row-col-pls: "<<row<<" "<<col<<" "<<pls<<std::endl;    
      
 //11-22-1      
  unsigned int vfpadiid=(col/5);
  
  int therelcol=(col*24-(vfpadiid*120))/24; //from 0 to 4
      
  //  std::cout<<"Pad vfat relative column="<<therelcol<<std::endl;
  
  // if(vfpadiid>6)
  /*
  if(pls!=0)
    {
      therelcol=4-therelcol;
    }
  */
	 
  int relativepadnumb=row+therelcol*24;
  vfchannel=3+relativepadnumb; //FIRST DATA CHANNEL IS 3 IN THE MONITOR, 4 IN THE MANUAL

  toret.vfatiid=(vfpadiid+2);
  toret.channel=vfchannel;
	
	//vfatsId.push_back(vfpadnumb);
      /*
      
      if (channelNR<4 || channelNR>123) {
      // std::cout << " ERROR: Channel (pad) number " << channelNR << " out of limits! " << std::endl;
      return false;
      }
      elementFlag = fPad;
      referencePad = (vfatIID - 2) * 120 + 97;
      columnOfSector = (Int_t)TMath::Floor((channelNR - 4) / 24);
      stepInColumn = (channelNR - 4) % 24;
      elementNR = referencePad - 24 * columnOfSector + stepInColumn;
    
      //
      //NOTE:
      //Output elementNR goes from 1 to 1560. (First pad column go from 1 to 24 )
      //

      */

  
      
  return toret;
}

std::vector<T2GeometryUtil::vfatid_channel> T2GeometryUtil::PadVfatsIdsFromPadVect(T2Hit hit)
{

  std::vector<vfatid_channel> toRets;

  std::vector<cluster_entry> entriespadcl= hit.ClusterPad_entries;
  unsigned int size=entriespadcl.size();
  hit.GetHitY();
  uint32_t detid=(hit.GetHitDetRawId());
  HV_Ysign(detid);
//  double signYHitHV=(Y*signhv);
  // std::cout << " a";
  for(unsigned int i=0;i<size;i++)
    {
      cluster_entry cle=entriespadcl.at(i);
      //std::cout << " b";
      int row=cle.rad_coord;
      int col=cle.ang_coord;
      vfatid_channel onevfid_ch=PadVfatsIdFromRowCol(row,col,detid);
      toRets.push_back(onevfid_ch);
    }

    std::vector<vfatid_channel> NodoubleCount;
  
  if(toRets.size()>0)
    NodoubleCount.push_back(toRets.at(0));
 
  if(toRets.size()>1)
    for(unsigned int i=1;i<toRets.size();i++)
      {

	bool alreadyinfinal=false;
	for(unsigned int h=0;h<NodoubleCount.size();h++)
	  {

	    if((NodoubleCount.at(h).vfatiid=toRets.at(i).vfatiid))
	      alreadyinfinal=true;
	    
	  }

	if(alreadyinfinal==false)
	  NodoubleCount.push_back(toRets.at(i));
	
      }

  return NodoubleCount;


}


T2GeometryUtil::vfatid_channel T2GeometryUtil::PadVfatsIdFromRecoHit(T2Hit hit)
{

  T2GeometryUtil::vfatid_channel toret;

  T2Geometry geom;
  geom.setPlane(hit.GetHitDetRawId());
  double R=hit.GetHitR();
  double Phi=hit.GetHitPhi();//Note: Hit Phi is by default in degree.//*3.141592/180.0;
  //  double Y=hit.GetHitY();
  int row= geom.getNearestPadRow_(R, Phi);
  int col= geom.getNearestPadCol_(R,Phi,hit.GetHitDetRawId());
  //std::cout<<"Pad Hit R-Phi: "<<hit.GetHitR()<<"-"<<hit.GetHitPhi()<<std::endl;

  if((row>=0)&&(col>=0))
    {
      
      toret=PadVfatsIdFromRowCol(row,col,(hit.GetHitDetRawId()));
    } 
  else
    {
      string ds= DetInfoFromRawId(hit.GetHitDetRawId());
      std::cout<<"! Error in T2GeometryUtil::PadVfatsIdFromRecoHit row="<<row<<" col="<<col<<std::endl;
      std::cout<<"Pad Hit R-Phi: "<<hit.GetHitR()<<"-"<<hit.GetHitPhi()<<" Detector :"<<ds.c_str()<<"  row-col "<<row<<"-"<<col<<" vfat:"<<toret.vfatiid<<" Channel:"<<toret.channel<<" Plane:"<<RawtoSymb((hit.GetHitDetRawId()))<<" PadCLS:"<<hit.GetHitNumPad()<<" StripCLS:"<<hit.GetHitNumStrip()<<std::endl;

      toret.vfatiid=-1;
      toret.channel=-1;
			 
    }

   
/*
  string ds= DetInfoFromRawId(hit.GetHitDetRawId());
  std::cout<<"Hit R-Phi: "<<hit.GetHitR()<<"-"<<hit.GetHitPhi()<<" Detector :"<<ds.c_str()<<" vfat:"<<toret.vfatiid<<" Channel:"<<toret.channel<<std::endl;
*/
  return toret;
}




T2GeometryUtil::vfatid_channel  T2GeometryUtil::StripVfatsIdFromRowCol(int row, int col, uint32_t detid)
{
  
  T2GeometryUtil::vfatid_channel toret;
  unsigned int vfchannel;
//  unsigned int pls=(RawtoSymb(detid)%2);
  //std::cout<<"StripVfatsIdFromRowCol row-col-pls: "<<row<<" "<<col<<" "<<pls<<std::endl;    
      
      if(row<128)//Can be 0 or 16
	{
	  if(col==0)
	    {	      
	      toret.vfatiid=0;
	      vfchannel=127-(row%127);        //This is like the monitor and the T2GeoWork
	      /*
	      if(pls==0)                     //if(signHitHV>0)
		vfchannel=row%127;
	      else
		vfchannel=127-(row%127);
	      */
	    }
	  else
	    {
	      toret.vfatiid=16;
	      vfchannel=row%127;
	      
	      /*
	      if(pls==0)  
		vfchannel=row%127;
	      else
		vfchannel=127-(row%127);	   
	      */
	    }
	}
      else
	{
	  if(col==0)
	    {
	     toret.vfatiid=1;
	     vfchannel=127-(row%128);
	     
	     /*
	      if(pls==0)  
		vfchannel=row%128;
	      else
		vfchannel=128-(row%128);
	     */
	    }
	  else
	    {
	      vfchannel=row%128;
	      toret.vfatiid=15;
	      /*
	      if(pls==0)  
		vfchannel=row%128;
	      else
		vfchannel=128-(row%128);
	      */
	    }
	}
      
      //I want vfat channel as in the MONITOR data (vfat controller): strip:0-127; Pad:3-122//->All 0-127 
     
      // vfchannel++;
      toret.channel=vfchannel;
   

      return toret;
}




T2GeometryUtil::vfatid_channel T2GeometryUtil::StripVfatsIdFromRecoHit(T2Hit hit)
{

  T2GeometryUtil::vfatid_channel toret;
  T2Geometry geom;
  geom.setPlane(hit.GetHitDetRawId());
  double R=hit.GetHitR();
  double Phi=hit.GetHitPhi();//*3.141592/180.0;
  hit.GetHitY();
  //int row= geom.getNearestStripRow_(R,hit.GetHitDetRawId());
  int row= geom.getNearestStripRow0_256(R,hit.GetHitDetRawId());
  int col= geom.getNearestStripCol_(Phi,hit.GetHitDetRawId());
  // std::cout<<"$$$ row col returned here: "<<row<<"  "<<col<<std::endl;
   if((row>=0)&&(col>=0))
    {
      // uint32_t detid=(hit.GetHitDetRawId());
      /*
      double signHitHV=(Y*HV_Ysign(detid));
      
      if(signHitHV<0)
	if(col==0)
	  std::cout<<"Warning T2GeometryUtil: strip (col=0) not compatible with Y*HV_Ysign<0!"<<std::endl;
      */
      toret=StripVfatsIdFromRowCol(row,col,hit.GetHitDetRawId());
    }
   else
     {
     std::cout<<"! Error in T2GeometryUtil::StripVfatsIdFromRecoHit row= "<<row<<"   col="<<col<<std::endl;
     std::cout<<"Strip Hit R-Phi: "<<hit.GetHitR()<<"-"<<hit.GetHitPhi()<</*" Detector :"<<ds.c_str()<<*/"  row-col "<<row<<"-"<<col<<" vfat:"<<toret.vfatiid<<" Channel:"<<toret.channel<<std::endl;
     toret.vfatiid=-1;
     }

   //string ds= DetInfoFromRawId(hit.GetHitDetRawId());
 
    
  
  return toret;
}




std::vector<T2GeometryUtil::vfatid_channel> T2GeometryUtil::StripVfatsIdsFromStripVect(T2Hit hit)
{
  std::vector<cluster_entry> entriesstripcl= hit.ClusterStrip_entries;
  unsigned int size=entriesstripcl.size();
  std::vector<vfatid_channel> toRets;
  //double Y=hit.GetHitY();
  uint32_t detid=(hit.GetHitDetRawId());

  //double signHitHV=(Y*HV_Ysign(detid));
  for(unsigned int i=0;i<size;i++)
    {
      cluster_entry cle=entriesstripcl.at(i);
      
      int row=cle.rad_coord;
      int col=cle.ang_coord;
      vfatid_channel onevfid_ch=StripVfatsIdFromRowCol(row,col,detid);
      toRets.push_back(onevfid_ch);      
    }

  std::vector<vfatid_channel> NodoubleCount;
  
  if(toRets.size()>0)
    NodoubleCount.push_back(toRets.at(0));
 
  if(toRets.size()>1)
    for(unsigned int i=1;i<toRets.size();i++)
      {

	bool alreadyinfinal=false;
	for(unsigned int h=0;h<NodoubleCount.size();h++)
	  {

	    if(NodoubleCount.at(h).vfatiid=toRets.at(i).vfatiid)
	      alreadyinfinal=true;
	    
	  }

	if(alreadyinfinal==false)
	  NodoubleCount.push_back(toRets.at(i));
	
      }

  return NodoubleCount;
}




//Return +-1, depending on the HV Y sign;
//Never yet used, to check the PLSmultipliefactor.
int T2GeometryUtil::HV_Ysign(uint32_t thedetID)
{
  int hvpos=0;
  unsigned int symb=RawtoSymb(thedetID);
  unsigned int quarterid=symb/10;
//  unsigned int planeid=symb%10;
  unsigned int planeside=symb%2;
  
  int PLSmultiplierfactor;
  if(planeside==0)
    PLSmultiplierfactor=1;
  if(planeside==1)
    PLSmultiplierfactor=-1;

  switch (quarterid)
    {
    case 0:
      hvpos=(-1)*PLSmultiplierfactor; //near plus has HS0 with HV down. -> plane0 is -1
      break;
    case 1:
      hvpos=(+1)*PLSmultiplierfactor;
      break;
    case 2:
      hvpos=(+1)*PLSmultiplierfactor; //near minus has HS0 with HV up. -> plane0 is +1
      break;   
    case 3:
      hvpos=(-1)*PLSmultiplierfactor;
      break;   
      
    }
  return hvpos;
 

}

/**
 * Class T2RecHit   Gives same monitor results for Class1 hits!
 *
 * Author: Mirko Berretti / University of Siena 
 * Email:  mirko.berretti@gmail.com
 * Date:   2007-12-08
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"


#include "RecoTotemT1T2/T2RecHit/interface/T2DetHitReconst.h"
#include "RecoTotemT1T2/T2RecHit/interface/T2RecHit.h"


#include "DataFormats/T2Digi/interface/T2StripDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2PadDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2StripDigi.h"
#include "DataFormats/T2Digi/interface/T2PadDigi.h"
#include "DataFormats/T2Cluster/interface/T2StripClusterCollection.h"
#include "DataFormats/T2Cluster/interface/T2PadClusterCollection.h"
#include "DataFormats/T2Hit/interface/T2Hit.h"
#include "DataFormats/T2Hit/interface/T2HitCollection.h"
#include "DataFormats/T2Hit/interface/T2HitCollectionMapping.h"
#include "DataFormats/T2Hit/interface/T2PadStripAssociator.h"
#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"
#include <iostream>
#include <fstream>

double T2RecHit::myround(double d)
{
  return floor(d + 0.5);
}

void T2RecHit::beginJob()
{
  ParRminV.clear();
  ParRmaxV.clear();

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

  /*
    Explicit recalculation from T2R0Geometry
  44.6770	42.5500
  47.0190	44.7770
  49.4830	47.1190
  52.0760	49.5830
  54.8060	52.1760
  57.6770	54.9060
  60.7000	57.7770
  63.8800	60.8000
  67.2270	63.9800
  70.7260	67.3270
  
  74.4550	70.8260
  78.3560	74.5550
  82.4160	78.4560
  86.7810	82.5160

  91.3270	86.8810
  96.1110	91.4270
  101.1450	96.2110
  106.4430	101.2450
  112.0110	106.5430
  117.8780	112.1110
  124.0500	117.9780
  130.5470	124.1500
  137.3840	130.6470
  144.5800	137.4840
  */



  checkdispsize=false;

  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo;



  if((DXdisp.size()==40)&&(DYdisp.size()==40)&&(InsertAlignmentbyCFG))//Corrections from cfg
    {
      checkdispsize=true;
      //std::cout<<"Some displacements ... "<<std::endl;
    }
  else
    {
      if(InsertAlignmentbyCFG==true)
	if(useTXTfile==false)
	  std::cout<<"Wrong size of displacements: DXdisp.size()-DYdisp.size()"<<DXdisp.size()<<"-"<<DYdisp.size()<<std::endl;
    }

  if(useTXTfile)
    checkdispsize=true;


if(checkdispsize)
    {

      if(InsertAlignmentbyCFG)//convert vector to a  map
	{            
	  if(useTXTfile==false)
	    {
	      for(unsigned int k=0;k<DYdisp.size();k++)
		{
		  T2AlignmentCorrection thet2alcorr(DXdisp.at(k),DYdisp.at(k),0.);
		  T2CorrectionMap.Add(k,thet2alcorr);
		}
	    }
	  else
	    {


	      //Read directly from alignment output. Addeed 26 12 09	      
	      std::string line;
	      char *cmsswPath = getenv("CMSSW_BASE");
	      if (cmsswPath && inputFileNameMisal[0] != '/')      
		inputFileNameMisal= string(cmsswPath) + string("/src/") + inputFileNameMisal;
	      ifstream myfile(inputFileNameMisal.c_str());
	      int detid=-1;
	      double dx=0.,dy=0.;
	      double qtiltx=0., qtilty=0.;     
	      double qX0shift=0., qY0shift=0.;
	      int quarterid=-1;
      
	      if (myfile.is_open())
		{
		  while (! myfile.eof() )
		    {
		      getline (myfile,line);
		      //cout << line << endl;
		      bool founddet=false;
		      bool founddx=false;
		      bool founddy=false;
		      bool foundTiltX=false;
		      bool foundTiltY=false;
		      bool foundShiftX0=false;
		      bool foundShiftY0=false;
		      bool foundquarter=false;

		      istringstream in(line);
		      
		      std::vector<std::string> tokens; // Create vector to hold our words
		      std::string oneword;
		      //scan the line, look at your Tags
		      while(in>>oneword)
			{	
			  tokens.push_back(oneword);
			}
		      
		      
		      for(unsigned int m=0;m<tokens.size();m++)
			{
			  //"parola" each time is overwritten
			  //expected format: DetId:<spaces>i<spaces>DX:<spaces>a<spaces>DY:<spaces>b<spaces> .....
			  
			  oneword=tokens.at(m);
			  
			  bool verbosity=false;
			  if(verbosity)
			    std::cout<<oneword<<std::endl;
			  
			  if(oneword=="DetId:")
			    {
			      detid = atoi(tokens.at(m+1).c_str());
			      if((detid>=0)&&(detid<40))
				founddet=true;
			    }
			  if(oneword=="DX:")
			    {
			      //dx=strtod(tokens.at(m+1).c_str());
			      dx= boost::lexical_cast<double>(tokens.at(m+1));
			      if(verbosity)
				std::cout<<"dx= "<<dx<<std::endl;
			      
			      founddx=true;
			    }
			  if(oneword=="DY:")
			    {	
			      dy= boost::lexical_cast<double>(tokens.at(m+1));
			      //dy=strtod(tokens.at(m+1).c_str());
			      if(verbosity)
				  std::cout<<"dy= "<<dy<<std::endl;
			      
			      founddy=true;
			    }
			  
			  //Part for quarter tilt/shift
			  if(oneword=="QuarterId:")
			    {
			      quarterid= atoi(tokens.at(m+1).c_str());
			      if((quarterid>=0)&&(quarterid<4))
				foundquarter=true;
			      else
				std::cout<<"Invalid quarter number: id must be from 0 to 3"<<std::endl;
			      
			      if(verbosity)
				std::cout<<"quarterId= "<<quarterid<<std::endl;
			    }

			  if(oneword=="TiltX:")
			    {
			      qtiltx= boost::lexical_cast<double>(tokens.at(m+1));	    			
			      if(verbosity)
				std::cout<<"TiltX= "<<qtiltx<<std::endl;
			      foundTiltX=true;
			    }

			  if(oneword=="TiltY:")
			    {
			      qtilty= boost::lexical_cast<double>(tokens.at(m+1));	    			
			      if(verbosity)
				std::cout<<"TiltY= "<<qtilty<<std::endl;
			      foundTiltY=true;
			    }

		    
			  if(oneword=="ShiftX0:")
			    {
			      qX0shift= boost::lexical_cast<double>(tokens.at(m+1));	    			
			      if(verbosity)
				std::cout<<"ShiftX0= "<<qX0shift<<std::endl;
			      foundShiftX0=true;
			    }
		    
			  if(oneword=="ShiftY0:")
			    {
			      qY0shift= boost::lexical_cast<double>(tokens.at(m+1));	    			
			      if(verbosity)
				std::cout<<"ShiftY0= "<<qY0shift<<std::endl;
			      foundShiftY0=true;
			    }
	
			  
			}
		      
		      // Here it is assumed that DX and DY are the quantity to add to the x-y position of Psimhit.
		      
		      if((founddet)&&(founddx)&&(founddy))
			{
			  std::cout<<"Correcting for Internal Det-Dx-Dy= "<<detid<<"|"<<dx<<"|"<<dy<<std::endl;
			  //T2AlignmentCorrection thet2alcorr(-dx,-dy,0.);
			  //minus sign because I look on how it works T2Alignment code corrections
			  //Alignment found -DetDisplacement, while in order to correct you have to add to the Hit the Detector Displ.
			  T2AlignmentCorrection thet2alcorr(-dx,-dy,0.); 
			  T2CorrectionMap.Add(detid,thet2alcorr);			  
			}

		      if(foundquarter)
			{
			  unsigned int DETid_;
			  double zdet;
			  double firstplanez;
			  double toinsert;
			  for(unsigned int j=0;j<10;j++)
			    {
			      planeinfo=conv.GetT2Info(quarterid*10);
			      firstplanez=planeinfo.Zdet;
			      DETid_=j+quarterid*10;
			      planeinfo=conv.GetT2Info(DETid_);
			      zdet=planeinfo.Zdet;
			      //Note T2AlignmentCorrections::Add, add a translation to the previous one, 
			      //if the detId has been already inserted.
			      
			      
			      if(foundShiftX0)
				{
				  toinsert=qX0shift*(-1.0);
				  T2AlignmentCorrection thet2alcorr(toinsert,0.,0.);
				  T2CorrectionMap.Add(DETid_,thet2alcorr);
				  if(verbosity)
				    std::cout<<"Added the X common shift= "<<toinsert<<" in detId"<<DETid_<<std::endl;
				}
			      //fabs(zdet-firstplanez) inserted on 16 June 2010
			      if(foundTiltX)
				{
				  toinsert=(-1.0)*qtiltx*fabs(zdet-firstplanez);
				  T2AlignmentCorrection thet2alcorr(toinsert,0.,0.);
				  T2CorrectionMap.Add(DETid_,thet2alcorr);	
				  if(verbosity)
				    {
				      std::cout<<DETid_<<" | Tiltdx= "<<toinsert<<std::endl;
				      
				      T2AlignmentCorrections::const_iterator it =T2CorrectionMap.find(DETid_);
				      if (it != T2CorrectionMap.end()) 
					{
					  DDTranslation forTest=it->second.Translation(); 
					  std::cout<<"Check : total x shift"<<forTest.x()<<std::endl;
					}
				      else
					std::cout<<"Error in Hit shifter"<<std::endl;
				    }

				}
			      
			      if(foundShiftY0)
				{
				  toinsert=qY0shift*(-1.0);
				  T2AlignmentCorrection thet2alcorr(0.,toinsert,0.);
				  T2CorrectionMap.Add(DETid_,thet2alcorr);	
				  if(verbosity)
				    std::cout<<"Added the Y common shift= "<<toinsert<<" in detId"<<DETid_<<std::endl;
				}

			      if(foundTiltY)
				{
				  toinsert=(-1.0)*qtilty*fabs(zdet-firstplanez);
				  T2AlignmentCorrection thet2alcorr(0.,toinsert,0.);
				  T2CorrectionMap.Add(DETid_,thet2alcorr);	
				  if(verbosity)
				    std::cout<<DETid_<<" | Tiltdy= "<<toinsert<<std::endl;
				}
			    }
			}		      		    
		      

		    }
		  myfile.close();
		}
	      else cout << "Warning T2RecHit: Unable to open misalignent file"; 			      	      
	      
	    }
	}

    }

}



T2RecHit::T2RecHit(const edm::ParameterSet& paraSet):RecHitProdParSet_(paraSet) 
{
  
  includeClass0=paraSet.getParameter<bool>("IncludeClass0Hits");
  
  InsertAlignmentbyCFG=paraSet.getParameter<bool>("InsertAlignmentbyCFG");

  Cl1MaxPad=paraSet.getParameter<unsigned int>("Cl1MaxPad");
  Cl1MaxStrip=paraSet.getParameter<unsigned int>("Cl1MaxStrip");

  DXdisp=paraSet.getParameter<std::vector<double> >("DXdisp");
  DYdisp=paraSet.getParameter<std::vector<double> >("DYdisp");
  useTXTfile=paraSet.getParameter<bool>("useTXTfile");
  inputFileNameMisal= paraSet.getUntrackedParameter<std::string>("inputFileNameMisal");
  CorrectWithResolution=paraSet.getParameter<bool>("CorrectWithResolution");
  verbosity= paraSet.getUntrackedParameter<bool>("verbosity",false);
  LabelproductInstanceName=paraSet.getParameter<std::string>("LabelproductInstanceName");
  ModuleLabelInput=paraSet.getParameter<std::string>("ModuleLabelInput");

  produces<T2HitCollection>(LabelproductInstanceName);
  produces<T2HitCollectionMapping>(LabelproductInstanceName+"Map");
  produces<T2PadStripAssociator>(LabelproductInstanceName+"MapPadonStrip");
  //  produces<T2HitCollection>("T2Hits");
  //produces<T2HitCollectionMapping>("T2HitsMap");
  //produces<T2PadStripAssociator>("T2HitsMapPadonStrip");    

}


 

void T2RecHit::CorrectHitPosition(T2Hit* ahit){

  planeinformation=conv.GetT2Info((*ahit).GetHitDetRawId());
  unsigned int t2symbId=planeinformation.symb;
  //  unsigned int curarm=planeinformation.arm;
  //unsigned int curht=planeinformation.ht;

  double r;
  double phi;
  //Check if the key is a symb ID!!
  //  std::cout<<"Correct det "<<t2symbId<<std::endl;
  T2AlignmentCorrections::const_iterator it =T2CorrectionMap.find(t2symbId);
  if (it != T2CorrectionMap.end()) 
    {

      R_m = it->second.Rotation(); 
      S_m = it->second.Translation();
      DDTranslation v((*ahit).GetHitX(), (*ahit).GetHitY(), 0);  //Initial position
      
      //std::cout<<"Before X-Y "<<(*ahit).GetHitX()<<"||"<<(*ahit).GetHitY()<<"  arm,ht= "<<curarm<<","<<curht<<std::endl;
      //
      // Check this formula under some assumption.
      //The assumption is: Rot and shift are displ of detector respect to its nominal position
      //The order is important: here rotation is assumed respect to the non shifted (shifted-corrected) position 
      v = v + S_m; 
      //v=(R_m.Inverse())*v;


      if(CorrectWithResolution==false){      
	(*ahit).SetHitX(v.x());
	(*ahit).SetHitY(v.y());
	r=sqrt(v.x()*v.x()+v.y()*v.y());
	(*ahit).SetHitR(r);
	//In T2DetHitRec: x-y=Rcos-sin phi 
	//Not consistent with CMS but only with Hit X definition starting from phi
	
	double expY_2=v.y();
	double expX_2=v.x();
	
	double expPhi_2=atan(fabs(expY_2)/fabs(expX_2));
	expPhi_2=expPhi_2*180.0/3.14159265;				    
	if((expY_2<0)&&(expX_2>0))
	  expPhi_2=360.0-expPhi_2;

	if((expY_2>0)&&(expX_2<0))
	  expPhi_2=180.0-expPhi_2;
	//expPhi_2=expPhi_2+90.;
	
	if((expY_2<0)&&(expX_2<0))
	  expPhi_2=expPhi_2+180.;
	
	phi=expPhi_2;
	(*ahit).SetHitPhi(phi);
      }
      else{

	r=sqrt(v.x()*v.x()+v.y()*v.y());


	double drchange=r-(*ahit).GetHitR();
	double expY_2=v.y();
	double expX_2=v.x();
	
	double expPhi_2=atan(fabs(expY_2)/fabs(expX_2));
	expPhi_2=expPhi_2*180.0/3.14159265;	
			    
	if((expY_2<0)&&(expX_2>0))
	  expPhi_2=360.0-expPhi_2;
	
	if((expY_2>0)&&(expX_2<0))
	  expPhi_2=180.0-expPhi_2;
	
	if((expY_2<0)&&(expX_2<0))
	  expPhi_2=expPhi_2+180.;
	
	phi=expPhi_2;

	//Initialization
	double newdigiR=(*ahit).GetHitR();
	double newdigiphi=(*ahit).GetHitPhi();

	//if(false)
	//PAD POSITION CORRECTION 

	bool reproduceDigi=false;
	
	if(reproduceDigi){
	if((*ahit).GetHitNumPad()>0)
	  {

	    double phichange=phi-(*ahit).GetHitPhi();
	    T2ROGeometry t2rogeo((*ahit).GetHitDetRawId());
	
	    //A- extract Pad Row from the two vector ParRminV, PadRmaxV
	    int rindex=0;
	    bool RIndexFound=false;
	    for(unsigned int ii=0;ii<ParRmaxV.size();ii++){
	      if(((*ahit).GetHitR()<ParRmaxV.at(ii))&&((*ahit).GetHitR()>ParRminV.at(ii))){
		rindex=ii;
		RIndexFound=true;
		break;
	      }
	    }
	    if(RIndexFound==false){// hitR in the spacer, check it.. could happen for cls>1 hits.
	      for(unsigned int ii=0;ii<ParRmaxV.size();ii++){
		if((ii+1)<ParRminV.size())
		  if(((*ahit).GetHitR()<ParRminV.at(ii+1))&&((*ahit).GetHitR()>ParRmaxV.at(ii))){
		    rindex=ii;
		    RIndexFound=true;
		    break;
		  }
	      }
	      
	    }

	    if(RIndexFound==false){
	      std::cout<<"Warning T2RecHitAlign correction -> Pad Row for "<<(*ahit).GetHitR()<<" not found. Correction not implemented"<<std::endl;	  
	    }else{
	      //B- extract Pad Col from the two vector ParRminV, PadRmaxV
	      bool PhiIndexFound=false;
	      int colindex=0;
	      for(unsigned int coll=0;coll<65;coll++){
		//double myphimin=t2rogeo.GetPadPhiMin(rindex,coll);
		//double myphimax=t2rogeo.GetPadPhiMax(rindex,coll);
		if(((*ahit).GetHitPhi()<t2rogeo.GetPadPhiMax(rindex,coll))&&((*ahit).GetHitPhi()>t2rogeo.GetPadPhiMin(rindex,coll))){
		  colindex=coll;
		  PhiIndexFound=true;	
		  break;
		}
		if(((*ahit).GetHitPhi()<t2rogeo.GetPadPhiMin(rindex,coll))&&((*ahit).GetHitPhi()>t2rogeo.GetPadPhiMax(rindex,coll))){
		  colindex=coll;
		  PhiIndexFound=true;	
		  break;
		}
		
	      }


	      if(PhiIndexFound==false)   //LOOK in the GAPS
		for(unsigned int coll=0;coll<65;coll++){
		  if(coll+1<65){
		    if(((*ahit).GetHitPhi()<t2rogeo.GetPadPhiMin(rindex,coll))&&((*ahit).GetHitPhi()<t2rogeo.GetPadPhiMax(rindex,coll)))
		      if(((*ahit).GetHitPhi()>t2rogeo.GetPadPhiMin(rindex,coll+1))&&((*ahit).GetHitPhi()>t2rogeo.GetPadPhiMax(rindex,coll+1)))
			{
			  colindex=coll;
			  PhiIndexFound=true;
			  break;
			  //std::cout<<"There is a gap between A: |"<<((*ahit).GetHitPhi()<t2rogeo.GetPadPhiMin(rindex,coll))<<" "<<((*ahit).GetHitPhi()<t2rogeo.GetPadPhiMax(rindex,coll))<<"|   and  !"<<((*ahit).GetHitPhi()>t2rogeo.GetPadPhiMin(rindex,coll+1))<<" "<<((*ahit).GetHitPhi()>t2rogeo.GetPadPhiMax(rindex,coll+1))<<". FIXED"<<std::endl;

			}

		    if(((*ahit).GetHitPhi()>t2rogeo.GetPadPhiMin(rindex,coll))&&((*ahit).GetHitPhi()>t2rogeo.GetPadPhiMax(rindex,coll)))
		      if(((*ahit).GetHitPhi()<t2rogeo.GetPadPhiMin(rindex,coll+1))&&((*ahit).GetHitPhi()<t2rogeo.GetPadPhiMax(rindex,coll+1)))
			{
			  colindex=coll;
			  PhiIndexFound=true;
			  break;
			  //std::cout<<"There is a gap between A: |"<<((*ahit).GetHitPhi()<t2rogeo.GetPadPhiMin(rindex,coll))<<" "<<((*ahit).GetHitPhi()<t2rogeo.GetPadPhiMax(rindex,coll))<<"|   and  !"<<((*ahit).GetHitPhi()>t2rogeo.GetPadPhiMin(rindex,coll+1))<<" "<<((*ahit).GetHitPhi()>t2rogeo.GetPadPhiMax(rindex,coll+1))<<". FIXED"<<std::endl;
			}

		    
		  }
		}



	      //Pad Col for 359.974  Pad Col for 359.976 359.996 0.0259705

	     
	      if(PhiIndexFound==false){
		std::cout<<"Warning T2RecHitAlign correction -> Pad Col for "<<(*ahit).GetHitPhi()<<" not found. Correction not implemented"<<std::endl;	  
		
		
	      }else{
		double modulephichange=fabs(phichange);
		newdigiphi=(*ahit).GetHitPhi();

		double PadphiDistance=fabs(t2rogeo.GetPadPhiMax(rindex,colindex)-t2rogeo.GetPadPhiMin(rindex,colindex));
		if(colindex>0)
		  PadphiDistance=fabs(t2rogeo.GetPadPhiMax(rindex,colindex)-t2rogeo.GetPadPhiMax(rindex,colindex-1));
		else
		  PadphiDistance=fabs(t2rogeo.GetPadPhiMax(rindex,colindex)-t2rogeo.GetPadPhiMax(rindex,colindex+1));

		double numSemistepPhiDouble=modulephichange/(PadphiDistance/2.0);

		unsigned int numSemistepPhi=0;
		if(numSemistepPhiDouble>1){
		  numSemistepPhi=(unsigned int) myround(numSemistepPhiDouble);
		  //if semid0<1->0 if semiD>1->0
		  //numSemistepPhi=2*(numSemistepPhiDouble/2);
		  
		}

		double phisign=0.;
		if(modulephichange>0)
		  phisign=phichange/modulephichange;
		

		if(numSemistepPhi>0){
		  //First Way:
		  //newdigiphi=newdigiphi+phisign*numSemistepPhi*(PadphiDistance/2.0);
		  //Second Way:
		  newdigiphi=newdigiphi+phisign*myround(modulephichange/PadphiDistance)*PadphiDistance;
		  //std::cout<<"IdealCorrection: "<<phichange<<"Ideal Phi-R:"<<phi<<" "<<r<<". Correcting from "<<(*ahit).GetHitPhi()<<" to "<<newdigiphi<<" Actual Dphi:"<<PadphiDistance<<std::endl;
		}

		if((*ahit).GetHitNumStrip()==0){
		  //update ALSO pad-only-cluster R
		  double moduleRchange=fabs(drchange);
		  newdigiR=(*ahit).GetHitR();
		  //This time is different: DR change with R
		  double dr=0.;
		  
		  //unsigned int numSemistepR=0;
		  dr=ParRmaxV.at(rindex)-ParRminV.at(rindex);
		  
		  if(moduleRchange>dr/2.){
		    double rsign=0.;
		    if(moduleRchange>0)
		      rsign=drchange/moduleRchange;

		    unsigned int runningRowIndex=rindex;

		    if(rsign>0){
		      
		      while((runningRowIndex<24)&&(moduleRchange>dr/2.)){			
			newdigiR=newdigiR+dr/2.;
			moduleRchange=moduleRchange-dr/2.;
			runningRowIndex++;
			dr=ParRmaxV.at(runningRowIndex)-ParRminV.at(runningRowIndex);
		      }

		    }else{
		      while((runningRowIndex>0)&&(moduleRchange>dr/2.)){			
			newdigiR=newdigiR-dr/2.;
			moduleRchange=moduleRchange-dr/2.;
			runningRowIndex++;
			dr=ParRmaxV.at(runningRowIndex)-ParRminV.at(runningRowIndex);
		      }
		    }


		  }
		}

	      }
	    }

	  }



	if((*ahit).GetHitNumStrip()>0){
	  //digitize the strip R coordinate
	  
	  //	  double StripWidth = 0.080; 
	  double StripPitch = 0.400;
	  //double StripMinR = 42.460;
	  double moduleRchange=fabs(drchange);//drchange=hitoldr-r	  	  
	  
	  newdigiR=(*ahit).GetHitR();

	  if(moduleRchange>0.5*StripPitch)
	    {
	      double signdr=0.;
	      if(drchange>0)
		signdr=1.;
	      else
		signdr=-1.;
	      

	      //First Way
	      // newdigiR=newdigiR+numRsteps*signdr*(StripPitch/2.0);
	      //Second Way
	      newdigiR=newdigiR + signdr*myround(moduleRchange /StripPitch)* StripPitch;   
	      //newdigiR=newdigiR+numRsteps*signdr*(StripPitch/2.0);
	    }
	  
	}
	}

	//(*ahit).SetHitR(r);
	//(*ahit).SetHitPhi(newdigiphi);
	//(*ahit).SetHitX(r*cos(newdigiphi*3.14159265/180.0));
	//(*ahit).SetHitY(r*sin(newdigiphi*3.14159265/180.0));
	

	//	(*ahit).SetHitR(newdigiR);
	//	(*ahit).SetHitPhi(newdigiphi);
	//	(*ahit).SetHitX(newdigiR*cos(newdigiphi*3.14159265/180.0));
	//	(*ahit).SetHitY(newdigiR*sin(newdigiphi*3.14159265/180.0));
	
	//(*ahit).SetHitX(newdigiR*cos((*ahit).GetHitPhi()*3.14159265/180.0));
	//(*ahit).SetHitY(newdigiR*sin((*ahit).GetHitPhi()*3.14159265/180.0));


	
	(*ahit).SetHitR(r);
	(*ahit).SetHitX(r*cos((*ahit).GetHitPhi()*3.14159265/180.0));
	(*ahit).SetHitY(r*sin((*ahit).GetHitPhi()*3.14159265/180.0));
	
	  /*

	  
	  -according to phichange sign decide if look at bigger col (j) or smaller col
	  bool loopcontinue=true;
	  internalloopcount=0;
	  while(loopcontinue)
{
internalloopcount++;
	  evalaute closerDigiPadStep=(phi(internalloopcount+j)-phi(internalloopcount+j pm 1))/2;
	  To fire on this pad phichange should be > closerDigiPadStep/2
	  if( fabs(phichange)> closerDigiPadStep/2)
	  {
	  digiphichange=digiphichange+sign(phichange)*closerDigiPadStep/2;
	  //calculate residual phi step
	  phichange=digiphichange-fabs(closerDigiPadStep/2);
	  }
	  else
	  loopcontinue=false;
}
	 //Update the hit with phi=digiphichange
	  
	  */


	//Only R feel the misalignment

	   
      }
 
      /*
	CorrectWithResolution
      if((v.y())>=0)
	phi = acos ((v.x())/r) * 180.0 / 3.14159265;          //acos: Principal arc cosine of x, in the interval [0,pi] radians.
      else
	{
	  if((v.x())>=0)
	    phi = 360.0 - acos ((v.x())/r) * 180.0 / 3.14159265;
	  else
	    phi = 180.0 + acos ((v.x())/r) * 180.0 / 3.14159265;
	}
	//phi = 180.0 + acos ((v.x())/r) * 180.0 / 3.14159265;     
	*/
      
	
      /*
      //if((phi<200.)&&(phi>150.))
      //	{
      //	  std::cout<<"Warn: phi= "<<phi<<" expX_2 | expY_2 = "<<expX_2<<" | "<<expY_2<<std::endl;
      //	}
      */
      //  std::cout<<"After X-Y "<<(*ahit).GetHitX()<<"||"<<(*ahit).GetHitY()<<"  arm,ht= "<<curarm<<","<<curht<<std::endl;
    }
  // symbdetid=ahit.GetHitPlane()*2+ahit.GetHitPlaneSide()+ahit.GetHitHalftele()*10+20*ahit.GetHitArm();
  //ahit.SetHitX(ahit.GetHitX()+DXdisp.at(symbdetid));
  //ahit.SetHitY(ahit.GetHitY()+DYdisp.at(symbdetid));
}





void T2RecHit::produce(edm::Event& ev, const edm::EventSetup& evSet) {

  //get the digitization information from the event
  edm::Handle<T2StripClusterCollection> t2strclcoll;
  //ev.getByLabel("T2MCl","T2StripClusters",t2strclcoll);
  ev.getByLabel(ModuleLabelInput,"T2StripClusters",t2strclcoll);
  
  edm::Handle<T2PadClusterCollection> t2padclcoll;
  // ev.getByLabel("T2MCl","T2PadClusters",t2padclcoll);
  ev.getByLabel(ModuleLabelInput,"T2PadClusters",t2padclcoll);
  


 edm::ESHandle<T2AlignmentCorrections> alignments; //smart pointer to T2AlignmentCorrections


  unsigned long int IndexPositionIntheT2Hits=0;

 
  t2padclcoll.product();
  t2strclcoll.product();

  
  auto_ptr<T2HitCollection> theT2Hits (new T2HitCollection);  //Container for all the event hits
 
  //cout << "T2RecHit2 Clusters Take ... OK" << endl;
  auto_ptr<T2HitCollectionMapping> theT2HitsMap (new T2HitCollectionMapping);  
  
  auto_ptr<T2PadStripAssociator> theT2MapPadonStrip (new T2PadStripAssociator);
  
  T2PadClusterCollection::const_iterator itpad0 = t2padclcoll->begin();
  // Case for plane with only pads ON
  if(includeClass0==true)
   while(itpad0 != t2padclcoll->end())
    {
      
      T2DetId myT2Det0 = itpad0->first;
      if((t2strclcoll->count(myT2Det0))==0)
	{
	  //  std::cout<<"A "<<myT2Det0.arm()<<myT2Det0.halfTelescope()<<myT2Det0.plane()<<myT2Det0.planeSide()<<"   | ";
	  std::vector<T2Cluster> padClv0 = itpad0->second;
	  T2Hit ahit;
	  // boost::shared_ptr<T2DetHitReconst> theT2ClusterMatchingA(new T2DetHitReconst(itpad0->first.rawId(),Cl1MaxPad,Cl1MaxStrip));
	  //change on 3/9/2010
	  boost::scoped_ptr<T2DetHitReconst> theT2ClusterMatchingA(new T2DetHitReconst(itpad0->first.rawId(),Cl1MaxPad,Cl1MaxStrip));
	  
	    //theT2ClusterMatching = new T2DetHitReconst(itpad0->first.rawId(),Cl1MaxPad,Cl1MaxStrip);
	   theT2ClusterMatchingA->AddPadClusters(padClv0);	   	   
	   theT2ClusterMatchingA->FindHits();
	  
           T2HitCollection t2planehits;
	   t2planehits=theT2ClusterMatchingA->TakePlaneHits();
	  
	   
	   for(unsigned int k=0;k<t2planehits.size();k++)
	     {
	       
	       ahit=t2planehits[k];
	       //Alignment correction
	       if(checkdispsize==true)
		{
		  CorrectHitPosition(&ahit);		  	
		}
	       	       
	       IndexPositionIntheT2Hits=theT2Hits->size();	       

	       ahit.Hit_PosInCollection=IndexPositionIntheT2Hits;
	       // theT2HitsMap[ahit.HitUniqueId]=IndexPositionIntheT2Hits;
	       theT2HitsMap->insert(std::pair< std::pair<long int,long int>, unsigned long int>(ahit.HitUniqueId,IndexPositionIntheT2Hits) );
	       theT2Hits->push_back(ahit);		
	     }

	}
      

      itpad0++;
    }

  // Case for plane with only strips ON
  T2StripClusterCollection::const_iterator itstrip0 = t2strclcoll->begin();
  if(includeClass0==true)
   while(itstrip0 != t2strclcoll->end())
    {
           
      T2DetId myT2Det0 = itstrip0->first;
      if((t2padclcoll->count(myT2Det0))==0)
	{
	  //  std::cout<<"B "<<myT2Det0.arm()<<myT2Det0.halfTelescope()<<myT2Det0.plane()<<myT2Det0.planeSide()<<"   | ";
	  std::vector<T2Cluster> stripClv0 = itstrip0->second;
	  T2Hit ahit;
	  // boost::shared_ptr<T2DetHitReconst> 
	  boost::scoped_ptr<T2DetHitReconst> theT2ClusterMatchingB(new T2DetHitReconst(itstrip0->first.rawId(),Cl1MaxPad,Cl1MaxStrip));
	  //theT2ClusterMatching = new T2DetHitReconst(itstrip0->first.rawId(),Cl1MaxPad,Cl1MaxStrip);
	   theT2ClusterMatchingB->AddStripClusters(stripClv0);	   	   
	   theT2ClusterMatchingB->FindHits();
	  
           T2HitCollection t2planehits;
	   t2planehits=theT2ClusterMatchingB->TakePlaneHits();
	  
	   
	   for(unsigned int k=0;k<t2planehits.size();k++)
	     {
	       ahit=t2planehits[k];
	      
	       //Alignment correction
	       if(checkdispsize==true)
		 {
		    CorrectHitPosition(&ahit);	
		 }	
	
	       IndexPositionIntheT2Hits=theT2Hits->size();
	      
	       ahit.Hit_PosInCollection=IndexPositionIntheT2Hits;
	       //theT2HitsMap[ahit.HitUniqueId]=IndexPositionIntheT2Hits;
	       theT2HitsMap->insert(std::pair< std::pair<long int,long int>, unsigned long int>(ahit.HitUniqueId,IndexPositionIntheT2Hits) );
	       theT2Hits->push_back(ahit);		
	     }

	}

      itstrip0++;
    }








  // Other cases
     
 for(T2StripClusterCollection::const_iterator itstrip = t2strclcoll->begin();
       itstrip != t2strclcoll->end(); itstrip++)
   {
     

     std::vector<T2Cluster> stripClv = itstrip->second; 
     T2DetId myT2Det = itstrip->first;
     
     
     T2PadClusterCollection::const_iterator itpad = t2padclcoll->begin();
     
     if((t2padclcoll->count(myT2Det)>0)&&(t2strclcoll->count(myT2Det)>0))//Added on 3/9/2010. Other case already resolved above.
     while(itpad != t2padclcoll->end())
       {
	 if(itpad->first.rawId()== myT2Det.rawId())
	 {
	   //   std::cout<<"C "<<myT2Det.arm()<<myT2Det.halfTelescope()<<myT2Det.plane()<<myT2Det.planeSide()<<"   | ";
	   //T2Hit ahit;
	   
	   std::vector<T2Cluster> padClv = itpad->second;
	   //  boost::shared_ptr<T2DetHitReconst> 
	   boost::scoped_ptr<T2DetHitReconst> theT2ClusterMatchingC(new T2DetHitReconst(itpad->first.rawId(),Cl1MaxPad,Cl1MaxStrip));
	    //theT2ClusterMatching = new T2DetHitReconst(itpad->first.rawId(),Cl1MaxPad,Cl1MaxStrip);
	   //  cout << " T2RecHit process plane at Z: " << theT2ClusterMatching->GetDetZ() << endl;
	   theT2ClusterMatchingC->AddStripClusters(stripClv);
	   theT2ClusterMatchingC->AddPadClusters(padClv);	   	   
	   theT2ClusterMatchingC->FindHits();
	  
           T2HitCollection t2planehits;
	   t2planehits=theT2ClusterMatchingC->TakePlaneHits();
	  	   	   	     

	   for(unsigned int k=0;k<t2planehits.size();k++)
	     {
	       
	       T2Hit ahit=t2planehits[k];
	       //cout << " Class bo hit of rawid "<< ahit.GetHitDetRawId()<<" saved"<<std::endl;
	       
	       //Alignment correction
	       if(checkdispsize==true)
		{
		   CorrectHitPosition(&ahit);		   
		}
		 
	       
	       
	       if(ahit.GetHitClass()==1) 
		 {


		   IndexPositionIntheT2Hits=theT2Hits->size();
		   (*theT2MapPadonStrip)[ahit.HitUniqueId.first].push_back(ahit.HitUniqueId.second);
		   
		   ahit.Hit_PosInCollection=IndexPositionIntheT2Hits;
		   //theT2HitsMap[ahit.HitUniqueId]=IndexPositionIntheT2Hits;
		   theT2HitsMap->insert(std::pair< std::pair<long int,long int>, unsigned long int>(ahit.HitUniqueId,IndexPositionIntheT2Hits) );
		   theT2Hits->push_back(ahit);	
		   //cout << " Class 1 hit of rawid "<< ahit.GetHitDetRawId()<<" saved"<<std::endl;
		   //cout<<"Cl1Hit:   Z:"<<ahit.GetHitZ()<<"R: "<<ahit.GetHitR()<<" Phi:"<<ahit.GetHitPhi()<<std::endl;
		   
		   //cout << " Class 1 hit in plane at Z: " << theT2ClusterMatching->GetDetZ() <<endl;
		 }
	       else
		 if(includeClass0==true)		    
		   {
		     IndexPositionIntheT2Hits=theT2Hits->size();
		     ahit.Hit_PosInCollection=IndexPositionIntheT2Hits;
		     //theT2HitsMap[ahit.HitUniqueId]=IndexPositionIntheT2Hits;
		     theT2HitsMap->insert(std::pair< std::pair<long int,long int>, unsigned long int>(ahit.HitUniqueId,IndexPositionIntheT2Hits) );
		     theT2Hits->push_back(ahit);		     
		   }
	     }
	     
	 }  
	 itpad++;

       }

     

   }

 // std::cout<< " $$$$$" <<std::endl;

 //bool ret=theT2HitsMap->insert(std::pair< std::pair<int,int>, unsigned int>(ahit.HitUniqueId,IndexPositionIntheT2Hits) ); 
 //if (ret.second==false)
 //   cout << "element 'z' already existed";

 //produces<T2HitCollection>(LabelproductInstanceName);
 //produces<T2HitCollectionMapping>(LabelproductInstanceName+"Map");
 //produces<T2PadStripAssociator>(LabelproductInstanceName+"MapPadonStrip");

 //ev.put(theT2Hits, "T2Hits");
 //ev.put(theT2HitsMap, "T2HitsMap");
 // ev.put(theT2MapPadonStrip,"T2HitsMapPadonStrip");

ev.put(theT2Hits, LabelproductInstanceName);
ev.put(theT2HitsMap, LabelproductInstanceName+"Map");
ev.put(theT2MapPadonStrip,LabelproductInstanceName+"MapPadonStrip");

} // produce

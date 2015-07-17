/**
 * Class PSimHitMisaligner
 *
 * Author: Mirko Berretti / University of Siena 
 * Email:  mirko.berretti@gmail.com
 * Date:   2007-12-08
 */

#include "FWCore/Framework/interface/ESHandle.h"
#include <iostream>
#include "SimTotem/T2Digitizer/interface/PSimHitMisaligner.h"
#include <fstream>
#include <string>

#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"


PSimHitMisaligner::PSimHitMisaligner(const edm::ParameterSet& conf,const edm::EventSetup &iSetup) 
{
  
  
  inputFileNameMisal= conf.getUntrackedParameter<std::string>("inputFileNameMisal");
  simulatemisalign= conf.getUntrackedParameter<bool>("simulatemisalign");
  generaterandom= conf.getUntrackedParameter<bool>("generaterandom");

  sigmaDispl= conf.getUntrackedParameter<double>("sigmaDispl"); //dx dy internal misal
  sigmaPhiDispl= conf.getUntrackedParameter<double>("sigmaPhiDispl"); //dphi internal misal
  sigmaGlobalThetaXY= conf.getUntrackedParameter<double>("sigmaGlobalThetaXY");
  sigmaGlobalShiftXY= conf.getUntrackedParameter<double>("sigmaGlobalShiftXY");
  verbosity= conf.getUntrackedParameter<bool>("verbosity");

 
  unsigned int expectedrawsize=11;
  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo;
  
  std::vector<double> onedetdisp;
  for(unsigned int u=0;u<expectedrawsize;u++)
    {
      onedetdisp.push_back(0.);
    }

  //Initialization
  if(simulatemisalign)
    {
      for (unsigned int j=0;j<40;j++)
	{
	  int detid=j;
	  onedetdisp.at(0)=detid/20;//arm
	  onedetdisp.at(1)=(detid%20)/10;//ht
	  onedetdisp.at(2)=(detid%10)/2;//pl
	  onedetdisp.at(3)=(detid%10)%2;//pls
	  onedetdisp.at(10)=detid;

	  onedetdisp.at(4)=0.;//Rx
	  onedetdisp.at(5)=0.;//Ry
	  onedetdisp.at(6)=0.;//Rz
	  onedetdisp.at(7)=0.;//dx
	  onedetdisp.at(8)=0.;//dy
	  onedetdisp.at(9)=0.;//Dz

	  dispvect.push_back(onedetdisp);
	  if(verbosity)
	    std::cout<<"arm-ht-pl-pls "<<onedetdisp.at(0)<<"-"<<onedetdisp.at(1)<<"-"<<onedetdisp.at(2)<<"-"<<onedetdisp.at(3)<<std::endl;
	}
    }

      //--------------------------------------
      //Read Displacemnt from a file To put also in RecHit?

  if((simulatemisalign)&&(!generaterandom)&&(inputFileNameMisal!=""))
    {

      char *cmsswPath = getenv("CMSSW_BASE");
      inputFileNameMisal= string(cmsswPath) + string("/src/") + inputFileNameMisal;

      std::string line;
      ifstream myfile (inputFileNameMisal.c_str());
      int detid=-1;
      double dx=0.,dy=0., dphi=0.;
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
		    
		    if(oneword=="DPhi:")
		      {	
		        dphi= boost::lexical_cast<double>(tokens.at(m+1));
			//dy=strtod(tokens.at(m+1).c_str());
			if(verbosity)
			  std::cout<<"dphi= "<<dphi<<std::endl;


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
		// Here DPhi will change  position of Psimhit as: v'=R(phi)v

		if((founddet)&&(founddx)&&(founddy))
		  {
		    //Arm-Ht-Pl-Pls-Rx-Ry-Rz-Dx-Dy-Dz-SymbId
		    
		    onedetdisp.at(0)=detid/20;//arm
		    onedetdisp.at(1)=(detid%20)/10;//ht
		    onedetdisp.at(2)=(detid%10)/2;//pl
		    onedetdisp.at(3)=(detid%10)%2;//pls
		    onedetdisp.at(4)=0.;//Rx
		    onedetdisp.at(5)=0.;//Ry
		    onedetdisp.at(6)=dphi;//Rz 
		    onedetdisp.at(7)=dx;
		    onedetdisp.at(8)=dy;
		    onedetdisp.at(9)=0.;//Dz
		    onedetdisp.at(10)=detid;		    
		    dispvect.at(detid)=onedetdisp;
		  }

		if(foundquarter)
		  {
		    unsigned int DETid_;
		    double zdet;
		    double firstplanez;
		    for(unsigned int j=0;j<10;j++)
		      {
			planeinfo=conv.GetT2Info(quarterid*10);
			firstplanez=planeinfo.Zdet;
			DETid_=j+quarterid*10;
			planeinfo=conv.GetT2Info(DETid_);
			zdet=planeinfo.Zdet;


			//modification of for tilt fabsDZ made on
			if(foundShiftX0)
			  {
			    dispvect.at(DETid_).at(7)=dispvect.at(DETid_).at(7)+qX0shift;
			    if(verbosity)
			     std::cout<<"Det "<<DETid_<<" digiShiftX Added, now dx="<<dispvect.at(DETid_).at(7)<<std::endl;
			  }
			if(foundTiltX)
			  {
			    dispvect.at(DETid_).at(7)=dispvect.at(DETid_).at(7)+qtiltx*fabs(zdet-firstplanez);
			    if(verbosity)
			      std::cout<<"Det "<<DETid_<<" digiTiltX Added, now dx="<<dispvect.at(DETid_).at(7)<<std::endl;
			  }

			if(foundShiftY0)
			  {
			    dispvect.at(DETid_).at(8)=dispvect.at(DETid_).at(8)+ qY0shift;
			    if(verbosity)
			      std::cout<<"Det "<<DETid_<<" digiShiftY Added, now dy="<<dispvect.at(DETid_).at(8)<<std::endl;
			  }
			if(foundTiltY)
			  {
			    dispvect.at(DETid_).at(8)=dispvect.at(DETid_).at(8)+qtilty*fabs(zdet-firstplanez);
			    if(verbosity)
			      std::cout<<"Det "<<DETid_<<" digiTiltY Added, now dy="<<dispvect.at(DETid_).at(8)<<std::endl;
			  }
		      }
		  }

	    }
	  myfile.close();
	}
      else
	{ 
		cout << "Warning PSimHitMisaligner: Unable to open file with name:"; 
		std::cout<<" "<<inputFileNameMisal.c_str()<<std::endl;
	}
    }
 




  //--------------------------------------
  //Fill random Dispacements

 if((simulatemisalign)&&(generaterandom))
   {
    
      if(verbosity)
	std::cout<<"Displacement Random Generation"<<std::endl;
     
     
     for(unsigned int arm=0;arm<2;arm++)
       for(unsigned int ht=0;ht<2;ht++)
	 // if(countdisp<NumDispInHalf)
	 //  {
	     for(unsigned int pl=0;pl<5;pl++)
	       for(unsigned int pls=0;pls<2;pls++)
		 {
		   std::vector<double> onedetdisp;
		   for(unsigned int u=0;u<expectedrawsize;u++)
		     {
		       onedetdisp.push_back(0.);
		     }
	 

		   onedetdisp.at(0)=arm;
		   onedetdisp.at(1)=ht;
		   onedetdisp.at(2)=pl;
		   onedetdisp.at(3)=pls;
		
		   onedetdisp.at(6)=CLHEP::RandGaussQ::shoot(0,sigmaPhiDispl);
		   onedetdisp.at(7)=CLHEP::RandGaussQ::shoot(0,sigmaDispl);
		   onedetdisp.at(8)=CLHEP::RandGaussQ::shoot(0,sigmaDispl);

		   if(verbosity)
		     std::cout<<"dx|dy|dphi= "<<onedetdisp.at(7)<<"|"<<onedetdisp.at(8)<<"|"<<onedetdisp.at(6)<<std::endl;
		
		   int detid = arm*20+ht*10+pl*2+pls;
		   onedetdisp.at(10)=detid;
		
		   if(verbosity)
		     std::cout<<"Saved detid= "<<detid<<std::endl;

		   dispvect.at(detid)=onedetdisp;

		 }    
	     // countdisp++;
	     //}
   }





//--------------------------------------
//Create maps displacements<->rawId

//If you don't read from DB
//if(!((simulatemisalign)&&(!generaterandom)&&(inputFileNameMisal=="")))
if((simulatemisalign)&&(!generaterandom)&&(inputFileNameMisal!=""))
  {
    if(simulatemisalign)   //Added on 1-Sept-09
      {
	T2DetId detconverter;
	uint32_t rawiddet;
	TMatrix RotmatrixX(3,3);
	TMatrix RotmatrixY(3,3);
	TMatrix RotmatrixZ(3,3);
	TMatrix RotmatrixAll(3,3);
	TMatrix Translation(3,3);


	//rot:4-5-6 shift 7-8-9  (x-y-z)
	//Arm-Ht-Pl-Pls-Dx-Dy-Dz-Rx-Ry-Rz-SymbId

  

	for(unsigned int l=0;l<dispvect.size();l++)
	  if(dispvect.at(l).size()==expectedrawsize)
	    {
	      if(verbosity)
		std::cout<<"Det "<< dispvect.at(l).at(10) <<" taken for displacements"<<" pls:"<<dispvect.at(l).at(3)<<std::endl;
	      Translation.Zero();
	      RotmatrixAll.Zero();
	      RotmatrixZ.UnitMatrix();
	      RotmatrixY.UnitMatrix();
	      RotmatrixX.UnitMatrix();
      
	      RotmatrixX(1,1)=cos(dispvect.at(l).at(4));
	      RotmatrixX(1,2)=-sin(dispvect.at(l).at(4));
	      RotmatrixX(2,1)=sin(dispvect.at(l).at(4));
	      RotmatrixX(2,2)=cos(dispvect.at(l).at(4));

	      
	      RotmatrixY(0,0)=cos(dispvect.at(l).at(5));
	      RotmatrixY(0,2)=sin(dispvect.at(l).at(5));
	      RotmatrixY(2,0)=-sin(dispvect.at(l).at(5));
	      RotmatrixY(2,2)=cos(dispvect.at(l).at(5));

	      RotmatrixZ(0,0)=cos(dispvect.at(l).at(6));
	      RotmatrixZ(0,1)=-sin(dispvect.at(l).at(6));
	      RotmatrixZ(1,0)=sin(dispvect.at(l).at(6));
	      RotmatrixZ(1,1)=cos(dispvect.at(l).at(6));


	      RotmatrixAll=RotmatrixX*RotmatrixY*RotmatrixZ;
	      Translation(0,0)=dispvect.at(l).at(7); //X
	      Translation(1,1)=dispvect.at(l).at(8); //Y	
	      Translation(2,2)=dispvect.at(l).at(9); //Z

	      

	      rawiddet=detconverter.calculateRawId((unsigned int) dispvect.at(l).at(0),(unsigned int) dispvect.at(l).at(1),(unsigned int) dispvect.at(l).at(2),(unsigned int) dispvect.at(l).at(3)); //in t2detid arm ht pl pls --- in cfi arm pl pls ht

	      T2DetMatrixRotation.insert(pair<uint32_t,TMatrix>(rawiddet,RotmatrixAll));
	      T2DetMatrixTranslation.insert(pair<uint32_t,TMatrix>(rawiddet,Translation));
	    }  

  
  

      }  
  }

  
} // END PSimHitMisaligner::PSimHitMisaligner(....)



PSimHitMisaligner::~PSimHitMisaligner() 
{
  //TFile *f = TFile::Open("./test.root", "recreate");
  //if( !f || !f->IsWritable() ){
  // std::cout << "Output file not opened correctly !!" << std::endl;
  // }
  //TestDx->Write("");

}






Local3DPoint PSimHitMisaligner::DisplacePoint(const Local3DPoint &p,uint32_t rawiddet)
{

  T2DetId detconverter;
  T2GeometryUtil conv3;
  T2GeometryUtil::T2DetInfo planeinfo3;
   
  // std::cout<<" Try to Disalign T2 hit..."<<std::endl;
  TVector hitxycoord(3);
  hitxycoord(0)=p.x();   
  hitxycoord(1)=p.y(); 
  hitxycoord(2)=0.;  

  TVector hitzcoord(3);
  hitzcoord(0)=0.;
  hitzcoord(1)=0;
  hitzcoord(2)=p.z();//You cannot get this info from local t2h->GetHitZ();

  //rawiddet=detconverter.calculateRawId(t2h->GetHitArm(),t2h->GetHitHalftele(),t2h->GetHitPlane(),t2h->GetHitPlaneSide());

  T1T2DetMatrixMap::iterator it1=T2DetMatrixRotation.find(rawiddet);
  T1T2DetMatrixMap::iterator it2=T2DetMatrixTranslation.find(rawiddet);

  //TMatrix matrixrot=T2DetMatrixRotation[rawiddet];
  //TMatrix matrixtrsl=T2DetMatrixTranslation[rawiddet];
  TVector finalhitcoord(3);
  // std::cout<<"----- Plane "<<rawiddet<<" --------------------"<<std::endl;
  if((it1!=T2DetMatrixRotation.end())&&(it2!=T2DetMatrixTranslation.end()))
    {

      planeinfo3=conv3.GetT2Info(rawiddet);
      unsigned int curht=planeinfo3.ht;
      unsigned int curarm=planeinfo3.arm;
      unsigned int curpls=planeinfo3.plside;
      unsigned int curplane=planeinfo3.pl;
      bool toprint=false;
      if((curht==0)&&(curarm==0))
	toprint=true;

      TMatrix matrixrot=(it1->second);
      TMatrix matrixtrsl=(it2->second);
  
      //Change sign of the input trasl in order to be compatible with the Erik PsimHit Conv      
      ChangeDXDYSignAccordingPsimHitErikConv(&matrixtrsl,rawiddet);
      ChangePhiSignAccordingPsimHitErikConv(&matrixrot,rawiddet);
     
      //hitxycoord.Print();
      finalhitcoord=matrixrot*hitxycoord;
      if((verbosity)&&(toprint))
	 {
	   std::cout<<"planeside-plane  "<<curpls<<"-"<<curplane<<"  Local Point Geant Hit Before"<<std::endl;
	   finalhitcoord.Print();
	 }
     
      TVector one(3);
      one(0)=1.;
      one(1)=1.;
      one(2)=1.;
      // matrixtrsl.Print();
      // hitxycoord.Print();
      finalhitcoord=finalhitcoord+hitzcoord+(matrixtrsl*one); 
      //If it assumed that the input file is MINUS the detector Displacement: The Geant hit should have this displacement.



      //hitzcoord.Print();
      //   std::cout<<"After"<<std::endl;
      //  matrixtrsl.Print();
      //  finalhitcoord.Print();
      //TestDx->Fill(hitxycoord(0)-finalhitcoord(0));

      /*
	t2h->SetHitX(finalhitcoord(0));
	t2h->SetHitY(finalhitcoord(1));
	t2h->SetHitR(sqrt(finalhitcoord(0)*finalhitcoord(0)+finalhitcoord(1)*finalhitcoord(1)));
	float phih;
	if((t2h->GetHitY()>0)&&(t2h->GetHitX()>0))
	phih=acos(t2h->GetHitX()/t2h->GetHitR())*(180.0/3.14159265);

	if((t2h->GetHitY()>0)&&(t2h->GetHitX()<0))
	phih=acos(t2h->GetHitX()/t2h->GetHitR())*(180.0/3.14159265);

	if((t2h->GetHitY()<0)&&(t2h->GetHitX()>0))
	phih=360.0-acos(fabs(t2h->GetHitX())/t2h->GetHitR())*(180.0/3.14159265);

	if((t2h->GetHitY()<0)&&(t2h->GetHitX()<0))
	phih=180.0+acos(fabs(t2h->GetHitX())/t2h->GetHitR())*(180.0/3.14159265);
      
	double r2=t2h->GetHitX();
	double drr=r1-r2;
	unsigned int symbolic=t2h->GetHitPlane()*4+t2h->GetHitHalftele()*2+t2h->GetHitPlaneSide()+20*t2h->GetHitArm();
      
	t2h->SetHitPhi(phih);

      
	t2h->SetHitZ(finalhitcoord(2));
      */
      if((verbosity)&&(toprint))
	 {
	   std::cout<<"planeside-plane  "<<curpls<<"-"<<curplane<<"  Local Point Geant Hit After"<<std::endl;
	   finalhitcoord.Print();   
	 }
      // std::cout<<"After Phi: "<<t2h->GetHitPhi() <<std::endl;
      //std::cout<<t2h->GetHitHalftele()<<t2h->GetHitPlane()<<t2h->GetHitPlaneSide()<<t2h->GetHitArm()<<std::endl;
    }

  //std::cout<<"After Phi: "<<t2h->GetHitPhi() <<std::endl;
  // else
  // {
  //   std::cout<<" Trasf. not implemented for this detector"<<std::endl;
  // }

  return Local3DPoint(finalhitcoord(0),  finalhitcoord(1), finalhitcoord(2));

}
 

PSimHit PSimHitMisaligner::DisplacePSimHit(const PSimHit &input,uint32_t rawdetid)
{

  // Local3DPoint ep = input.entryPoint(), xp = DisplacePoint(input.exitPoint(),rawdetid);
  PSimHit ahit = PSimHit(DisplacePoint(input.entryPoint(),rawdetid), DisplacePoint(input.exitPoint(),rawdetid), input.pabs(), input.tof(), input.energyLoss(), input.particleType(),input.detUnitId(), input.trackId(), input.thetaAtEntry(), input.phiAtEntry(), input.processType());
  
  return ahit;

  //return PSimHit(DisplacePoint(input.entryPoint(),rawdetid), DisplacePoint(input.exitPoint(),rawdetid), 
  //	 input.pabs(), input.tof(), input.energyLoss(), input.particleType(),
  //	 input.detUnitId(), input.trackId(), input.thetaAtEntry(), input.phiAtEntry(), input.processType());


}


void PSimHitMisaligner::ChangePhiSignAccordingPsimHitErikConv(TMatrix* matrixrot,uint32_t rawiddet){
  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinformation;
  planeinformation=conv.GetT2Info(rawiddet);
  /*
    RotmatrixZ(0,0)=cos(dphi);
    RotmatrixZ(0,1)=-sin(dphi);
    RotmatrixZ(1,0)=sin(dphi);
    RotmatrixZ(1,1)=cos(dphi);
    dphi=d(atan(Y/X))
  */

  //Look to  ChangeDXDYSignAccordingPsimHitErikConv for sign changing
  if((planeinformation.ht==0)&&(planeinformation.plside==0)) //type 2
    {
      (*matrixrot)(0,1)=((*matrixrot)(0,1))*(-1.0);
      (*matrixrot)(1,0)=((*matrixrot)(1,0))*(-1.0);
    }

  
  /*
    This is the reference plane 
  if((planeinformation.ht==0)&&(planeinformation.plside==1)) //type 1
    {
      
    }
  */

  if((planeinformation.ht==1)&&(planeinformation.plside==0)) //type 3
    {
      (*matrixrot)(0,1)=((*matrixrot)(0,1))*(-1.0);
      (*matrixrot)(1,0)=((*matrixrot)(1,0))*(-1.0);      
    }

  /*
    Double sign change. dPhi unchanged
  if((planeinformation.ht==1)&&(planeinformation.plside==1)) //type 4
    {
      
    }
  */  
}
     

//Transform CMS-coordinate in PSIM-Hit coordinate
void PSimHitMisaligner::ChangeDXDYSignAccordingPsimHitErikConv(TMatrix* matrixtrsl,uint32_t rawiddet)
{
  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo detInfo;
  detInfo=conv.GetT2Info(rawiddet);
  //if(verbosity)
  //std::cout<<"Arm: "<<planeinformation.arm<<" Plane "<<planeinformation.pl<<"  | PlSide "<<planeinformation.plside<<"  | HT "<<planeinformation.ht<<std::endl;
 
  //ht   0==Inner 1==Outer
  //pls  0==Front 1==Back (Front = "active area looking at IP5")
  //(*matrixtrsl)(0,0)=X
  //(*matrixtrsl)(1,1)=Y

  if (detInfo.ht == 1 && detInfo.plside == 1){
   
    (*matrixtrsl)(0,0) = -(*matrixtrsl)(0,0);      
     if (detInfo.arm == 0){
      (*matrixtrsl)(1,1) = -(*matrixtrsl)(1,1);       
       }
  } else if (detInfo.ht == 1 && detInfo.plside == 0){
    (*matrixtrsl)(0,0) = -(*matrixtrsl)(0,0);
    if (detInfo.arm == 1){
      (*matrixtrsl)(1,1) = -(*matrixtrsl)(1,1);       
        }
  } else if (detInfo.ht == 0 && detInfo.plside == 1){  //REFERENCE
     if (detInfo.arm == 1){
     (*matrixtrsl)(1,1) = -(*matrixtrsl)(1,1);       
      }
  } else if (detInfo.ht == 0 && detInfo.plside == 0){
    if (detInfo.arm == 0){
      (*matrixtrsl)(1,1) = -(*matrixtrsl)(1,1);       
     }
  }   

  
   /*
  if((planeinformation.ht==0)&&(planeinformation.plside==0)) //type 2
    {
      (*matrixtrsl)(1,1)= ((*matrixtrsl)(1,1))*(-1.0);
    }

 
    // This is the reference plane 
  if((planeinformation.ht==0)&&(planeinformation.plside==1)) //type 1
    {
      (*matrixtrsl)(0,0)=((*matrixtrsl)(0,0));
      (*matrixtrsl)(1,1)= ((*matrixtrsl)(1,1));
    }
  

  if((planeinformation.ht==1)&&(planeinformation.plside==0)) //type 3
    {
      (*matrixtrsl)(0,0)=((*matrixtrsl)(0,0))*(-1.0);      
    }

  if((planeinformation.ht==1)&&(planeinformation.plside==1)) //type 4
    {
      (*matrixtrsl)(0,0)=((*matrixtrsl)(0,0))*(-1.0);
      (*matrixtrsl)(1,1)= ((*matrixtrsl)(1,1))*(-1.0);
    }
   */ 

}

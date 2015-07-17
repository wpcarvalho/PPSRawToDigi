/**
 * Class to simulate the electronics of the T2 strip
 * detector.
 *
 * Author: Erik Brücken / University of Helsinki
 * Email:  brucken@cc.helsinki.fi
 * Date:   2007-11-26
 *         Mirko Berretti / University of Siena & Pisa INFN
 *         mirko.berretti@cern.ch
*/

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Utilities/General/interface/FileInPath.h"
#include "SimTotem/T2Digitizer/interface/T2StripElectronicsSim.h"
#include "SimTotem/T2Digitizer/interface/T2DetectorHit.h"
#include "DataFormats/T2Digi/interface/T2StripDigi.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandExponential.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <iostream>
#include <fstream>
#include <string>

#include "TotemCondFormats/DAQInformation/interface/AnalysisMask.h"
#include "TotemCondFormats/DataRecord/interface/TotemDAQMappingRecord.h"

/**
 *
 */


T2StripElectronicsSim::T2StripElectronicsSim(const edm::ParameterSet & parameterSet,const edm::EventSetup& iSet,CLHEP::HepRandomEngine* rndEngine)
  : bins_(parameterSet.getParameter<std::vector<int> >("bins")),
    sigmaExtraNoise_(parameterSet.getParameter<std::vector<int> >("sigmaExtraNoise")), 
    simpleThreshold_(parameterSet.getParameter<std::vector<int> >("simpleThreshold")),
    capaNoiseFactorStrip_(parameterSet.getParameter<std::vector<double> >("capaNoiseFactorStrip")), 
    UseCFGInfo(parameterSet.getParameter<bool>("UseCFGInfo")),
    UseVFATs(parameterSet.getParameter<bool>("UseVFATs")),
    SetVfatEfficiency(parameterSet.getParameter<bool>("SetVfatEfficiency")), 
    inputFileNameEffi(parameterSet.getParameter<std::string>("inputFileNameEffi")),
    inputFileNameCorrupt(parameterSet.getParameter<std::string>("inputFileNameCorrupt")),
    inputFileNameDeadSect(parameterSet.getParameter<std::string>("inputFileNameDeadSect")),
    inputFileNameDeadChannels(parameterSet.getParameter<std::string>("inputFileNameDeadChannels")),

    inputFileNameNoiseCorrelation(parameterSet.getParameter<std::string>("inputFileNameNoiseCorrelation")),
    inputFileNameNoisyChannels(parameterSet.getParameter<std::string>("inputFileNameNoisyChannels"))
{
 
  theVFatEfficiency=boost::shared_ptr<VFatEfficiency>(new VFatEfficiency(SetVfatEfficiency,inputFileNameEffi,inputFileNameCorrupt,inputFileNameDeadChannels,inputFileNameNoisyChannels));

  readFile("SimTotem/T2Digitizer/data/StripNoiseHalf_0.dat", capaNoiseSigma[0]);
  readFile("SimTotem/T2Digitizer/data/StripNoiseHalf_1.dat", capaNoiseSigma[1]);
  readFile("SimTotem/T2Digitizer/data/StripThresholdHalf_0.dat", threshold[0]);
  readFile("SimTotem/T2Digitizer/data/StripThresholdHalf_1.dat", threshold[1]);
  

  rndEngineStr=rndEngine;

  LoadDeadSector(inputFileNameDeadSect,VectDeadSect_Plane,VectDeadSect_Sector);
  
  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinformation;
  
  if(UseCFGInfo==false)
    {
      std::cout<<"Warning forced to use cfg strip threshold information Used"<<std::endl;      
      UseCFGInfo=true;
    }
  
  //if(UseCFGInfo==true)
  //  std::cout<<"I will use cfg"<<std::endl;
  //else
  //  std::cout<<"I will use db"<<std::endl;

  // VFATS Map Thresholds  And Efficiency Construction and initialization ----------------------------------------------------
  if(UseVFATs)
    {

      T2DigiVfat theVfat;

      if(theVFatEfficiency->SetVfatEfficiency)
	{
	  char *cmsswPath = getenv("CMSSW_BASE");
	  if(inputFileNameNoiseCorrelation!=""){
	    inputFileNameNoiseCorrelation= string(cmsswPath) + string("/src/") + inputFileNameNoiseCorrelation;
	    MCfileNoise = boost::shared_ptr<TFile> (new TFile(inputFileNameNoiseCorrelation.c_str())); 	  
	    EvtVfat_Strip_WhenCompletelyOn=(THnSparseD*)MCfileNoise->Get("VFAT_Monitoring/EvtVfat_Strip_WhenCompletelyOn");
	  }
	}
      
     if(UseCFGInfo==true) //Set thresholds from cfg
       {
       for(unsigned int d=0;d<bins_.size();d++)
	{

	  // std::cout<<"In stripElectSim looks";
	  planeinformation=conv.GetT2Info(d);       
	  //  std::cout<<"io gan4a2 "<<planeinformation.cmsswid<<" "<<planeinformation.symb<<std::endl;
	  T2DetId thet2det(planeinformation.cmsswid);
	  //  std::cout<<"OK"<<std::endl;
	  uint32_t thedetid=planeinformation.cmsswid;
	  //int symbid=(int) planeinformation.symb;
	  std::map<unsigned int, T2DigiVfat> them1;



	  for(unsigned int m=0;m<2;m++)
	    {
	      //Same conv used in xml reader and Dqm	   
	      theVfat= T2DigiVfat(thedetid,m,0);

	      theVfat.SetThreshold(-1,bins_.at(d));
	      if(theVFatEfficiency->SetVfatEfficiency)
		{
		  unsigned int absvfatid=d*100+m;		  
		  double eff=0.8;

		  if(((theVFatEfficiency)->EffiMap).find(absvfatid)!=((theVFatEfficiency)->EffiMap).end())
		    eff=((theVFatEfficiency)->EffiMap)[absvfatid];
		  else
		    std::cout<<"Efficiency not found for VFAT"<< m <<" in plane "<<d<<". Set to 0.8"<<std::endl;

		  /*
		  if((d<20)&&(d>=10))
		    if(eff>=0.04)
		    eff=eff-0.04;
		  */
		  
		  theVfat.Efficiency_=eff;
		  //std::cout<<"Messo "<<theVfat.Efficiency_<<std::endl;
		  //The higher is the efficiency, the lower is the equivalent threshold.		      
		  
		  double thresholddependentEffi= EqThrFromEffi(eff)/*bins_.at(d)*/ ; 
		  //double thresholddependentEffi=bins_.at(d);

		  theVfat.SetThreshold(-1,thresholddependentEffi);
		  
		  int numdeadchannel=0;
		  if(((theVFatEfficiency)->VfatID_ToDeadChannelList).find(absvfatid)!=((theVFatEfficiency)->VfatID_ToDeadChannelList).end())
		    numdeadchannel=(theVFatEfficiency->VfatID_ToDeadChannelList)[absvfatid].size();
		  
		 
		  if(numdeadchannel>0)
		    {
		      std::vector<unsigned int> listofDeadchannel=(theVFatEfficiency->VfatID_ToDeadChannelList)[absvfatid];
		      theVfat.SetDeadChannels(listofDeadchannel);
		    }
		}

	     
	      them1.insert(std::pair<unsigned int, T2DigiVfat>(m,theVfat));

	    }

	   for(unsigned int m=15;m<17;m++)
	    {
	      //std::auto_ptr<T2DigiVfat> theVfat(new T2DigiVfat(thedetid,m,0));
	      //(*theVfat).SetThreshold(-1,bins_.at(d));
	      theVfat= T2DigiVfat(thedetid,m,0);
	      theVfat.SetThreshold(-1,bins_.at(d));
	      if(theVFatEfficiency->SetVfatEfficiency)
		{
		  unsigned int absvfatid=d*100+m;		  
		  double eff=0.8;

		  if(((theVFatEfficiency)->EffiMap).find(absvfatid)!=((theVFatEfficiency)->EffiMap).end())
		    eff=((theVFatEfficiency)->EffiMap)[absvfatid];
		  else
		    std::cout<<"Efficiency not found for VFAT"<< m <<" in plane "<<d<<". Set to 0.8"<<std::endl;



		  theVfat.Efficiency_=eff;
		  //The higher is the efficiency, the lower is the equivalent threshold.
		  double thresholddependentEffi= EqThrFromEffi(eff); /*+ bins_.at(d)*(1-eff)*/; 
		  //double thresholddependentEffi=bins_.at(d);
		  theVfat.SetThreshold(-1,thresholddependentEffi);
		   

		   int numdeadchannel=0;
		   if(((theVFatEfficiency)->VfatID_ToDeadChannelList).find(absvfatid)!=((theVFatEfficiency)->VfatID_ToDeadChannelList).end())
		     numdeadchannel=(theVFatEfficiency->VfatID_ToDeadChannelList)[absvfatid].size();

		  
		  if(numdeadchannel>0)
		    {
		      std::vector<unsigned int> listofDeadchannel=(theVFatEfficiency->VfatID_ToDeadChannelList)[absvfatid];
		      theVfat.SetDeadChannels(listofDeadchannel);
		    }
		}
	      
	      
	      them1.insert(std::pair<unsigned int, T2DigiVfat>(m,theVfat));	    	   
	    }

	   StripVFats.insert(std::pair<int, std::map<unsigned int, T2DigiVfat> >(d,them1));

	}
       // std::cout<<"VFMAP loaded with size "<<StripVFats.size()<<std::endl;
       }

     // Check the size of VFATS  Map Thresholds from DB acquisition and put default values------------------------------------

     if(UseCFGInfo==true)
       {

	 // for(unsigned int vfnumb =0;vfnumb<17;d++)
	 // if((vfnumb==0)||(vfnumb==1)||(vfnumb==15)||(vfnumb==16)) 
	 for(unsigned int d=0;d<40;d++)       
	   if(StripVFats.find(d)==StripVFats.end())
		{
		  std::cout<<"Warning: Strip Thresholds in plane "<<d<<" not found. Put to default values."<<std::endl;
		  T2GeometryUtil conv;
		  T2GeometryUtil::T2DetInfo planeinformation;
		  planeinformation=conv.GetT2Info(d); 	       

		  std::map<unsigned int, T2DigiVfat> them1;
		  //T2DigiVfat theVfat;

		  for(unsigned int m=0;m<2;m++)
		    {
		      //Same conv used in xml reader and Dqm

		      //std::auto_ptr<T2DigiVfat> theVfat(new T2DigiVfat(planeinformation.cmsswid,m,0));
		      theVfat=T2DigiVfat(planeinformation.cmsswid,m,0);
		      theVfat.SetThreshold(-1,35);//DefaultThr 
		      theVfat.Efficiency_=0.8;
		      them1.insert(std::pair<unsigned int, T2DigiVfat>(m,theVfat));

		    }

		for(unsigned int m=15;m<17;m++)
		  {	      
		    // std::auto_ptr<T2DigiVfat> theVfat(new T2DigiVfat(planeinformation.cmsswid,m,0));
		    //(*theVfat).SetThreshold(-1,35);//DefaultThr
		    theVfat=T2DigiVfat(planeinformation.cmsswid,m,0);
		    theVfat.SetThreshold(-1,35);//DefaultThr 
		    theVfat.Efficiency_=0.8;
		    them1.insert(std::pair<unsigned int, T2DigiVfat>(m,theVfat));	    	   
		  }

		StripVFats.insert(std::pair<int, std::map<unsigned int, T2DigiVfat> >(d,them1));

	     }
	   else
	     {
	       //T2DigiVfat theVfat;
	       if((StripVFats[d]).size()!=4)
		 {
		   std::cout<<"Warning: Only "<<(StripVFats[d]).size()<<"strip vfat found for threshold setting:" <<std::endl;
		   T2GeometryUtil conv;
		   T2GeometryUtil::T2DetInfo planeinformation;
		   planeinformation=conv.GetT2Info(d);   

		   for(unsigned int m=0;m<17;m++)
		     {
		       if((m==0)||(m==1)||(m==15)||(m==16))
			 {
			   if((StripVFats[d]).find(m)==(StripVFats[d]).end())
			     {
			       //std::auto_ptr<T2DigiVfat> theVfat(new T2DigiVfat(planeinformation.cmsswid,m,0));
				theVfat=T2DigiVfat(planeinformation.cmsswid,m,0);
				theVfat.SetThreshold(-1,35);//DefaultThr
				theVfat.Efficiency_=0.8;
			       //std::map<unsigned int, T2DigiVfat> them1;
			       //them1.insert(std::pair<unsigned int, T2DigiVfat>(m,(*theVfat)));
			       (StripVFats[d]).insert(std::pair<unsigned int, T2DigiVfat>(m,theVfat));
			       std::cout<<m<<" Put to default values."<<std::endl;
			     }

			 }
		     }

		 }
	     }
	 //std::cout<<"VFMAP completed, final size "<<StripVFats.size()<<std::endl;
       }

      

      bool simulateNoisyChannels = false;
      if (parameterSet.exists("simulateNoisyChannels"))
      { 
        simulateNoisyChannels = parameterSet.getParameter<bool>("simulateNoisyChannels");
      }
      //Create Noisy channel mask
      if (simulateNoisyChannels) 
      {
        edm::ESHandle<AnalysisMask> analysisMask;
        iSet.get<TotemDAQMappingRecord>().get(analysisMask);

      
        //Treat noisy channels as dead channels
        unsigned int key;
        unsigned int planenumb;
        unsigned int vfnumb;
        for (map<TotemSymbID, VFATAnalysisMask>::const_iterator it = analysisMask->analysisMask.begin(); it != analysisMask->analysisMask.end(); it++)
        {
          key = it->first.symbolicID;
          planenumb=key/100;
          vfnumb=key%100;
          if((vfnumb==0)||(vfnumb==1)||(vfnumb==15)||(vfnumb==16)) //stripVFAT
          {
            if(StripVFats.find(planenumb)!=StripVFats.end())
            {
              if(StripVFats[planenumb].find(vfnumb)!=StripVFats[planenumb].end())
              {
                vector<unsigned int> vec;
                for (set<unsigned char>::iterator it2 = it->second.maskedChannels.begin(); it2 != it->second.maskedChannels.end(); it2++)
                {
                  vec.push_back(*it2);
                }
                StripVFats[planenumb][vfnumb].SetDeadChannels(vec);
              } else
                std::cout<<"Warning: one vfat-strip not initialized"<<std::endl;
            } else
             std::cout<<"Warning: one plane not initialized in strip digitization"<<std::endl;
          }
        }
      }


     

    }//end if useVFATs
   


   if(sigmaExtraNoise_.size()!=40) 
     {
       std::cout<<"T2Digi Strip Warning: Found "<<sigmaExtraNoise_.size()<<" detector sigmaExtraNoise_ insted of 40. Missing will be put to default values"<<std::endl;
      if(sigmaExtraNoise_.size()<40)
	{
	  unsigned int miss=40-sigmaExtraNoise_.size();
	  for (unsigned int i=0;i<miss;i++)
	    {
	      sigmaExtraNoise_.push_back(0);
	      // std::cout<<"Det "<<sigmaExtraNoise_.size()<<"th at 0."<<std::endl;
	    }
	}
    }


   if(capaNoiseFactorStrip_.size()!=40) 
     {
       std::cout<<"T2Digi Strip Warning: Found "<<capaNoiseFactorStrip_.size()<<" detector capaNoiseFactorStrip_ insted of 40. Missing will be put to default values"<<std::endl;
      if(capaNoiseFactorStrip_.size()<40)
	{
	  unsigned int miss=40-capaNoiseFactorStrip_.size();
	  for (unsigned int i=0;i<miss;i++)
	    {
	     capaNoiseFactorStrip_.push_back(0.);
	     //std::cout<<"Det "<<capaNoiseFactorStrip_.size()<<"th at 1."<<std::endl;
	    }
	}
    }



   if(simpleThreshold_.size()!=40) 
     {
        std::cout<<"T2Digi Strip Warning: Found "<<simpleThreshold_.size()<<" detector simpleThreshold_ insted of 40. Missing will be put to default values"<<std::endl;
      if(simpleThreshold_.size()<40)
	{
	  unsigned int miss=40-simpleThreshold_.size();
	  for (unsigned int i=0;i<miss;i++)
	    {
	      simpleThreshold_.push_back(600);
	      // std::cout<<"Det "<<simpleThreshold_.size()<<"th at 600"<<std::endl;
	    }
	}
    }



} // T2StripElectronicsSim

/**
 *
 */

double T2StripElectronicsSim::EqThrFromEffi(double effi_measured){
  
  /*
    Gaussian fit obtained with single muon at SimpleThr=400
    
    Simu have been used with fully alive and 0 corrupted vfats, no dead sector.

    Old Test
    StripWidth=0.14    diffCoeff=0.4    gain=15000    capaNoiseFactorStrip=1    simpleThreshold=400    sigmaExtraNoise=0
    Gaussian Parameter: C=1.016 mu=-62.82  sigma=284.3 where  f(x) = p0*exp(-0.5*((x-p1)/p2)^2))  
    double eqThr=0.;
    double mu=(-1.0)*3.967;
    double C=0.9905;
    double sigma=118.7;
  */

  /*
  //StripWidth=0.05    diffCoeff=0.24    gain=30000    capaNoiseFactorStrip=0.1    simpleThreshold=400    sigmaExtraNoise=0
  //Extracted Param:
  double eqThr=0.;
  double mu=(-1.0)*17.22;
  double C=0.9999;
  double sigma=138.1;
  */
  /*
  //StripWidth=0.05    diffCoeff=0.21    gain=25000    capaNoiseFactorStrip=0.1    simpleThreshold=400    sigmaExtraNoise=0
  //Extracted Param:
  double eqThr=0.;
  double mu=(-1.0)*12.13;
  double C=0.99;
  double sigma=123.3;
  */
  /*
//StripWidth=0.05    diffCoeff=0.27    gain=25000    capaNoiseFactorStrip=0.1    simpleThreshold=400P/550S    sigmaExtraNoise=0
  //Extracted Param:
  double eqThr=0.;
  double mu=(-1.0)*9.98;
  double C=0.99;
  double sigma=105.;
  double logar=0.;
  */
  /*
  //Extracted Param - CURVE E :
  double eqThr=0.;
  double mu=(-1.0)*9.495;
  double C=0.99;
  double sigma=99.5;
  double logar=0.;
  */
  
  //Extracted Param - CURVE G :
  double logar=0.;
  double eqThr=0.;
  double mu=(-1.0)*6.931;
  double C=0.99;
  double sigma=75.18;
  if(effi_measured>0.1)
    {
      if((effi_measured/C)<1.)//Normal Case, logar<0.
	{
	  logar=log(effi_measured/C);
	  eqThr=mu+sigma*sqrt(-2.0*logar);
	}
      else//Very High Effi case: put lower Thr limit.
	{
	  eqThr=5;// was 6 on 12 Jul ;
	}
    }
  else //Not very predictable at low effi, give back 3 sigma.
    {
      eqThr=3*sigma;
      if(effi_measured<0.05)
	eqThr=3000000.0*sigma;
    }

  if(eqThr<5)
    std::cout<<"Warning: low strip eq-threshold"<<eqThr<<std::endl;

  return eqThr;
  
  
}










void T2StripElectronicsSim::LoadDeadSector(std::string inputFileNameDeadSect,std::vector<unsigned int> &VectDeadSect_Plane,std::vector<unsigned int> &VectDeadSect_Sector)
{
  VectDeadSect_Plane.clear();
  VectDeadSect_Sector.clear();  
    
  char *cmsswPath = getenv("CMSSW_BASE");
  inputFileNameDeadSect= string(cmsswPath) + string("/src/") + inputFileNameDeadSect;  
  std::string line;
  ifstream myfile (inputFileNameDeadSect.c_str());
  
  
  int detid=-1;
  int sectordead=-1; 
 
  if (myfile.is_open())
	{
	while (! myfile.eof() )
	  {
	    getline (myfile,line);
	    bool founddet=false;
	    bool foundSect=false;
	    std::vector<std::string> tokens; // Create vector to hold our words
	    std::string oneword;
	    //scan the line, look at your Tags
	    istringstream in(line);
	    while(in>>oneword)
	      {	
		tokens.push_back(oneword);
	      }
	    
	    for(unsigned int m=0;m<tokens.size();m++)
	      {
		//"parola" each time is overwritten
		//expected format: DetId:<spaces>i<spaces>DX:<spaces>a<spaces>DY:<spaces>b<spaces> .....
		
		oneword=tokens.at(m);
		
		//bool verbosity=false;
		//if(verbosity)
		//std::cout<<oneword<<std::endl;

		if(oneword=="DetId:"){
		  detid = atoi(tokens.at(m+1).c_str());
		  if((detid>=0)&&(detid<40))
		    {
		      founddet=true;		      
		    }
		}

		if(oneword=="SectorDead:"){
		  sectordead = atoi(tokens.at(m+1).c_str());
		  if((sectordead>=0)&&(sectordead<=5))
		    {
		      foundSect=true;		    
		    }
		}
						
	      }
	    
	    if((foundSect)&&(founddet))
	      {
		VectDeadSect_Plane.push_back((unsigned int)detid);  
		VectDeadSect_Sector.push_back((unsigned int)sectordead);  
	      }
	    
	  }
	  
	}
  else
    std::cout<<"Dead Sect file not opened"<<std::endl;
}




unsigned int T2StripElectronicsSim::RawtoSymb(uint32_t thedet)
{
  T2DetId converter;
  unsigned int pl=converter.plane(thedet);
  unsigned int pls=converter.planeSide(thedet);
  unsigned int ht=converter.halfTelescope(thedet);
  unsigned int arm=converter.arm(thedet);
  unsigned int symbolic=pl*2+pls+ht*10+20*arm;
  return symbolic;
}


bool T2StripElectronicsSim::IsStripInDeadSector(unsigned int RVal,unsigned int symbdetid){
  bool toreturn=false;
  
  unsigned int Rsector=0;
  //RVal is 0..255
 
  if(RVal<19)
    {
      Rsector=0;
    }
  else
    {
      if(RVal<126)
	{
	  Rsector=1;
	}
      else
	{
	  if(RVal<198)
	    Rsector=2;
	  else
	    Rsector=3;
	}
    }
  
  for(unsigned int j=0;j<VectDeadSect_Plane.size();j++){
    if(VectDeadSect_Plane.at(j)==symbdetid) 
      if(VectDeadSect_Sector.at(j)==Rsector) 
	toreturn=true;
  }

  return toreturn;
}











void T2StripElectronicsSim::readFile(std::string fileName, double* dataFile) {
  
  std::string path(getenv( "CMSSW_SEARCH_PATH" ));

  FileInPath file(path, fileName);

  if (file() == 0 ) {
    
    edm::LogError("T2StripElectronicsSim") << "File with Capacitancenoise for VFat not found: " 
					   << file.name();
    throw cms::Exception(" Capacitance noise file not found!");
  } 

  LogDebug("T2StripElectronicsSim") <<"Reading in " << file.name();
  edm::LogInfo("T2StripElectronicsSim") << "-> Reading in " << file.name();
 
  std::ifstream & fin = *file();
  
  LogDebug("T2StripElectronicsSim") << "T2StripElectronicsSim: opening file " << file.name();

  if (fin == 0) {
    std::string errorMessage =  "Cannot open file " +  file.name();
    edm::LogError("T2StripElectronicsSim") << errorMessage;
    throw cms::Exception(errorMessage);
  }

  double var1;
  double var2;
  
  int channel = 0;
  
  while (fin.good() && channel<NUMBER_OF_CHANNELS) {
    
    fin >> var1 >> var2;
    
    dataFile[channel] = var2;
    
    channel++;
  }

  if (channel<NUMBER_OF_CHANNELS) {
    
    throw cms::Exception("Data file error") << "To little data in fine " << file.name();
  }
      
  fin.close();
  
} // readFile

/**
 *
 */


void T2StripElectronicsSim::simulate(std::map<int, int*> & chargeMap) {
  
  // std::cout<<"StripElectr. simulate Beg."<<std::endl;

  std::map<int, int*>::iterator i;
  for (i=chargeMap.begin(); i != chargeMap.end(); ++i) {

    unsigned int symbdetid= RawtoSymb((*i).first);
    //std::cout<<"symbdetid = "<<symbdetid<<std::endl;
    for (int sector=0; sector<2; sector++) {
      
      for(int k=0; k<NUMBER_OF_CHANNELS; k++) {

	// noise: Gauss with Sigma from capacitance noise. Mean from energy deposition

	//std::cout<<test<<" Added from strip noise"<<std::endl;
	
	(*i).second[k+sector*NUMBER_OF_CHANNELS] += static_cast<int> (
	  capaNoiseFactorStrip_.at(symbdetid)*(static_cast<int>(ceil(CLHEP::RandGaussQ::shoot( 0, capaNoiseSigma[sector][k]*simpleThreshold_.at(symbdetid)))))
	  // additional noise to fill up to meassured noise
	  +capaNoiseFactorStrip_.at(symbdetid)*static_cast<int>(ceil(CLHEP::RandGaussQ::shoot( 0, sigmaExtraNoise_.at(symbdetid)*simpleThreshold_.at(symbdetid) )))  );
	//std::cout<<" ACTUAL VALUE OF STRIP CHARGE= "<<(*i).second[k+sector*NUMBER_OF_CHANNELS]<<std::endl;
	
      }
    }
  }
  // std::cout<<"StripElectr. simulate End."<<std::endl;
} // simulate

/**
 *
 */

void T2StripElectronicsSim::fillDigis(T2StripDigiCollection & stripDigis,std::map<int, int*> & chargeMap) 
{
  //std::cout<<"Initial size in filldigi "<<StripVFats.size()<<std::endl;

  std::map<int, int*>::iterator i;

  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinformation;
  unsigned int VfatID;
  int detectorID;
  unsigned int vfatchannel;
  int vfatchannelThr;

  
  



  for (i=chargeMap.begin(); i != chargeMap.end(); ++i) {
    for (int sector=0; sector<2; sector++) {
      for(int k=0; k<NUMBER_OF_CHANNELS; k++) {

	//  (NUMBER_OF_CHANNELS 256)

	VfatID=(sector*15);
	if(sector==0)  
	  VfatID=VfatID+(k/128);
	if(sector==1)
	  if(k<=127)
	    VfatID=16;
	  else
	    VfatID=15;
	//Rember:0-16 inner; 1-15 outer.

	detectorID=(*i).first;
	vfatchannel=k%128;                  //Check if it is correct!!!!!!!! //Read data channel:4-123 only for pad
	planeinformation=conv.GetT2Info(detectorID);
	
	//	std::cout<<"Channel "<<vfatchannel<<" of Vfat n° "<<VfatID<<"In plane "<<planeinformation.symb<<" has threshold "<<vfatchannelThr<<std::endl;

	int charge = (*i).second[k+ sector*NUMBER_OF_CHANNELS];        
	//unsigned int symbdetid= RawtoSymb((*i).first);
	unsigned int symbdetid=planeinformation.symb;

	if(UseVFATs)
	  vfatchannelThr=StripVFats[symbdetid][VfatID].GetThreshold(vfatchannel);
	else
	  vfatchannelThr=bins_.at(symbdetid);



	//double electronicEfficiency=1.0; 
	double randomValue=1.;
	double corruptProb=0.;

	if(UseVFATs)
	  if(theVFatEfficiency->SetVfatEfficiency)
	    {	      
	      // electronicEfficiency= StripVFats[symbdetid][VfatID].Efficiency_; 
	      unsigned int Thisabsvfatid=symbdetid*100+VfatID;

	      if(((theVFatEfficiency)->CorruptMap).find(Thisabsvfatid)!=((theVFatEfficiency)->CorruptMap).end())
		corruptProb=((theVFatEfficiency)->CorruptMap)[Thisabsvfatid];

	      randomValue=(double)(CLHEP::RandFlat::shoot(0.,1.));	      
	      //std::cout<<"ReadEffi="<<electronicEfficiency<<" random:"<<randomValue<<std::endl;
	    }

	if(StripVFats[symbdetid][VfatID].IsChannelDead(vfatchannel)==false)	  
	  {

	     bool stripindeadsector=IsStripInDeadSector((unsigned int) k,symbdetid);
	     // Threshold = number of bins * simpleThr(e-)	
	     // 
	     //if(charge > bins_.at(symbdetid) * simpleThreshold_.at(symbdetid))	     
	     //&&(electronicEfficiency>randomValue)
	    if((charge > vfatchannelThr * simpleThreshold_.at(symbdetid)))
	      if(randomValue>corruptProb)//First electronic gate which take into account corrupt-prob and dead vfats. 
		if(stripindeadsector==false)
		  { 
		    
		    //std::cout<<"here in filldigi1 "<<std::endl;
		    T2DetId planeId((*i).first);
		    //std::cout<<"here in filldigi2 "<<std::endl;
		    T2StripDigi sDigi(0, k, sector, charge); //sector means col in rawdata and k is the row
		    stripDigis.insertDigi(planeId, sDigi);
		    // std::cout<<"Setting StripVFats["<<symbdetid<<"]["<<VfatID<<"].SetChannel("<<vfatchannel<<","<<1<<")"<<std::endl;
		    StripVFats[symbdetid][VfatID].SetChannel(vfatchannel,1);
		    //std::cout<<"symbdetid-charge-thr "<<symbdetid<<" "<<charge<<" "<<(bins_.at(symbdetid)* simpleThreshold_.at(symbdetid))<<std::endl;

		  }
	  }
//	else
//	  std::cout<<"Warn: vfat pad-channel simulated dead"<<std::endl;
      }
    }
  }    


  //----------------------------------------------------------------------
  // COMPLETELY-ON NOISE SIMULATION
  //----------------------------------------------------------------------

   unsigned int symbdetid_=0;unsigned int VfatID_=0;unsigned int absidvfat=0;
   unsigned int sect=0; unsigned int channelInPlane=0;

//   std::cout<<"File name "<<inputFileNameNoiseCorrelation.c_str()<<std::endl;
  if(UseVFATs)
    if((inputFileNameNoiseCorrelation!="")&&(EvtVfat_Strip_WhenCompletelyOn->GetEntries()>0))
      if(theVFatEfficiency->SetVfatEfficiency)
	{
	  
	  /*	  
	   edm::Service<edm::RandomNumberGenerator> rng;
	   if ( ! rng.isAvailable()) {
	     throw cms::Exception("Configuration")
	       << "This class requires the RandomNumberGeneratorService\n"
	       "which is not present in the configuration file.  You must add the service\n"
	       "in the configuration file or remove the modules that require it.";
	   }
  
	   rndEngine = &(rng->getEngine());
	   
	  */
	  unsigned int rseed=0;
	  rseed=(unsigned int)( CLHEP::RandFlat::shoot(rndEngineStr) * std::numeric_limits<unsigned int>::max() );
	  gRandom->SetSeed(rseed);

	   // std::cout<<"RSEED="<<rseed<<std::endl;
	   /*

	   //m_randomEngine = rng->getEngine();//m_rndmSvc->GetEngine(m_randomEngineName);
	   unsigned int rseed=0;
	   rseed=(unsigned int)( CLHEP::RandFlat::shoot(rndEngine) * std::numeric_limits<unsigned int>::max() );
	   gRandom->SetSeed(rseed);
	   */

	  int correlation_dimension= EvtVfat_Strip_WhenCompletelyOn->GetNdimensions();	  
	  Double_t *rand = new Double_t[correlation_dimension];
	  EvtVfat_Strip_WhenCompletelyOn->GetRandom(rand,false);
	  /*
	    std::cout<<"Plotting strip correlation, size:"<<correlation_dimension<<std::endl;
	    for(unsigned int uu=0;uu<correlation_dimension;uu++)
	    {
	    std::cout<<rand[uu]<<"-";
	    }
	    std::cout<<std::endl;
	  */
	  
	  
   	for(int oo=0;oo<correlation_dimension;oo++){
	  if(rand[oo]>0){

	    //    TAxis* currentAx=  EvtVfat_Strip_WhenCompletelyOn->GetAxis(g);
	    //const char* thetile = currentAx->GetTitle();
	    // std::cout<<"Axes Title: "<<g<<"-"<<thetile<<std::endl;
	    
	    TAxis* currentAx=  EvtVfat_Strip_WhenCompletelyOn->GetAxis(oo);
	    const char* thetile = currentAx->GetTitle();
	    // std::cout<<"!!Ax Title: "<<oo<<"-"<<thetile<<std::endl;
	    absidvfat = atoi(thetile);
	    //  unsigned int absidvfat=tostring(currentAx->GetTitle(oo));
	    //unsigned int absidvfat=1010;
	    symbdetid_=absidvfat/100;
	    VfatID_=absidvfat%100;
	    StripVFats[symbdetid_][VfatID_].SetChannel(-1,1);
	  
	    //------Convert symbdetid_ to cmsswID
	    planeinformation=conv.GetT2Info(symbdetid_);
	    T2DetId planeId(planeinformation.cmsswid);
	    //Rember:0-16 inner; 1-15 outer.
	    if((VfatID_==15)||(VfatID_==16))
	      {
		sect=1;
			  
		for(int k=0; k<128; k++)//NUMBER_OF_CHANNELS=256 in a sector
		  {
		    channelInPlane=k;
		    if(VfatID_==15)
		      channelInPlane=channelInPlane+128;
		    T2StripDigi sDigi(0, channelInPlane, sect, 10000); //sector means col in rawdata and k is the row
		    stripDigis.insertDigi(planeId, sDigi);
		     
		  }
	      }
	    else{
	      sect=0;
	      
	      for(int k=0; k<128; k++)//NUMBER_OF_CHANNELS=256 in a sector
		  {
		    channelInPlane=k;
		    if(VfatID_==1)
		      channelInPlane=channelInPlane+128;
		    T2StripDigi sDigi(0, channelInPlane, sect, 10000); //sector means col in rawdata and k is the row
		    stripDigis.insertDigi(planeId, sDigi);
		     
		  }
	    }
	    
	    //std::cout<<"Setting StripVFats["<<symbdetid_<<"]["<<VfatID_<<"]"<<" arm-ht"<<planeId.arm()<<" "<<planeId.halfTelescope()<<" completely ON"<<std::endl;

	    
	    //std::cout<<"Setting StripVFats["<<symbdetid<<"]["<<VfatID<<"].SetChannel("<<vfatchannel<<","<<1<<")"<<std::endl;
	    
	    //std::cout<<"symbdetid-charge-thr "<<symbdetid<<" "<<charge<<" "<<(bins_.at(symbdetid)* simpleThreshold_.at(symbdetid))<<std::endl;
	  }
	}
     
      delete[] rand;

      
    }

 
 // std::cout<<"Final size in filldigi "<<StripVFats.size()<<std::endl;


} // fillDigis

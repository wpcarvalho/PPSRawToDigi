#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Utilities/General/interface/FileInPath.h"
#include "SimTotem/T2Digitizer/interface/T2PadElectronicsSim.h"
#include "SimTotem/T2Digitizer/interface/T2DetectorHit.h"
#include "DataFormats/T2Digi/interface/T2PadDigi.h"
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
 



T2PadElectronicsSim::T2PadElectronicsSim(const edm::ParameterSet & parameterSet,const edm::EventSetup& iSet)
  : bins_(parameterSet.getParameter<std::vector<int> >("bins")),
    sigmaExtraNoise_(parameterSet.getParameter<std::vector<int> >("sigmaExtraNoise")),
    simpleThreshold_(parameterSet.getParameter<std::vector<int> >("simpleThreshold")),
    capaNoiseFactorPad_(parameterSet.getParameter<std::vector<double> >("capaNoiseFactorPad")), 


    SetVfatEfficiency(parameterSet.getParameter<bool>("SetVfatEfficiency")), 
    inputFileNameEffi(parameterSet.getParameter<std::string>("inputFileNameEffi")),
    inputFileNameCorrupt(parameterSet.getParameter<std::string>("inputFileNameCorrupt")),
    inputFileNameDeadSect(parameterSet.getParameter<std::string>("inputFileNameDeadSect")),
    inputFileNameDeadChannels(parameterSet.getParameter<std::string>("inputFileNameDeadChannels")),

    inputFileNameNoiseCorrelation(parameterSet.getParameter<std::string>("inputFileNameNoiseCorrelation")),
    inputFileNameNoisyChannels(parameterSet.getParameter<std::string>("inputFileNameNoisyChannels")), 
    UseCFGInfo(parameterSet.getParameter<bool>("UseCFGInfo")),
    UseVFATs(parameterSet.getParameter<bool>("UseVFATs"))
{
  theVFatEfficiency=boost::shared_ptr<VFatEfficiency>(new VFatEfficiency(SetVfatEfficiency,inputFileNameEffi,inputFileNameCorrupt,inputFileNameDeadChannels,inputFileNameNoisyChannels));


  readFile("SimTotem/T2Digitizer/data/PadsENC.dat", capaNoiseSigma[0]);
  //Pad capacitance from padindex: 0 to 120 (5 periods)


  //  readFile("SimTotem/T2Digitizer/data/PadThreshold.dat", threshold[0]);
  //std::cout<<"NumDeadSect:"<<VectDeadSect_Plane.size()<<std::endl;
  
  LoadDeadSector(inputFileNameDeadSect,VectDeadSect_Plane,VectDeadSect_Sector);
  

  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinformation;
  
  if(UseCFGInfo==false)
    {
      std::cout<<"Warning forced to use cfg pad threshold information Used"<<std::endl;      
      UseCFGInfo=true;
    }

  
  //std::cout<<"NumDeadSect:"<<VectDeadSect_Plane.size()<<std::endl;

  // VFATS Map Thresholds And Efficiency Construction and initialization ----------------------------------------------------
  if(UseVFATs)
    {
      
      T2DigiVfat theVfat;

      if(UseCFGInfo==true) //Set thresholds from cfg
	{
	  if(bins_.size()!=40)
	    std::cout<<"Error in T2PadElectronicSim.cc (a)"<<std::endl;
	  else
	  for(unsigned int d=0;d<40;d++)
	    {
	      //std::cout<<"In stripElectSim looks";
	      planeinformation=conv.GetT2Info(d);       
	      T2DetId thet2det(planeinformation.cmsswid);
	      // std::cout<<"OK"<<std::endl;
	      uint32_t thedetid=planeinformation.cmsswid;
	      std::map<unsigned int, T2DigiVfat> them1;
	      
	      for(unsigned int m=2;m<15;m++)
		{
		  theVfat= T2DigiVfat(thedetid,m,0);
		  
		  //Initialization of the common threshold (-1)
		  theVfat.SetThreshold(-1,bins_.at(d));
		  //Final Thr bin*simplethr where bin~40, simplethr~400
		  if(theVFatEfficiency->SetVfatEfficiency)
		    {
		      unsigned int absvfatid=d*100+m;		      
		      double eff=0.8;

		      if(((theVFatEfficiency)->EffiMap).find(absvfatid)!=((theVFatEfficiency)->EffiMap).end())
			eff=((theVFatEfficiency)->EffiMap)[absvfatid];
		      else
			std::cout<<"Efficiency not found for VFAT"<< m <<" in plane "<<d<<". Set to 0.8"<<std::endl;
		      
		      //if((absvfatid-1300<100)&&(absvfatid-1300>=0))
		      //std::cout<<"Plane"<<d<<"with absvfatid "<<absvfatid<<" efficiency to"<<eff<<std::endl;
		      /*
		      if((d<20)&&(d>=10))
			if(eff>=0.20)
			  eff=eff-0.20;
		      */

		      theVfat.Efficiency_=eff;		     

		      //The higher is the efficiency, the lower is the equivalent threshold.
		      //double thresholddependentEffi= bins_.at(d)*(eff); 
		      
		      double thresholddependentEffi= EqThrFromEffi(eff); 
		      // double thresholddependentEffi=bins_.at(d);
		      if(thresholddependentEffi<0)
			std::cout<<"Error in thresholddependentEffi Convertion: Thr:"<<thresholddependentEffi<<" Effi:"<<eff<<std::endl;
		      if(thresholddependentEffi<5){
			std::cout<<"Warning: low pad eq-threshold "<<thresholddependentEffi<<" for vfat:"<<absvfatid<<". Set Thr to 5"<<std::endl;
			thresholddependentEffi=5.0;
		      }

		       theVfat.SetThreshold(-1,thresholddependentEffi);
		       //std::cout<<"Test before: "<<theVfat.GetThreshold(13)<<std::endl;

		       int numdeadchannel=0;
		       if(((theVFatEfficiency)->VfatID_ToDeadChannelList).find(absvfatid)!=((theVFatEfficiency)->VfatID_ToDeadChannelList).end())
			 numdeadchannel=(theVFatEfficiency->VfatID_ToDeadChannelList)[absvfatid].size();

		      if(numdeadchannel>0)
			{
			  std::vector<unsigned int> listofDeadchannel=(theVFatEfficiency->VfatID_ToDeadChannelList)[absvfatid];			  		    theVfat.SetDeadChannels(listofDeadchannel);
			}
		      
		    }

		  
		  them1.insert(std::pair<unsigned int, T2DigiVfat>(m,theVfat));	     	      
		}    
	      PadVFats.insert(std::pair<int, std::map<unsigned int, T2DigiVfat> >(d,them1));	   
	    }

	  //  std::cout<<"VFMAP loaded with size "<<PadVFats.size()<<std::endl;
	}

    // Check the size of VFATS  Map Thresholds from DB acquisition and put default values------------------------------------

    if(UseCFGInfo==true)
      {

	for(unsigned int d=0;d<40;d++)
	  if(PadVFats.find(d)==PadVFats.end())
	    {
	      std::cout<<"Warning: Pad Thresholds in plane "<<d<<" not found. Put to default values."<<std::endl;
	      T2GeometryUtil conv;
	      T2GeometryUtil::T2DetInfo planeinformation;
	      planeinformation=conv.GetT2Info(d);   
	      std::map<unsigned int, T2DigiVfat> them1;

	      for(unsigned int m=2;m<15;m++)
		{
		   theVfat=T2DigiVfat(planeinformation.cmsswid,m,0);
		  //	std::auto_ptr<T2DigiVfat> theVfat(new T2DigiVfat(planeinformation.cmsswid,m,0));
		   theVfat.SetThreshold(-1,35);//DefaultThr 
		   theVfat.Efficiency_=0.8;
		  them1.insert(std::pair<unsigned int, T2DigiVfat>(m,theVfat));	   	   
		}
	      PadVFats.insert(std::pair<int, std::map<unsigned int, T2DigiVfat> >(d,them1));
	    }
	  else//Plane of vfat pad found
	    {
	      if((PadVFats[d]).size()!=13)
		{
		  std::cout<<"Warning: Only "<<(PadVFats[d]).size()<<"pad vfat found for threshold setting:" <<std::endl;
		  T2GeometryUtil conv;
		  T2GeometryUtil::T2DetInfo planeinformation;
		  planeinformation=conv.GetT2Info(d);   

		  for(unsigned int m=0;m<17;m++)
		    {
		      if((m>=2)&&(m<=14))
			{
			  if((PadVFats[d]).find(m)==(PadVFats[d]).end())
			    {

			      theVfat=T2DigiVfat(planeinformation.cmsswid,m,0);
			      theVfat.SetThreshold(-1,35);//DefaultThr
			      theVfat.Efficiency_=0.8;
			      //std::map<unsigned int, T2DigiVfat> them1;
			      //them1.insert(std::pair<unsigned int, T2DigiVfat>(m,(*theVfat)));
			      (PadVFats[d]).insert(std::pair<unsigned int, T2DigiVfat>(m,theVfat));
			      std::cout<<m<<" Put to default values."<<std::endl;
			    }

			}
		    }

		}
	    }
	//  std::cout<<"VFMAP completed, final size "<<PadVFats.size()<<std::endl;

       
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
            planenumb = key/100;
            vfnumb = key%100;
            if((vfnumb>=2)&&(vfnumb<=14)) //padVFAT
            {
              if(PadVFats.find(planenumb)!=PadVFats.end())
              {
                if(PadVFats[planenumb].find(vfnumb)!=PadVFats[planenumb].end())
                {
                  vector<unsigned int> vec;
                  for (set<unsigned char>::iterator it2 = it->second.maskedChannels.begin(); it2 != it->second.maskedChannels.end(); it2++)
                  {
                    vec.push_back(*it2);
                  }
                  PadVFats[planenumb][vfnumb].SetDeadChannels(vec);
                } else
                  std::cout<<"Warning: one vfat-strip not initialized"<<std::endl;
              } else
                std::cout<<"Warning: one plane not initialized in strip digitization"<<std::endl;
            }
          }
        } 
  
      }
    } //if useVfats END
  

   if(sigmaExtraNoise_.size()!=40) 
     {
       std::cout<<"T2Digi Warning: Found "<<sigmaExtraNoise_.size()<<" detector sigmaExtraNoise_ insted of 40. Missing will be put to default values"<<std::endl;
      if(sigmaExtraNoise_.size()<40)
	{
	  unsigned int miss=40-sigmaExtraNoise_.size();
	  for (unsigned int i=0;i<miss;i++)
	    {
	      sigmaExtraNoise_.push_back(0.);
	      //    std::cout<<"Det "<<sigmaExtraNoise_.size()<<"th at 0."<<std::endl;
	    }
	}
    }


   if(capaNoiseFactorPad_.size()!=40) 
     {
       std::cout<<"T2Digi Warning: Found "<<capaNoiseFactorPad_.size()<<" detector capaNoiseFactorPad_ insted of 40. Missing will be put to default values"<<std::endl;
      if(capaNoiseFactorPad_.size()<40)
	{
	  unsigned int miss=40-capaNoiseFactorPad_.size();
	  for (unsigned int i=0;i<miss;i++)
	    {
	     capaNoiseFactorPad_.push_back(0.);
	     // std::cout<<"Det "<<capaNoiseFactorPad_.size()<<"th at 1."<<std::endl;
	    }
	}
    }



   if(simpleThreshold_.size()!=40) 
     {
       std::cout<<"T2Digi Warning: Found "<<simpleThreshold_.size()<<" detector simpleThreshold_ insted of 40. Missing will be put to default values"<<std::endl;
      if(simpleThreshold_.size()<40)
	{
	  unsigned int miss=40-simpleThreshold_.size();
	  for (unsigned int i=0;i<miss;i++)
	    {
	      simpleThreshold_.push_back(600);
	      //    std::cout<<"Det "<<simpleThreshold_.size()<<"th at 600"<<std::endl;
	    }
	}
    }






  
} // T2PadElectronicsSim

/**
 *
 */

double T2PadElectronicsSim::EqThrFromEffi(double effi_measured){
  
  /*
    Gaussian fit obtained with single muon at SimpleThr=400
    StripWidth=0.14    diffCoeff=0.4    gain=15000    capaNoiseFactorStrip=1    simpleThreshold=400    sigmaExtraNoise=0
    Simu have been used with fully alive and 0 corrupted vfats, no dead sector.
    Gaussian Parameter: C=1.016 mu=-62.82  sigma=284.3 where  f(x) = p0*exp(-0.5*((x-p1)/p2)^2))  
    
    Old Test
    StripWidth=0.14    diffCoeff=0.4    gain=15000    capaNoiseFactorStrip=1    simpleThreshold=400    sigmaExtraNoise=0
    Extracted Param:
    double eqThr=0.;
    double mu=(-1.0)*62.82;
    double C=1.016;
    double sigma=284.3;
  */
  /*
  //StripWidth=0.05    diffCoeff=0.24    gain=30000    capaNoiseFactorStrip=0.1    simpleThreshold=400    sigmaExtraNoise=0
  //Extracted Param:
  double eqThr=0.;
  double mu=(-1.0)*50.28;
  double C=1.002;
  double sigma=334.7;
  */
  /*
  //StripWidth=0.05    diffCoeff=0.21    gain=25000    capaNoiseFactorStrip=0.1    simpleThreshold=400    sigmaExtraNoise=0
  //Extracted Param:
  double eqThr=0.;
  double mu=(-1.0)*47.59;
  double C=1.004;
  double sigma=296.4;
  */
  
  /*
//StripWidth=0.05    diffCoeff=0.27    gain=25000    capaNoiseFactorStrip=0.1    simpleThreshold=400P/550S    sigmaExtraNoise=0
  //Extracted Param:
  double eqThr=0.;
  double mu=(-1.0)*206.6;
  double C=0.997;
  double sigma=948.9;
  */

  double logar=0.;
  //Extracted Param: Curve G
  double eqThr=0.;
  double mu=(-1.0)*148.1;
  double C=1.001;
  double sigma=693.3;

  if(effi_measured>0.1)
    {
      if((effi_measured/C)<1.)//Normal Case, logar<0.
	{
	  logar=log(effi_measured/C);
	  eqThr=mu+sigma*sqrt(-2.0*logar);
	  if(eqThr<0)//Can happen if mu<0. Use a minimum value.
	    eqThr=5;
	}
      else//Very High Effi case: put lower Thr limit.
	{
	  eqThr=5;//Was 10. il 20/7/2012
	}
    }
  else //Not very predictable at low effi, give back 4 sigma.
    {
      eqThr=3*sigma;
      if(effi_measured<0.05)
	eqThr=300000.0*sigma;
    }

  


  return eqThr;
  
  
}



void T2PadElectronicsSim::LoadDeadSector(std::string inputFileNameDeadSect,std::vector<unsigned int> &VectDeadSect_Plane,std::vector<unsigned int> &VectDeadSect_Sector)
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
		VectDeadSect_Plane.push_back(((unsigned int)detid));  
		VectDeadSect_Sector.push_back(((unsigned int)sectordead));  
	      }
	    
	  }
    }
  else
    std::cout<<"Dead Sect file not opened"<<std::endl;
}



unsigned int T2PadElectronicsSim::RawtoSymb(uint32_t thedet)
{
  T2DetId converter;
  unsigned int pl=converter.plane(thedet);
  unsigned int pls=converter.planeSide(thedet);
  unsigned int ht=converter.halfTelescope(thedet);
  unsigned int arm=converter.arm(thedet);
  unsigned int symbolic=pl*2+pls+ht*10+20*arm;
  return symbolic;
}


bool T2PadElectronicsSim::IsPadInDeadSector(unsigned int RVal,unsigned int symbdetid){
  bool toreturn=false;
  //R is from 0 to 24

  /*
    rMinPlane =  42.46; // all dimensions should be in cm (CMS convention)
    // for pads 42.55 and for strips 42.46 TO CHECK!!!!
    rMaxPlane = 144.58; // for pads 144.58 and for strips 144.54 TO CHECK!!!!

    Pad max R{44.677,47.019,49.483,52.076,54.806,57.677,60.700,63.880,
    67.227,70.726,74.455,78.356,82.461,86.781,91.327,96.111,
    101.145,106.443,112.011,117.878,124.050,130.547,137.384,144.580}
if(RVal<50.) Rsector=0; ->2/3
 else if(RVal<92.5) ->14 //91.327
        Rsector=1;
      else if(RVal<121.5)->19/20
        Rsector=2;
      else Rsector=3;

  */

  unsigned int Rsector=0;
  if(RVal<3)
    {
      Rsector=0;
    }
  else
    {
      if(RVal<15)
	{
	  Rsector=1;
	}
      else
	{
	  if(RVal<20)
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


void T2PadElectronicsSim::readFile(std::string fileName, double* dataFile) {
  
  std::string path(getenv( "CMSSW_SEARCH_PATH" ));

  FileInPath file(path, fileName);

  if (file() == 0 ) {
    
    edm::LogError("T2PadElectronicsSim") << "File with Capacitancenoise for VFat not found: " 
					 << file.name();
    throw cms::Exception(" Capacitance noise file not found!");
  } 

  LogDebug("T2PadElectronicsSim") <<"Reading in " << file.name();
  edm::LogInfo("T2PadElectronicsSim") << "-> Reading in " << file.name();
 
  std::ifstream & fin = *file();
  
  LogDebug("T2PadElectronicsSim") << "T2PadElectronicsSim: opening file " << file.name();

  if (fin == 0) {
    std::string errorMessage =  "Cannot open file " +  file.name();
    edm::LogError("T2PadElectronicsSim") << errorMessage;
    throw cms::Exception(errorMessage);
  }

  double var1;
  double var2;
  
  int channel = 0;
  
  while (fin.good() && channel<NUMBER_OF_PCHANNELS) {
    
    fin >> var1 >> var2;
    
    dataFile[channel] = var2;
    
    channel++;
  }

  if (channel<NUMBER_OF_PCHANNELS) {
    
    throw cms::Exception("Data file error") << "To little data in file " << file.name();
  }
      
  fin.close();
  
} // readFile

/**
 *
 */

void T2PadElectronicsSim::simulate(std::map<int, int*> & chargeMap) {
  
  //  std::cout<<"PadElectr. simulate Beg."<<std::endl;
  std::map<int, int*>::iterator i;
  for (i=chargeMap.begin(); i != chargeMap.end(); ++i) {

   unsigned int symbdetid= RawtoSymb((*i).first);
   int additionalnoise=0.;
   int capacitancenoise=0.;
    
    for (int sector=0; sector<13; sector++) {
      
      for(int k=0; k<NUMBER_OF_PCHANNELS; k++) {
	
	// Warning:
	// noise: Gauss with Sigma from capacitance noise. Mean from energy deposition
	//il file di dati va da aree + grandi a + piccole per Erik mentre per Eraldo da + piccole a + grandi (la formula di Erik funziona con la convenzione di Eraldo)
	//(*i).second[k+sector*NUMBER_OF_PCHANNELS] += capaNoiseFactorPad_.at(symbdetid)*(static_cast<int>(ceil(RandGaussQ::shoot( 0, capaNoiseSigma[0][k]))))
	
	capacitancenoise= static_cast<int>(capaNoiseFactorPad_.at(symbdetid)*(static_cast<int>(ceil(CLHEP::RandGaussQ::shoot( 0, simpleThreshold_.at(symbdetid)*capaNoiseSigma[0][k])))));
	
	additionalnoise= static_cast<int>(capaNoiseFactorPad_.at(symbdetid)*(static_cast<int>(ceil(CLHEP::RandGaussQ::shoot( 0, simpleThreshold_.at(symbdetid)*sigmaExtraNoise_.at(symbdetid) )))));
	
	
	(*i).second[k+sector*NUMBER_OF_PCHANNELS] =(*i).second[k+sector*NUMBER_OF_PCHANNELS]+capacitancenoise+additionalnoise ;

	//std::cout<<"capac-addit-PadCharge=   "<<capacitancenoise<<" "<<additionalnoise<<" "<<(*i).second[k+sector*NUMBER_OF_PCHANNELS]<<std::endl;

	//(*i).second[k+sector*NUMBER_OF_PCHANNELS] += capaNoiseFactorPad_.at(symbdetid)*(static_cast<int>(ceil(RandGaussQ::shoot( 0, simpleThreshold_.at(symbdetid)*capaNoiseSigma[0][k]))))
	  	  // additional noise to fill up to meassured noise
	//+capaNoiseFactorPad_.at(symbdetid)*static_cast<int>(ceil(RandGaussQ::shoot( 0, simpleThreshold_.at(symbdetid)*sigmaExtraNoise_.at(symbdetid) )));

	/*
	if(symbdetid==0)
	  {
	    if((k==22)&&(sector==0))
	      {
		std::cout<<"k=22  "<<(*i).second[k+sector*NUMBER_OF_PCHANNELS]<<std::endl;
		std::cout<<"capa "<<capaNoiseSigma[0][k]<<std::endl;
	      }
	    //  if((k==2)&&(sector==0))
	    // std::cout<<"k=2  "<<(*i).second[k+sector*NUMBER_OF_PCHANNELS]<<std::endl; 
	    
	  }
	*/  

      }
    }
  }
  // std::cout<<"PadElectr. simulate End."<<std::endl;
} // simulate

/**
 *
 */

void T2PadElectronicsSim::fillDigis(T2PadDigiCollection & padDigis,
				    std::map<int, int*> & chargeMap) {
  
  std::map<int, int*>::iterator i;
  
  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinformation;
  unsigned int VfatID;
  int detectorID;
  unsigned int vfatchannel;
  int vfatchannelThr;
  unsigned int stepincolumn;
  //int VfatEfficiency;

  


  for (i=chargeMap.begin(); i != chargeMap.end(); ++i) {
    for (int c=0; c<65; c++) {
      for(int r=0; r<24; r++) {
	
	VfatID=2+c/5;                    //each vfat has 5 col //2..14
	stepincolumn=c-(VfatID-2)*5;     //c-(c/5)*5 
	vfatchannel=stepincolumn*25+r;  
	if(vfatchannel>128)
	  std::cout<<"Error --- vfat data pad channel badly calculated: "<<vfatchannel<<std::endl;

	detectorID=(*i).first;
                  
	//Check if it is correct!!!!!!!! //Read data channel:4-123 only for pad
	planeinformation=conv.GetT2Info(detectorID);
	unsigned int symbdetid= planeinformation.symb;

	if(UseVFATs)
	  {
	    vfatchannelThr=PadVFats[symbdetid][VfatID].GetThreshold(vfatchannel);
	    if(vfatchannelThr<0)
	      std::cout<<"Error In T2PadElectronicSim: vfat channel thr:"<<vfatchannelThr<<" Symb:"<<symbdetid<<" VfatID:"<<VfatID<<" Channel:"<<vfatchannel<<std::endl;
	   
	    // std::cout<<"Take Thr"<<vfatchannelThr<<std::endl;
	    // VfatEfficiency=1.0;
	  }
	else
	  {
	    vfatchannelThr=bins_.at(symbdetid);
	    //VfatEfficiency=1.0;
	  }

	

	int charge = (*i).second[r + c*24];       
	
	//double electronicEfficiency=1.0; 
	double randomValue=1.;
	double corruptProb=0.;

	
	if(UseVFATs)
	  if(theVFatEfficiency->SetVfatEfficiency)
	    {	      
	      unsigned int Thisabsvfatid=symbdetid*100+VfatID;

	      if(((theVFatEfficiency)->CorruptMap).find(Thisabsvfatid)!=((theVFatEfficiency)->CorruptMap).end())
		corruptProb=((theVFatEfficiency)->CorruptMap)[Thisabsvfatid];
	      
	      //electronicEfficiency= PadVFats[symbdetid][VfatID].Efficiency_; 
	      randomValue=(double)(CLHEP::RandFlat::shoot(0.,1.));
	    }
	//std::cout<<"Channel "<<vfatchannel<<" of Vfat n?? "<<VfatID<<"In plane "<<planeinformation.symb<<" has threshold "<<vfatchannelThr<<std::endl;
	//unsigned int symbdetid= RawtoSymb((*i).first);
	//	if(charge > threshold[c][r]) { 
	//if(r==22)
	//  if(c==0)
	//  std::cout<<"Final charge and threshold: "<<charge<<" "<<(bins_.at(symbdetid)* simpleThreshold_.at(symbdetid))<<std::endl;
	

	if(PadVFats[symbdetid][VfatID].IsChannelDead(vfatchannel)==false)
	  {
	    bool padIndeadsector=IsPadInDeadSector((unsigned int) r,symbdetid);
	   
	    //if((charge > vfatchannelThr * simpleThreshold_.at(symbdetid))&&(electronicEfficiency>randomValue))
	    //	    std::cout<<"Symb plane:"<<symbdetid<<" Thr:"<<vfatchannelThr<<" ThrConv:"<<vfatchannelThr * simpleThreshold_.at(symbdetid)<<" Ch:"<<charge<<std::endl;
	    //Threshold = number of bins * around 1200(e-)
	    

	    if(charge > vfatchannelThr * simpleThreshold_.at(symbdetid))
	      if(randomValue>corruptProb)//First electronic gate which take into account corrupt-prob and dead vfats. 
		if(padIndeadsector==false)
		  { 
		    //std::cout<<"symbdetid-charge-thr "<<symbdetid<<" "<<charge<<" "<<(bins_.at(symbdetid)* simpleThreshold_.at(symbdetid))<<std::endl;
		    T2DetId planeId((*i).first);
		    T2PadDigi sDigi(0, r, c, charge);
                    
		    padDigis.insertDigi(planeId, sDigi);	  	 
		    PadVFats[symbdetid][VfatID].SetChannel(vfatchannel,1);
		    //		    std::cout<<"One shitty pad ON at r "<<r<<" c "<<c<<" chrg "<<charge;
		  }
	  }
	/*
	else{
	  //No further action taken at the moment in order to correct the cluster.
	  //std::cout<<"Warn: vfat pad-channel simulated dead"<<std::endl;
	  T2PadDigi sDigi(0, r, c, charge);
	  T2DetId planeId((*i).first);
	  deadpadDigisActivated.insertDigi(planeId, sDigi);	  
	}
	*/
      }
    }
  }    


  //Now You need to correct for the Dead Channel.


  
} // fillDigis

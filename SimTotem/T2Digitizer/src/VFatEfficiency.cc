/**
 * Class VFatEfficiency
 *
 * Author: Mirko Berretti / University of Siena 
 * Email:  mirko.berretti@gmail.com
 * Date:   2007-12-08
 */

#include "FWCore/Framework/interface/ESHandle.h"
#include <iostream>
#include "SimTotem/T2Digitizer/interface/VFatEfficiency.h"
#include <fstream>
#include <string>
#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"


VFatEfficiency::VFatEfficiency(bool SetVfatEfficiency_,std::string inputFileNameEffi_,std::string inputFileNameCorrupt_,std::string inputFileNameDeadChannels_,std::string inputFileNameNoisyChannels_) 
{
  
  
  inputFileNameEffi= inputFileNameEffi_;//conf.getUntrackedParameter<std::string>("inputFileNameEffi");
  inputFileNameCorrupt= inputFileNameCorrupt_;
  inputFileNameDeadChannels = inputFileNameDeadChannels_;
  inputFileNameNoisyChannels = inputFileNameNoisyChannels_;
  

  SetVfatEfficiency= SetVfatEfficiency_;//conf.getUntrackedParameter<bool>("SetVfatEfficiency");
  bool verbosity= false;//conf.getUntrackedParameter<bool>("verbosity");

  unsigned int expectedrawsize=5;
  unsigned int expectedrawsizeCorr=3;
  T2GeometryUtil conv;
 

  // Expected input Efficiency file format:
  // DetId: 0  VFat_iid: 0   Eff: 1        pm 0 ( stat: 1 )

  // Expected input Corruption file format:
  // DetId: 0  VFat_iid: 0   CorruptProb: 0

  std::vector<double> onevfateffi;
  for(unsigned int u=0;u<5;u++)
    {
      onevfateffi.push_back(0.);
    }

  std::vector<double> onevfatcorrupt;
  for(unsigned int u=0;u<3;u++)
    {
      onevfatcorrupt.push_back(0.);
    }

  
  //Initialization of the maps
  if(SetVfatEfficiency)
    {
      for (unsigned int j=0;j<40;j++)
	for (unsigned int vf=0;vf<17;vf++)
	{
	  int detid=j;
	  unsigned int VfatAbsID=100*detid+vf; //DQM convention

	  onevfateffi.at(0)=detid;
	  onevfateffi.at(1)=vf;// Vfat_iid
	  onevfateffi.at(2)=0.8;//Eff, average initialization
	  onevfateffi.at(3)=0.;//pm
	  onevfateffi.at(4)=0.;//stat

	  onevfatcorrupt.at(0)=detid;
	  onevfatcorrupt.at(1)=vf;// Vfat_iid
	  onevfatcorrupt.at(2)=0.;  //CorruptProb

	  std::vector<unsigned int> EmptydeadChannelvect;	 
	  EmptydeadChannelvect.clear();
	  VfatID_ToDeadChannelList.insert(pair<unsigned int,std::vector<unsigned int> >(VfatAbsID,EmptydeadChannelvect));
	  
	  effvect.push_back(onevfateffi);
	  corruptvect.push_back(onevfatcorrupt);
	  
	}
      
       if(verbosity)
	 std::cout<<"Vfat Efficiency-Corruption Init done "<<effvect.size()<<std::endl;
       
    }

  
  //-----------------------------------------------------------------------
  //Read Noise channel list from a file         >> And Put IT AS DEAD <<
  
  if((SetVfatEfficiency)&&(inputFileNameNoisyChannels!=""))
    {
      char *cmsswPath = getenv("CMSSW_BASE");
      inputFileNameNoisyChannels= string(cmsswPath) + string("/src/") + inputFileNameNoisyChannels;
      std::string line;
      ifstream myfile(inputFileNameNoisyChannels.c_str()); 

      int Vfat_Absiid=-1; int DetId=-1; int deadChannel=-1;
      if (myfile.is_open())
	{
	  while (! myfile.eof() )
	    {
	      getline (myfile,line);
	      //cout << line << endl;
	      bool foundDetId=false;
	      bool foundVfat_Absiid=false;	      
	      bool founddeadChannel=false;
	     
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
		    //expected format: DetId:<spaces>i<spaces>Vfat_iid:<spaces>a<spaces> .....
		    
		    oneword=tokens.at(m);
		    
		    if(verbosity)
		      std::cout<<oneword<<std::endl;

		    if(oneword=="Plane:")
		      {
			DetId = atoi(tokens.at(m+1).c_str());
			if((DetId>=0)&&(DetId<40))
			  foundDetId=true;
		      }

		    if(oneword=="AbsVfatID:")
		      {		
			Vfat_Absiid= boost::lexical_cast<double>(tokens.at(m+1));			
			foundVfat_Absiid=true;
		      }

		    if(oneword=="NoisyChannel:"){
		      deadChannel= boost::lexical_cast<double>(tokens.at(m+1));
		      founddeadChannel=true;
		    }
		    
		  }

		if((foundDetId)&&(foundVfat_Absiid)&&(founddeadChannel))		  
		  {
		    //--------------------------------------
		    //Fill the maps AbsId<->Vector DeadChannels
		    VfatID_ToDeadChannelList[Vfat_Absiid].push_back(deadChannel);
		   		    
		  }			  
		else
                 if(tokens.size()>0)
		   std::cout<<"Vfat Noisy Channels not read as expected: "<<DetId<<" "<<Vfat_Absiid<<" "<<deadChannel<<" |"<<(foundDetId)<<(foundVfat_Absiid)<<(founddeadChannel)<<std::endl;

	    }	  
	  myfile.close();
	}
      else
	{ 
	  cout << "Warning : VfatEfficiency Unable to open file with name:"; 
	  std::cout<<" "<<inputFileNameNoisyChannels.c_str()<<std::endl;
	}

    }

  //--------------------------------------
  //Read Dead channel list from a file 
  if((SetVfatEfficiency)&&(inputFileNameDeadChannels!=""))
    {
      int totdeads=0;
       char *cmsswPath = getenv("CMSSW_BASE");
       inputFileNameDeadChannels= string(cmsswPath) + string("/src/") + inputFileNameDeadChannels;
       std::string line;
       ifstream myfile (inputFileNameDeadChannels.c_str());
       int Vfat_Absiid=-1; int DetId=-1; int deadChannel=-1;
       //Format Expected 
       // outputFileDeadChannels_AliveVFAT<<"  Plane: "<<(absvfatId/100)<<"  AbsVfatID: "<<absvfatId<<"  DeadChannel: "<<iBin<<"\n";
        if (myfile.is_open())
	  {
	  while (! myfile.eof() )
	    {
	      getline (myfile,line);
	      //cout << line << endl;
	      bool foundDetId=false;
	      bool foundVfat_Absiid=false;	      
	      bool founddeadChannel=false;
	     
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
		    //expected format: DetId:<spaces>i<spaces>Vfat_iid:<spaces>a<spaces> .....
		    
		    oneword=tokens.at(m);
		    
		    if(verbosity)
		      std::cout<<oneword<<std::endl;

		    if(oneword=="Plane:")
		      {
			DetId = atoi(tokens.at(m+1).c_str());
			if((DetId>=0)&&(DetId<40))
			  foundDetId=true;
		      }

		    if(oneword=="AbsVfatID:")
		      {		
			Vfat_Absiid= boost::lexical_cast<double>(tokens.at(m+1));			
			foundVfat_Absiid=true;
		      }

		    if(oneword=="DeadChannel:"){
		      deadChannel= boost::lexical_cast<double>(tokens.at(m+1));
		      founddeadChannel=true;
		    }
		    
		  }

		if((foundDetId)&&(foundVfat_Absiid)&&(founddeadChannel))		  
		  {
		    //--------------------------------------
		    //Fill the maps AbsId<->Vector DeadChannels
		    VfatID_ToDeadChannelList[Vfat_Absiid].push_back(deadChannel);
		    totdeads++;		    
		  }			  
		else
		  std::cout<<"Vfat Dead Channels not read as expected: "<<DetId<<" "<<Vfat_Absiid<<" "<<deadChannel<<" |"<<(foundDetId)<<(foundVfat_Absiid)<<(founddeadChannel)<<std::endl;

	    }	  
	  myfile.close();
	  }
	else
	  { 
	    cout << "Warning : VfatEfficiency Unable to open file with name:"; 
	    std::cout<<" "<<inputFileNameDeadChannels.c_str()<<std::endl;
	  }

	//std::cout<<"Total Dead Channel Loaded: "<<totdeads<<std::endl;
  
    }

  





 
  

  //--------------------------------------
  //Read Efficiency from a file (Efficiency is normalized per Dead Vfat, Gem sector and dead channel (in future))
  if((SetVfatEfficiency)&&(inputFileNameEffi!=""))
    {

      char *cmsswPath = getenv("CMSSW_BASE");
      inputFileNameEffi= string(cmsswPath) + string("/src/") + inputFileNameEffi;

      std::string line;
      ifstream myfile (inputFileNameEffi.c_str());
      int DetId=-1;
      double Vfat_iid=0.,Eff=0.;
      double pm=0., stat=0.;     
      
      
      if (myfile.is_open())
	{
	  while (! myfile.eof() )
	    {
	      getline (myfile,line);
	      //cout << line << endl;
	      bool foundDetId=false;
	      bool foundVfat_iid=false;
	      bool foundEff=false;
	     
	      
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
		    //expected format: DetId:<spaces>i<spaces>Vfat_iid:<spaces>a<spaces> .....
		    
		    oneword=tokens.at(m);
		    
		    if(verbosity)
		      std::cout<<oneword<<std::endl;

		    if(oneword=="DetId:")
		      {
			DetId = atoi(tokens.at(m+1).c_str());
			if((DetId>=0)&&(DetId<40))
			  foundDetId=true;
			else
			  std::cout<<"Vfat efficiency Wrong DET-ID:"<<DetId<<std::endl;
		      }

		     
		    if(oneword=="VFat_iid:")
		      {
		
			Vfat_iid= boost::lexical_cast<double>(tokens.at(m+1));
			if((Vfat_iid>=0)&&(Vfat_iid<=16))
			  foundVfat_iid=true;
			else
			  std::cout<<"Vfat efficiency Wrong Vfat_iid:"<<Vfat_iid<<std::endl;

			if(verbosity)
			  std::cout<<"VFat_iid= "<<Vfat_iid<<std::endl;
			
		      }


		    if(oneword=="Eff:")
		      {	
		        Eff= boost::lexical_cast<double>(tokens.at(m+1));
			
			if(verbosity)
			  std::cout<<"Eff= "<<Eff<<std::endl;

			if(Eff<=1.0)
			  foundEff=true;
			else
			  std::cout<<"Vfat efficiency Wrong Efficiency:"<<Eff<<std::endl;
			
		      }

		    if(oneword=="pm:")
		      {
			pm=  boost::lexical_cast<double>(tokens.at(m+1));
		

		
			if(verbosity)
			  std::cout<<"pm= "<<pm<<std::endl;
		      }

		    if(oneword=="stat:")
		      {
			stat= boost::lexical_cast<double>(tokens.at(m+1));	    			
			if(verbosity)
			  std::cout<<"stat= "<<stat<<std::endl;
		      }		  
		  }
		
		

		if((foundDetId)&&(foundVfat_iid)&&(foundEff)/*&&(foundstat)&&(foundpm)*/)
		  {
		    
		    onevfateffi.at(0)=DetId;
		    onevfateffi.at(1)=Vfat_iid;
		    onevfateffi.at(2)=Eff;
		    onevfateffi.at(3)=pm;
		    onevfateffi.at(4)=stat;
		    
		    unsigned int AbsposInTempVector=DetId*17+Vfat_iid;
		    //  std::cout<<"vfat status updated"<<std::endl;
		    effvect.at(AbsposInTempVector)=onevfateffi;
		  }			  
		else
		  std::cout<<"Vfat Effi not read as expected: "<<DetId<<" "<<Vfat_iid<<" "<<Eff<<" |"<<(foundDetId)<<(foundVfat_iid)<<(foundEff)<<std::endl;

	    }
	  myfile.close();
	}
      else
	{ 
		cout << "Warning : VfatEfficiency Unable to open file with name:"; 
		std::cout<<" "<<inputFileNameEffi.c_str()<<std::endl;
	}
  
      //--------------------------------------
      //Create maps AbsId<->VfatEffi

      //If you don't read from DB


      T2DetId detconverter;
      unsigned int VfatAbsID;
      
      for(unsigned int l=0;l<effvect.size();l++)
	if(effvect.at(l).size()==expectedrawsize)
	  {
	    if(verbosity)
	      std::cout<<"Det "<< effvect.at(l).at(0) <<" taken for Set Effi of Vfat"<<" :"<<effvect.at(l).at(1)<<std::endl;
	    
	    VfatAbsID=100*effvect.at(l).at(0)+effvect.at(l).at(1); //DQM convention
	    EffiMap.insert(pair<unsigned int,double>(VfatAbsID,effvect.at(l).at(2)));
	    
	  } 
	else
	  std::cout<<"Problem: Vfat expected size "<<expectedrawsize<<" instead of "<<effvect.at(l).size()<<std::endl;
	  
	  
    }




  


  //--------------------------------------
  //Read Corruption from a file 

  if((SetVfatEfficiency)&&(inputFileNameCorrupt!=""))
    {

      char *cmsswPath = getenv("CMSSW_BASE");
      inputFileNameCorrupt= string(cmsswPath) + string("/src/") + inputFileNameCorrupt;

      std::string line;
      ifstream myfile (inputFileNameCorrupt.c_str());
      int DetId=-1;
      int Vfat_iid=0;
      double CorrProb=0.;
      
      
      if (myfile.is_open())
	{
	  while (! myfile.eof() )
	    {
	      getline (myfile,line);
	      //cout << line << endl;
	      bool foundDetId=false;
	      bool foundVfat_iid=false;
	      bool foundCorrProb=false;
	      
	      
	      istringstream in(line);
	      
	      std::vector<std::string> tokens; // Create vector to hold our words
	      std::string oneword;
	      //scan the line, look at your Tags
		while(in>>oneword)
		  {	
		    tokens.push_back(oneword);
		  }
		
		
		DetId=-1;
		Vfat_iid=-1;
		CorrProb=-1.0;

		for(unsigned int m=0;m<tokens.size();m++)
		  {
		    //"parola" each time is overwritten
		    //expected format: DetId:<spaces>i<spaces>Vfat_iid:<spaces>a<spaces> .....
		    
		    oneword=tokens.at(m);
		    
		    if(verbosity)
		      std::cout<<oneword<<std::endl;

		    if(oneword=="DetId:")
		      {
			DetId = atoi(tokens.at(m+1).c_str());
			if((DetId>=0)&&(DetId<40))
			  foundDetId=true;
		      }

		     
		    if(oneword=="VFat_iid:")
		      {
		
			Vfat_iid=  atoi(tokens.at(m+1).c_str());//boost::lexical_cast<double>(tokens.at(m+1));
			if(verbosity)
			  std::cout<<"VFat_iid= "<<Vfat_iid<<std::endl;
			foundVfat_iid=true;
		      }


		    if(oneword=="CorruptProb:")
		      {	
		       CorrProb= boost::lexical_cast<double>(tokens.at(m+1));			
			if(verbosity)
			  std::cout<<"CorruptProb= "<<CorrProb<<std::endl;

			foundCorrProb=true;
		      }

		    	  
		  }
		
		

		if((foundDetId)&&(foundVfat_iid)&&(foundCorrProb)/*&&(foundstat)&&(foundpm)*/)
		  {		    		    
		     
		    int AbsposInTempVector=DetId*17+Vfat_iid;
		    onevfatcorrupt.at(0)=(double)DetId;
		    onevfatcorrupt.at(1)=(double)Vfat_iid;
		    onevfatcorrupt.at(2)=CorrProb;		    
		    //  std::cout<<"vfat status updated"<<std::endl;
		    if(AbsposInTempVector>=0)
		      corruptvect.at(AbsposInTempVector)=onevfatcorrupt;
		    else
		      std::cout<<"Error: Abspos "<<DetId<<"|"<<Vfat_iid<<" not saved"<<std::endl;
		  }			  
		else
		  std::cout<<"Corruption not read as expected: "<<DetId<<" "<<Vfat_iid<<" "<<CorrProb<<" |"<<(foundDetId)<<(foundVfat_iid)<<(foundCorrProb)<<std::endl;

	    }
	  myfile.close();
	}
      else
	{ 
		cout << "Warning : VfatEfficiency Unable to open file with name:"; 
		std::cout<<" "<<inputFileNameCorrupt.c_str()<<std::endl;
	}
  
      //--------------------------------------
      //Create maps AbsId<->VfatCorruption

  	
      T2DetId detconverter;
      unsigned int VfatAbsID;
      
      for(unsigned int l=0;l<corruptvect.size();l++)
	if(corruptvect.at(l).size()==expectedrawsizeCorr)
	  {
	    if(verbosity)
	      std::cout<<"Det "<< corruptvect.at(l).at(0) <<" taken for Set Corruption of Vfat"<<" :"<<corruptvect.at(l).at(1)<<std::endl;
	    
	    VfatAbsID=100*corruptvect.at(l).at(0)+corruptvect.at(l).at(1); //DQM convention
	    CorruptMap.insert(pair<unsigned int,double>(VfatAbsID,corruptvect.at(l).at(2)));
	    
	  } 
	else
	  std::cout<<"Problem Corrupt: expected size "<<expectedrawsizeCorr<<" instead of "<<corruptvect.at(l).size()<<std::endl;
	  
	  
    }
 

  //std::cout<<"Vfat Configs loaded"<<std::endl;

  
} // END VfatEfficiency::VfatEfficiency(....)



VFatEfficiency::~VFatEfficiency() 
{

 
  //std::cout<<EffiMap<<" "<<CorruptMap<<std::endl;
  /*
  std::cout<<"EFFICIENCY MAP READ:"<<std::endl;
  std::map<unsigned int,double>::iterator effiIter;
  for(effiIter=EffiMap.begin();effiIter!=EffiMap.end();effiIter++){
    std::cout<<"Plane:"<<((effiIter->first)/100)<<"  VFatId:"<<((effiIter->first)%100)<<"  Effi:"<<effiIter->second<<std::endl;   
  }
  */
  //TFile *f = TFile::Open("./test.root", "recreate");
  //if( !f || !f->IsWritable() ){
  // std::cout << "Output file not opened correctly !!" << std::endl;
  // }
  //TestDx->Write("");

}






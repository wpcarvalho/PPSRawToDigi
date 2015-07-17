/**
 * Class T2RoadProducer.
 *
 * Author: Mirko Berretti / University of Siena 
 * Email:  mirko.berretti@gmail.com
 * Date:   2010- 05 Dec.
 */

#include <new>                  // std::bad_alloc
#include <algorithm>            // for std::swap
#include <typeinfo>             // for std::bad_cast
#include <memory>               // std::allocator



#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoTotemT1T2/T2RoadPadFinder/interface/Tube.h"
#include "RecoTotemT1T2/T2RoadPadFinder/interface/T2RoadPadFinder.h"
#include "DataFormats/T2Road/interface/T2Road.h"
#include "DataFormats/T2Road/interface/T2RoadCollection.h"
#include "DataFormats/T2Hit/interface/T2Hit.h"
#include "DataFormats/T2Hit/interface/T2HitCollection.h"
#include "DataFormats/T2Hit/interface/T2Hit_to_Track_Map.h"
#include "DataFormats/T2Cluster/interface/T2Cluster.h"
#include "DataFormats/T2Cluster/interface/T2PadClusterCollection.h"
#include "DataFormats/T2Cluster/interface/cluster_entry.h"
#include <iostream>
#include <math.h>



T2RoadPadFinder::T2RoadPadFinder(const edm::ParameterSet &config){

  T2RoadCollProdName = config.getParameter<string>("T2RoadCollProdName");//T2RoadColl
  HitToRoadAssociatorName = config.getParameter<string>("HitToRoadAssociatorName");
  QuartedSelected=config.getParameter<std::vector<int> >("QuartedSelected");
  VplaneToExclude=config.getParameter<std::vector<int> >("VplaneToExclude");

  TolleranceThetaX= config.getParameter<double>("TolleranceThetaX");
  TolleranceThetaY= config.getParameter<double>("TolleranceThetaY");
  TolleranceDX= config.getParameter<double>("TolleranceDX");
  TolleranceDY= config.getParameter<double>("TolleranceDY");
  InefficiencyMaxJump= config.getParameter<int>("InefficiencyMaxJump");
  AllowsPadReAssociation= config.getParameter<bool>("AllowsPadReAssociation");
  AllowsConcurrentBranches= config.getParameter<bool>("AllowsConcurrentBranches");
  Nmin_padsFinal= config.getParameter<int>("Nmin_padsFinal");
  verbosity= config.getParameter<int>("verbosity");
  TwoPointsTubesAngularCollinearity=config.getParameter<double>("TwoPointsTubesAngularCollinearity");
  HitLabel= config.getParameter<string>("HitLabel");
  CluLabel= config.getParameter<string>("CluLabel");

  MinimumNumCl1Hit= config.getParameter<int>("MinimumNumCl1Hit");
  chi2XYProb_Thr= config.getParameter<double>("chi2XYProb_Thr");
  
  ResolveOverlapDoubleCount= config.getParameter<bool>("ResolveOverlapDoubleCount");
  OverlapDoubleCountDR= config.getParameter<double>("OverlapDoubleCountDR");
  OverlapDoubleCountDPhi= config.getParameter<double>("OverlapDoubleCountDPhi");
  OverlapDoubleCountDTheta= config.getParameter<double>("OverlapDoubleCountDTheta");   
  MinCluSize_considered_asBlobs = config.getParameter<int>("MinCluSize_considered_asBlobs"); 
    
  BiggestTubeAngleConsideredINPUT= config.getParameter<double>("BiggestTubeAngleConsidered");  
  NumSigma= config.getParameter<double>("NumSigma");

  useStraightPadTowers= config.getParameter<bool>("useStraightPadTowers");
  NumPadCluOccupancyAlert= config.getParameter<double>("NumPadCluOccupancyAlert");

  produces<T2RoadCollection>(T2RoadCollProdName);  
  produces<T2Hit_to_Track_Map>(HitToRoadAssociatorName);
}


T2RoadPadFinder::~T2RoadPadFinder(){
}


template <typename T> void T2RoadPadFinder::FreeVectorMemory(T & t){
    T tmp;
    t.swap( tmp );
}
 

void  T2RoadPadFinder::InitializeMapClusterAssociators(const edm::Event& iEvent)
{

  TubeGeneratorId=0;
  //std::vector<int> H0CluPadIdsMatrix[24][65];
  //std::vector< std::vector< std::vector<int> > > H0CluPadIdsMatrix;
  
  
  /* 
     for(unsigned int i=0;i<24;i++)
     for(unsigned int j=0;j<65;j++)
     H0CluPadIdsMatrix[i][j].clear();
  */

  //edm::Handle<T2PadClusterCollection> t2padclcoll;
  // iEvent.getByLabel("T2MCl","T2PadClusters",t2padclcoll);
  iEvent.getByLabel(/*"T2MCl"*/CluLabel,"T2PadClusters",t2padclcoll);
  unsigned int planesymb=40; 
  unsigned int ht=5;
  int uniqueclusterid=0;
  double planez=0.;
  occupadcluH0=0.;
  occupadcluH1=0.;
  occupadcluH2=0.;
  occupadcluH3=0.;

  
  Map_ClustuniqueId_CluRoadId.clear();
  Map_ClustuniqueId_CluStatus.clear();
  Map_ClustuniqueId_ISBLOB.clear();
  Map_ClustuniqueId_CluIndexPos.clear();
  Map_PlaneSymb_ClustersuniqueId.clear();
  Map_ClustuniqueId_ZPos.clear();
  Map_ClustuniqueId_PlaneSymb.clear();
  Map_ClustuniqueId_VectorCluRoadId.clear();
  Map_RealClusterID_ToSymbUsed.clear();
  Map_ConcurrentTubeId.clear();  

  Map_cloneRoadID_ToRoadPosIndex.clear();
 
  int rowpad=0; int colpad=0; int convertedcol=0;

  

  for(T2PadClusterCollection::const_iterator itpad= t2padclcoll->begin(); itpad != t2padclcoll->end(); itpad++){
    // vector<T2Cluster> padClv = itpad->second;
    //T2DetId *detID =new T2DetId(itpad->first);
    //std::cout<<"Det loop $$$$$"<<std::endl;
   
    planeinfo=convGeo.GetT2Info(itpad->first);
    planesymb=planeinfo.symb; 
    ht=planesymb/10;
    planez=planeinfo.Zdet;
    bool isblob=false;

    bool planeToExclude=false;
     if(std::find(VplaneToExclude.begin(), VplaneToExclude.end(), planesymb)!=VplaneToExclude.end())
       planeToExclude=true;

     if((t2padclcoll->begin()!=t2padclcoll->end())&&(planeToExclude==false))
    if(std::find(QuartedSelected.begin(), QuartedSelected.end(), ht)!=QuartedSelected.end())
      for(unsigned int k=0; k<itpad->second.size();k++)
	{
	  uniqueclusterid++;
	  Map_ClustuniqueId_CluRoadId[uniqueclusterid]=-1;
	  Map_ClustuniqueId_CluStatus[uniqueclusterid]=-1;//seeds
	  
	  isblob=false;
	  if(itpad->second.at(k).GetNoOfEntries()>=MinCluSize_considered_asBlobs)
	    {
	      isblob=true;
	    }
	  Map_ClustuniqueId_ISBLOB[uniqueclusterid]=isblob;
	  
	  //CluIndexPosition is usefull for save time to get clu information without reiterate. 
	  //since T2PadClusterCollection is <T2DetId, std::vector<T2Cluster> > I need a pair <uint32_t,index in vector>
	  std::pair <uint32_t,unsigned int> apair;
	  apair= std::make_pair(itpad->first,k);
	  
	  Map_ClustuniqueId_CluIndexPos[uniqueclusterid]= apair; 
	  

	  if(useStraightPadTowers){
	
	    if(isblob==false)
	  for(unsigned int hh=0;hh<itpad->second.at(k).GetNoOfEntries();hh++)
	    {
	      //std::cout<<" "<<itpad->second.at(k).GetEntries()[hh].rad_coord<<"-"<<itpad->second.at(k).GetEntries()[hh].ang_coord;
	      rowpad=itpad->second.at(k).GetEntries()[hh].rad_coord;
	      colpad=itpad->second.at(k).GetEntries()[hh].ang_coord;
	      convertedcol=colpad;
	      if(planesymb%2==1)
		convertedcol=64-colpad;
	      
	      if(ht==0){
		H0CluPadIdsMatrix[rowpad][convertedcol].push_back(uniqueclusterid);
		//if((rowpad==1)&&(convertedcol==13))
		//std::cout<<"DebugH0CluPadIdsMatrix. pushback id: "<<uniqueclusterid<<std::endl;
	      }
	     
	      if(ht==1){
		H1CluPadIdsMatrix[rowpad][convertedcol].push_back(uniqueclusterid);
	      }

	      if(ht==2)
		H2CluPadIdsMatrix[rowpad][convertedcol].push_back(uniqueclusterid);

	      if(ht==3)
		H3CluPadIdsMatrix[rowpad][convertedcol].push_back(uniqueclusterid);
	      	      
	    }

	    if(itpad->second.at(k).GetEntries().size()!=itpad->second.at(k).GetNoOfEntries())
	      std::cout<<"CRYTICAL ERROR"<<std::endl;

	    //std::cout<<std::endl;

	  }
	  
	  //if((itpad->second.at(k).GetClusterDR()<0.1)||(itpad->second.at(k).GetClusterDPhi()<1.))
	  //std::cout<<"Error clu dr-dphi:"<<itpad->second.at(k).GetClusterDR()<<" "<<itpad->second.at(k).GetClusterDPhi()<<std::endl;
	
	   if(verbosity>=1)
	    {
	      double rr=itpad->second.at(k).GetClusterR(); double phii=itpad->second.at(k).GetClusterPhi();
	      double xx=rr*cos(phii*3.14159/180.); double yy=rr*sin(phii*3.14159/180.);
	      std::cout<<"Plane symb:"<<planesymb<<" CluId_symb: "<<uniqueclusterid<</*" DetId: "<<apair.first<<*/" VectPos:"<<apair.second<</*". R-Phi:"<<itpad->second.at(k).GetClusterR()<<" "<<itpad->second.at(k).GetClusterPhi()<<*/" X|Y:"<<xx<<"|"<<yy<<" Hrdw_ID: "<<itpad->second.at(k).ClustId_unique<<" CLs:"<<itpad->second.at(k).GetNoOfEntries()<<" Row-Col:";
	      for(unsigned int hh=0;hh<itpad->second.at(k).GetNoOfEntries();hh++)
		{
		  std::cout<<" "<<itpad->second.at(k).GetEntries()[hh].rad_coord<<"-"<<itpad->second.at(k).GetEntries()[hh].ang_coord;

		  

		}
	      std::cout<<std::endl;
	    }


	  Map_PlaneSymb_ClustersuniqueId[planesymb].push_back(uniqueclusterid);

	  Map_RealClusterID_ToSymbUsed[itpad->second.at(k).ClustId_unique]=uniqueclusterid;
	  
	  //if(verbosity>=1)
	  //std::cout<<"Plane symb:"<<planesymb<<" CluId: "<<uniqueclusterid<<" inserted"<<std::endl;

	  Map_ClustuniqueId_ZPos[uniqueclusterid]=planez;
	  Map_ClustuniqueId_PlaneSymb[uniqueclusterid]=planesymb;

	  if(ht==0)
	    occupadcluH0=occupadcluH0+1.0;
	  if(ht==1)
	    occupadcluH1=occupadcluH1+1.0;
	  if(ht==2)
	    occupadcluH2=occupadcluH2+1.0;
	  if(ht==3)
	    occupadcluH3=occupadcluH3+1.0;



	}
  }
  //Angular range should be betweeen 0.05 and 0.3
  //When occupancy per quarter is very high (100 pad cluster per quarter) BiggestTubeAngleConsidered=0.05.
  //A formula could be: BiggestTubeAngleConsidered = 0.05 + (0.3-0.05)(100-numcluinQ)/100.
  
  if((occupadcluH0/10.)>=NumPadCluOccupancyAlert)
    BiggestTubeAngleConsideredH0=0.05;
  else
    BiggestTubeAngleConsideredH0=0.05 + (BiggestTubeAngleConsideredINPUT-0.05)*(100.-(occupadcluH0/10.))/100.;

  if((occupadcluH1/10.)>=NumPadCluOccupancyAlert)
    BiggestTubeAngleConsideredH1=0.05;
  else
    BiggestTubeAngleConsideredH1=0.05 + (BiggestTubeAngleConsideredINPUT-0.05)*(100.-(occupadcluH1/10.))/100.;

  if((occupadcluH2/10.)>=NumPadCluOccupancyAlert)
    BiggestTubeAngleConsideredH2=0.05;
  else
    BiggestTubeAngleConsideredH2=0.05 + (BiggestTubeAngleConsideredINPUT-0.05)*(100.-(occupadcluH2/10.))/100.;

  if((occupadcluH3/10.)>=NumPadCluOccupancyAlert)
    BiggestTubeAngleConsideredH3=0.05;
  else
    BiggestTubeAngleConsideredH3=0.05 + (BiggestTubeAngleConsideredINPUT-0.05)*(100.-(occupadcluH3/10.))/100.;

}

bool T2RoadPadFinder::SmallAngleTube(int candidateCluIdStart, int candidateCluIdStop, const T2PadClusterCollection* T2PadClusterCollectionPtr){

  
  bool toret=false;
  std::pair <uint32_t,unsigned int> T2Det_CluIndex1=Map_ClustuniqueId_CluIndexPos[candidateCluIdStart];
  std::pair <uint32_t,unsigned int> T2Det_CluIndex2=Map_ClustuniqueId_CluIndexPos[candidateCluIdStop];
  T2Cluster cl1;T2Cluster cl2;
  bool proceed1=false; bool proceed2=false;
   //std::cout<<"--B"<<std::endl;
  // std::cout<<"I'm here1"<<std::endl;
  T2DetId T2Det1;
  T2Det1=T2DetId((T2Det_CluIndex1.first));
  T2DetId T2Det2;
  T2Det2=T2DetId((T2Det_CluIndex2.first));
  //std::cout<<"I'm here2 Pad index1-2:"<<(T2Det_CluIndex1.first)<<"   "<<(T2Det_CluIndex2.first)<<std::endl;  
  T2PadClusterCollection::const_iterator AllClusPoint1 = T2PadClusterCollectionPtr->find(T2Det1);
  //T2PadClusterCollection::const_iterator AllClusPoint1 = T2PadClusterCollectionPtr->
  //std::cout<<"I'm here3"<<std::endl;
  T2PadClusterCollection::const_iterator AllClusPoint2 = T2PadClusterCollectionPtr->find(T2Det2);
  if(AllClusPoint1!=T2PadClusterCollectionPtr->end()){
    if(AllClusPoint1->second.size()>(T2Det_CluIndex1.second)){
       cl1=AllClusPoint1->second.at((T2Det_CluIndex1.second));
       proceed1=true;
    }else{std::cout<<"Error in T2RoadPadFinder::New_TwoPoints_Tube Ai"<<std::endl;}
  }else{std::cout<<"Error in T2RoadPadFinder::New_TwoPoints_Tube Bi"<<std::endl;}
  
  if(AllClusPoint2!=T2PadClusterCollectionPtr->end()){
    if(AllClusPoint2->second.size()>(T2Det_CluIndex2.second)){
       cl2=AllClusPoint2->second.at((T2Det_CluIndex2.second));
      proceed2=true;
    }else{std::cout<<"Error in T2RoadPadFinder::New_TwoPoints_Tube Af"<<std::endl;}
  }else{std::cout<<"Error in T2RoadPadFinder::New_TwoPoints_Tube Bf"<<std::endl;}


  if(cl2.GetNoOfEntries()>=MinCluSize_considered_asBlobs)
    proceed2=false;
  
  if(cl1.GetNoOfEntries()>=MinCluSize_considered_asBlobs)
    proceed1=false;
  

  if((proceed2)&&(proceed1))
    {
      double r1=(double)cl1.GetClusterR();    //double Dr1=(double)cl1.GetClusterDR()*(2/sqrt(12.0));
      double r2=(double)cl2.GetClusterR();    //double Dr2=(double)cl2.GetClusterDR()*(2/sqrt(12.0)); 

      double phi1=(double)cl1.GetClusterPhi()*3.14159265/180.;    //double Dphi1=(double)cl1.GetClusterDPhi()*(2/sqrt(12.0))*3.14159265/180.; 
      double phi2=(double)cl2.GetClusterPhi()*3.14159265/180.;    //double Dphi2=(double)cl2.GetClusterDPhi()*(2/sqrt(12.0))*3.14159265/180.; 
	
      double x1=r1*cos(phi1);
      double x2=r2*cos(phi2);

      double y1=r1*sin(phi1);
      double y2=r2*sin(phi2);
      
      double z1=Map_ClustuniqueId_ZPos[candidateCluIdStart];
      double z2=Map_ClustuniqueId_ZPos[candidateCluIdStop];  

    

      if((fabs((x1-x2)/(z1-z2))<BiggestTubeAngleConsidered)&&(fabs((y1-y2)/(z1-z2))<BiggestTubeAngleConsidered))
	toret=true;
      else
	toret=false;
    }
  
  return toret;
}





Tube T2RoadPadFinder::New_TwoPoints_Tube(int candidateCluIdStart, int candidateCluIdStop, bool isbackw,const T2PadClusterCollection* T2PadClusterCollectionPtr){

 //the bigger |Z| Point is stored in the second cluster

   //std::cout<<"--C"<<std::endl;
 
  Tube minitube;
  std::pair <uint32_t,unsigned int> T2Det_CluIndex1=Map_ClustuniqueId_CluIndexPos[candidateCluIdStart];
  std::pair <uint32_t,unsigned int> T2Det_CluIndex2=Map_ClustuniqueId_CluIndexPos[candidateCluIdStop];
  T2Cluster cl1;T2Cluster cl2;
  bool proceed1=false; bool proceed2=false;

  // std::cout<<"I'm here1"<<std::endl;
  T2DetId T2Det1;
  T2Det1=T2DetId((T2Det_CluIndex1.first));
  T2DetId T2Det2;
  T2Det2=T2DetId((T2Det_CluIndex2.first));
  //std::cout<<"I'm here2 Pad index1-2:"<<(T2Det_CluIndex1.first)<<"   "<<(T2Det_CluIndex2.first)<<std::endl;  
  T2PadClusterCollection::const_iterator AllClusPoint1 = T2PadClusterCollectionPtr->find(T2Det1);
  //T2PadClusterCollection::const_iterator AllClusPoint1 = T2PadClusterCollectionPtr->
  //std::cout<<"I'm here3"<<std::endl;
  T2PadClusterCollection::const_iterator AllClusPoint2 = T2PadClusterCollectionPtr->find(T2Det2);
  if(AllClusPoint1!=T2PadClusterCollectionPtr->end()){
    if(AllClusPoint1->second.size()>(T2Det_CluIndex1.second)){
       cl1=AllClusPoint1->second.at((T2Det_CluIndex1.second));
       proceed1=true;
    }else{std::cout<<"Error in T2RoadPadFinder::New_TwoPoints_Tube Ai"<<std::endl;}
  }else{std::cout<<"Error in T2RoadPadFinder::New_TwoPoints_Tube Bi"<<std::endl;}
  
  if(AllClusPoint2!=T2PadClusterCollectionPtr->end()){
    if(AllClusPoint2->second.size()>(T2Det_CluIndex2.second)){
       cl2=AllClusPoint2->second.at((T2Det_CluIndex2.second));
      proceed2=true;
    }else{std::cout<<"Error in T2RoadPadFinder::New_TwoPoints_Tube Af"<<std::endl;}
  }else{std::cout<<"Error in T2RoadPadFinder::New_TwoPoints_Tube Bf"<<std::endl;}

  double x1=0.;
  double x2=0.;

  double y1=0.;
  double y2=0.;

  if((proceed2)&&(proceed1))
    {
      TubeGeneratorId++;
      minitube.uniquetubeId=TubeGeneratorId;
      if(isbackw)//The seed is the bigger |Z| stored in the second cluster.
	{
	  minitube.seedId=candidateCluIdStop;
	  minitube.planef=Map_ClustuniqueId_PlaneSymb[candidateCluIdStop];
	  minitube.planei=Map_ClustuniqueId_PlaneSymb[candidateCluIdStart];
	}
      else
	{
	  minitube.seedId=candidateCluIdStart;
	  minitube.planei=Map_ClustuniqueId_PlaneSymb[candidateCluIdStart];
	  minitube.planef=Map_ClustuniqueId_PlaneSymb[candidateCluIdStop];
	}
      

      
      //I'm not make a special treatment of the blobs in this stage I take the baricentre.

      
	minitube.ClustV.push_back(cl1);   
	minitube.ClustV.push_back(cl2);            
      
	minitube.uniqueCluIdV.push_back(candidateCluIdStart);
	minitube.uniqueCluIdV.push_back(candidateCluIdStop);
      
	double r1=(double)cl1.GetClusterR();    double Dr1=(double)cl1.GetClusterDR()*(2/sqrt(12.0));
	double r2=(double)cl2.GetClusterR();    double Dr2=(double)cl2.GetClusterDR()*(2/sqrt(12.0)); 

	double phi1=(double)cl1.GetClusterPhi()*3.14159265/180.;    double Dphi1=(double)cl1.GetClusterDPhi()*(2/sqrt(12.0))*3.14159265/180.; 
	double phi2=(double)cl2.GetClusterPhi()*3.14159265/180.;    double Dphi2=(double)cl2.GetClusterDPhi()*(2/sqrt(12.0))*3.14159265/180.; 
	
	 x1=r1*cos(phi1);
	 x2=r2*cos(phi2);

	 y1=r1*sin(phi1);
         y2=r2*sin(phi2);
	
	double z1=Map_ClustuniqueId_ZPos[candidateCluIdStart];
	double z2=Map_ClustuniqueId_ZPos[candidateCluIdStop];  

	double Ex1=sqrt(Dr1*Dr1*cos(phi1)*cos(phi1) + Dphi1*Dphi1*r1*r1*sin(phi1)*sin(phi1));
	double Ey1=sqrt(Dr1*Dr1*sin(phi1)*sin(phi1) + Dphi1*Dphi1*r1*r1*cos(phi1)*cos(phi1));
  
	double Ex2=sqrt(Dr2*Dr2*cos(phi2)*cos(phi2) + Dphi2*Dphi2*r2*r2*sin(phi2)*sin(phi2));
	double Ey2=sqrt(Dr2*Dr2*sin(phi2)*sin(phi2) + Dphi2*Dphi2*r2*r2*cos(phi2)*cos(phi2));  

        
	minitube.Zf=z2;minitube.Xf=x2; minitube.Yf=y2; minitube.Rf=r2; minitube.Phif=phi2; 
	minitube.Zi=z1;minitube.Xi=x1; minitube.Yi=y1; minitube.Ri=r1; minitube.Phii=phi1; 

	minitube.VectR.push_back(r1);
	minitube.VectPhi.push_back(phi1);
	minitube.VectX.push_back(x1);
	minitube.VectY.push_back(y1);
	minitube.VectZ.push_back(z1);
	minitube.VectEX.push_back(Ex1); //Error is not yet handled.
	minitube.VectEY.push_back(Ey1);

	minitube.VectR.push_back(r2);
	minitube.VectPhi.push_back(phi2);
	minitube.VectX.push_back(x2);
	minitube.VectY.push_back(y2);
	minitube.VectZ.push_back(z2);
	minitube.VectEX.push_back(Ex2);
	minitube.VectEY.push_back(Ey2);

	//2 Points Tube parameter estimation.
	minitube.Ax=(minitube.Xf-minitube.Xi)/(minitube.Zf-minitube.Zi);
	minitube.Ay=(minitube.Yf-minitube.Yi)/(minitube.Zf-minitube.Zi);
	//Real Error propagation will be made in the MatchConditonAndTubPropagation Stage.

    }
  else
    std::cout<<"Tube not created correctly, FIX THE CODE"<<std::endl;
  
      if(verbosity>=1)
	std::cout<<" begin New_TwoPoints_Tube planes:"<<Map_ClustuniqueId_PlaneSymb[candidateCluIdStart]<<" "<<Map_ClustuniqueId_PlaneSymb[candidateCluIdStop]<<"   x1 y1: "<<x1<<" "<<y1<<"    x2 y2:"<<x2<<" "<<y2<<std::endl;

if(verbosity>=3)
  std::cout<<"Minitube    Zf-Zi:"<<minitube.Zf<<"-"<<minitube.Zi<<" plane i-f:"<<minitube.planei<<"-"<<minitube.planef<<" Slope AX:"<<minitube.Ax<<std::endl;

    
  return minitube;
}


bool T2RoadPadFinder::ZombieForwCompatibility(Tube &TubeForw,Tube &TubeZombie,const T2PadClusterCollection *T2PadClusterCollectionPtr){
  bool closedistance=false;
  
  /*
  std::vector<double> zombieproj=TubePredictedPoint(TubeZombie, TubeForw.Zi); //x,y,ex,ey.
  double xDist=zombieproj.at(0)-TubeForw.Xi;
  double yDist=zombieproj.at(1)-TubeForw.Yi;
  
  if(fabs(xDist)<sqrt(zombieproj.at(2)*zombieproj.at(2)+TolleranceDX*TolleranceDX))
      if(fabs(yDist)<sqrt(zombieproj.at(3)*zombieproj.at(3)+TolleranceDY*TolleranceDY))
	closedistance=true;
  */
  // if(verbosity>=1)
  //std::cout<<"Ckeck ZombieMerge: Coll Ax F-Z:"<<TubeForw.Ax<<" "<<TubeZombie.Ax<<" Coll Ay F-Z"<<TubeForw.Ay<<" "<<TubeZombie.Ay<<std::endl;
  
  if((fabs(TubeForw.Ax-TubeZombie.Ax)<TwoPointsTubesAngularCollinearity)&&(fabs(TubeForw.Ay-TubeZombie.Ay)<TwoPointsTubesAngularCollinearity))//Very row, to be review.
    closedistance=true;


  bool quantumchance=false;
  
  if(closedistance==false){
    //   std::cout<<"Quantum jump will be called from ZombieForwCompatibility"<<TubeForw.ClustV.at(0).GetNoOfEntries()<<std::endl;
    quantumchance= QuantumJumpCheck(TubeForw.ClustV.at(0),TubeForw.uniqueCluIdV.at(0),TubeZombie.uniqueCluIdV,T2PadClusterCollectionPtr);
  }
    
  if(quantumchance)
    closedistance=true;

  return closedistance;
}


//used during recovering of zombie hits
Tube T2RoadPadFinder::MergeSmallTubes(Tube &TubeForw,Tube &TubeZombie){
  
  //I suppose that the vector are ordered in Z
  
  Tube mergedTube;

  mergedTube.seedId=TubeForw.seedId;
  mergedTube.uniquetubeId=TubeForw.uniquetubeId;

  mergedTube.Zi=TubeZombie.Zi;   mergedTube.Yi=TubeZombie.Yi;     mergedTube.Xi=TubeZombie.Xi;
  mergedTube.Zf=TubeForw.Zf;     mergedTube.Yf=TubeForw.Yf;       mergedTube.Xf=TubeForw.Xf;
 
  mergedTube.Ri=TubeZombie.Ri;
  mergedTube.Phii=TubeZombie.Phii;
  mergedTube.Rf=TubeForw.Rf;
  mergedTube.Phif=TubeForw.Phif;

  mergedTube.planei=TubeZombie.planei;
  mergedTube.planef=TubeForw.planef;

  for(unsigned int j=0; j<TubeZombie.ClustV.size();j++)
    {
      mergedTube.ClustV.push_back(TubeZombie.ClustV.at(j));
      mergedTube.uniqueCluIdV.push_back(TubeZombie.uniqueCluIdV.at(j));
 
      mergedTube.VectR.push_back(TubeZombie.VectR.at(j));
      mergedTube.VectPhi.push_back(TubeZombie.VectPhi.at(j));

      mergedTube.VectX.push_back(TubeZombie.VectX.at(j)); 
      mergedTube.VectY.push_back(TubeZombie.VectY.at(j)); 
      mergedTube.VectZ.push_back(TubeZombie.VectZ.at(j)); 
      mergedTube.VectEX.push_back(TubeZombie.VectEX.at(j)); 
      mergedTube.VectEY.push_back(TubeZombie.VectEY.at(j));
    }

  for(unsigned int j=0; j<TubeForw.ClustV.size();j++)
    {

      if(find(mergedTube.uniqueCluIdV.begin(),mergedTube.uniqueCluIdV.end(), (TubeForw.uniqueCluIdV.at(j)) )==mergedTube.uniqueCluIdV.end())
	{
	  mergedTube.ClustV.push_back(TubeForw.ClustV.at(j));
	  mergedTube.uniqueCluIdV.push_back(TubeForw.uniqueCluIdV.at(j));
	  mergedTube.VectR.push_back(TubeForw.VectR.at(j));
	  mergedTube.VectPhi.push_back(TubeForw.VectPhi.at(j));

	  mergedTube.VectX.push_back(TubeForw.VectX.at(j)); 
	  mergedTube.VectY.push_back(TubeForw.VectY.at(j)); 
	  mergedTube.VectZ.push_back(TubeForw.VectZ.at(j)); 
	  mergedTube.VectEX.push_back(TubeForw.VectEX.at(j)); 
	  mergedTube.VectEY.push_back(TubeForw.VectEY.at(j));
	}
    }


return mergedTube;
}


Tube  T2RoadPadFinder::MatchConditonAndTubePropagation(Tube atube,int cluuniqueid,const T2PadClusterCollection* T2PadClusterCollectionPtr,double &distance, bool &ismatching){

  //    Tube atube=atube2;
  ismatching=false;
  //if(verbosity>=1)
  
  double matchdist=0.;
  
  std::vector<double> ActualPointCoordinates;
  double ClusterZ=Map_ClustuniqueId_ZPos[cluuniqueid];
  
  //Then Call the procedure that start from this tube parameters and predict the point  in the clu_uniqueid Plane;
  int numgoodhit=0;
  std::vector<double> propagatedPoint= TubePredictedPoint(atube,ClusterZ,numgoodhit);
  T2Cluster clusterToEventuallyInclude;
  
  ismatching=CheckMatchingDistance(propagatedPoint,cluuniqueid,clusterToEventuallyInclude,atube.uniqueCluIdV,T2PadClusterCollectionPtr,matchdist);
  distance=matchdist;
   double clur=0.;double cluphi=0.; double cluEr=0; double cluEphi=0.; double cluEx=0.; double cluEy=0.;
   double clux=0.;double cluy=0.;  double cluz=Map_ClustuniqueId_ZPos[cluuniqueid];
   if(ismatching)
     {
       if(verbosity>=3)
	 {
	   std::cout<<"Tube Parameters matching at distance "<<distance<<" update tube SIZE "<<atube.ClustV.size();
	   std::cout<<"Cluster "<<clusterToEventuallyInclude.GetClusterR()<<" "<<clusterToEventuallyInclude.GetClusterPhi()<<" is compatible"<<std::endl;
	 }
       atube.ClustV.push_back(clusterToEventuallyInclude);   
       atube.uniqueCluIdV.push_back(cluuniqueid);  
      
      clur=clusterToEventuallyInclude.GetClusterR();
      cluphi=clusterToEventuallyInclude.GetClusterPhi()*3.14159265/180.0;
      cluEr=clusterToEventuallyInclude.GetClusterDR()*(2/sqrt(12.0));
      cluEphi=clusterToEventuallyInclude.GetClusterDPhi()*(2/sqrt(12.0))*3.14159265/180.0;	
      clux=clur*cos(cluphi); cluy=clur*sin(cluphi);     
      //The error on X,Y are from propagation formula	 
      cluEx=sqrt(cluEr*cluEr*cos(cluphi)*cos(cluphi) + cluEphi*cluEphi*clur*clur*sin(cluphi)*sin(cluphi));
      cluEy=sqrt(cluEr*cluEr*sin(cluphi)*sin(cluphi) + cluEphi*cluEphi*clur*clur*cos(cluphi)*cos(cluphi));  
     
      atube.VectR.push_back(clur);
      atube.VectPhi.push_back(cluphi);
      atube.VectX.push_back(clux);
      atube.VectY.push_back(cluy);
      atube.VectZ.push_back(cluz);
      atube.VectEX.push_back(cluEx);
      atube.VectEY.push_back(cluEy);

      //update minimum-maximum plane of the tube
      if(Map_ClustuniqueId_PlaneSymb[cluuniqueid]<atube.planei)
	{
	  atube.planei=Map_ClustuniqueId_PlaneSymb[cluuniqueid];
	  atube.Xi=clux; atube.Yi=cluy; atube.Zi=cluz; atube.Ri=clur; atube.Phii= cluphi;	    
	}
    
      if(Map_ClustuniqueId_PlaneSymb[cluuniqueid]>atube.planef)
	{
	  atube.planef=Map_ClustuniqueId_PlaneSymb[cluuniqueid];
	  atube.Xf=clux; atube.Yf=cluy; atube.Zf=cluz; atube.Rf=clur; atube.Phif=cluphi; 
	}
      if(verbosity>=1)
	std::cout<<" to size "<<atube.ClustV.size()<<" Cluster with x y:"<<clux<<" "<<cluy<<"added"<<std::endl;
    }
   else
     if(verbosity>=3)
       std::cout<<"Cluster "<<clusterToEventuallyInclude.GetClusterR()<<" "<<clusterToEventuallyInclude.GetClusterPhi()<<" does not match the tube"<<std::endl;

  return atube;
}



std::vector<double> T2RoadPadFinder::TubePredictedPoint(Tube &atube, double planeZ,int &countergoodHit)
{
  std::vector<double> PointPrediction;

   unsigned int sizeHitv=atube.VectX.size();

  countergoodHit=0;
   //int hemisphere=0; 
   std::vector<double> vectZGravity;
   vectZGravity.clear();

   double zgrav=0.; double oneoversigmax2=0;
   double lastz=fabs(atube.VectZ[0]);
 
   //blob if size >= MinCluSize_considered_asBlobs
   for(unsigned int i=0;i<sizeHitv;i++)
     {
   
       if(atube.ClustV[i].GetNoOfEntries()<MinCluSize_considered_asBlobs){

	 if(fabs(atube.VectZ[i])>lastz)
	   lastz=fabs(atube.VectZ[i]);
	 
	 countergoodHit++;
	 oneoversigmax2+=1.0/(atube.VectEX[i]*atube.VectEX[i]);
	 zgrav=zgrav+atube.VectZ[i]/atube.VectEX[i]/atube.VectEX[i];
       //  std::cout<<" "<<atube.VectX[i]<<" "<<atube.VectY[i]<<"  Ex:  "<<atube.VectEX[i]<<"   | ";
       }
     }
  
   
   zgrav=zgrav/oneoversigmax2;//((double)sizeHitv);
   //  zgrav=lastz;

   double Sx=0.;
   double Sxz=0.;
   double Szz_x=0.;
   double Sz_x=0.; 
   double S0_x=0.; 

   double Sy=0.;
   double Syz=0.;
   double Szz_y=0.;
   double Sz_y=0.; 
   double S0_y=0.; 
 
   
   double a_xz=0.;
   double b_xz=0.;
   double a_yz=0.;
   double b_yz=0;

  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      if(atube.ClustV[jj].GetNoOfEntries()<MinCluSize_considered_asBlobs){
	//atube.VectEX[jj]=atube.VectEX[jj]*2/3.46;
	Sxz += atube.VectX[jj]*(atube.VectZ[jj]-zgrav)/atube.VectEX[jj]/atube.VectEX[jj];
	Szz_x += (atube.VectZ[jj]-zgrav)*(atube.VectZ[jj]-zgrav)/atube.VectEX[jj]/atube.VectEX[jj];
	Sz_x += (atube.VectZ[jj]-zgrav)/atube.VectEX[jj]/atube.VectEX[jj];
	Sx += atube.VectX[jj]/atube.VectEX[jj]/atube.VectEX[jj];
	S0_x += 1.0/atube.VectEX[jj]/atube.VectEX[jj];
	


	Syz += atube.VectY[jj]*(atube.VectZ[jj]-zgrav)/atube.VectEY[jj]/atube.VectEY[jj];
	Szz_y += (atube.VectZ[jj]-zgrav)*(atube.VectZ[jj]-zgrav)/atube.VectEY[jj]/atube.VectEY[jj];
	Sz_y += (atube.VectZ[jj]-zgrav)/atube.VectEY[jj]/atube.VectEY[jj];
	Sy += atube.VectY[jj]/atube.VectEY[jj]/atube.VectEY[jj];
	S0_y += 1.0/atube.VectEY[jj]/atube.VectEY[jj];
      }
      
    }


  //std::cout<<" Here 2 "<<std::endl;

a_xz = (Sxz*S0_x - Sz_x*Sx) / (Szz_x*S0_x - Sz_x*Sz_x);   // angular coefficient
b_xz = (Sx*Szz_x - Sz_x*Sxz) / (Szz_x*S0_x - Sz_x*Sz_x);  // intercept   X=(a_xz)Z + b_xz

a_yz = (Syz*S0_y - Sz_y*Sy) / (Szz_y*S0_y - Sz_y*Sz_y);   // angular coefficient
b_yz = (Sy*Szz_y - Sz_y*Syz) / (Szz_y*S0_y - Sz_y*Sz_y);  // intercept  Y=(a_yz)Z + b_yz 


/*
double e_a_xz = sqrt( S0_x / (S0_x*Szz_x - Sz_x*Sz_x) );
double e_b_xz = sqrt( Szz_x / (S0_x*Szz_x - Sz_x*Sz_x) );
*/
 double e_a_xz = sqrt( S0_x / ((S0_x*Szz_x - Sz_x*Sz_x)) );
 double e_b_xz = sqrt( Szz_x / ((S0_x*Szz_x - Sz_x*Sz_x)) );

 double e_a_yz = sqrt( S0_y / ((S0_y*Szz_y - Sz_y*Sz_y)) );
 double e_b_yz = sqrt( Szz_y / ((S0_y*Szz_y - Sz_y*Sz_y)) );
/*
//USE T1 Convenction
TVectorD vect(4);
 for(int oo=0; oo<4; oo++)
   vect[oo]=0.;

 vect[0] = b_xz;
 vect[1] = b_yz;
 vect[2] = a_xz;
 vect[3] = a_yz;
*/
double correlx=(-1.0)*Sz_x*(1.0/(S0_x*Szz_x-(Sz_x*Sz_x)));
double correly=(-1.0)*Sz_y*(1.0/(S0_y*Szz_y-(Sz_y*Sz_y)));

 
/*
 TMatrixD mat(4,4);
 for(unsigned int oo=0; oo<4; oo++)
   for(unsigned int ooo=0;ooo<4; ooo++)
     mat[oo][ooo]=0.;


 mat[0][0] = e_b_xz*e_b_xz;
 mat[1][1] = e_b_yz*e_b_yz;
 mat[2][2] = e_a_xz*e_a_xz;
 mat[3][3] = e_a_yz*e_a_yz;
*/

 double xpred=a_xz*(planeZ-zgrav)+b_xz;
 double ypred=a_yz*(planeZ-zgrav)+b_yz;
 double Expred=e_b_xz*e_b_xz + (planeZ-zgrav)*(planeZ-zgrav)*e_a_xz*e_a_xz + 2*(planeZ-zgrav)*correlx;
 double Eypred=e_b_yz*e_b_yz + (planeZ-zgrav)*(planeZ-zgrav)*e_a_yz*e_a_yz + 2*(planeZ-zgrav)*correly;

 

 if((Expred*Eypred)<0)
   std::cout<<"T2RoadPadFinder::TubePredictedPoint error: error on point prediction <0 "<<std::endl;
 
 Expred=sqrt(Expred);
 Eypred=sqrt(Eypred);
 /*
if(((Expred)>40)||((Eypred)>40)){
   std::cout<<"T2RoadPadFinder::TubePredictedPoint Warning: Expred-Eypred:"<<Expred<<"-"<<Eypred<<" NumGoodHit:"<<countergoodHit<<" EAx-EAy: "<<e_a_xz<<" "<<e_a_yz<<" Eb_x Eb_y: "<<e_b_xz<<" "<<e_b_yz<<" Distance:"<<(planeZ-zgrav)<<std::endl;
   for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      if(atube.ClustV[jj].GetNoOfEntries()<MinCluSize_considered_asBlobs){
	std::cout<<"X:"<<atube.VectX[jj]<<" Ex:"<<atube.VectEX[jj]<<" Y:"<<atube.VectY[jj]<<" Ey:"<<atube.VectEY[jj]<<std::endl;
      }
    }
 }
 */


   // x, y, ex, ey
 PointPrediction.push_back(xpred); 
 PointPrediction.push_back(ypred);
 PointPrediction.push_back(Expred); 
 PointPrediction.push_back(Eypred);
    
 //std::cout<<"Prediction X-Y | Ex-Ey "<<xpred<<"-"<<ypred<<" | "<<Expred<<"-"<<Eypred<<"  CorrX:"<<correlx<<"CorrY:"<<correly<<"  E_ax="<<e_a_xz<<"  E_bx="<<e_b_xz<<" Zjump:"<<(planeZ-zgrav)<<std::endl;
 return PointPrediction;
}




double T2RoadPadFinder::GetClusterRadialError(T2Cluster &clu)
{
  double toret=0.;
  
  
  return toret;
}




bool T2RoadPadFinder::CheckMatchingDistance(std::vector<double> &propagatedP,int cluuniqueid, T2Cluster &clusterchecked,std::vector<int> &tubeCluIds,const T2PadClusterCollection *T2PadClusterCollectionPtr, double &distance)
{
  //propagatedP is assumed x-y Ex-Ey.
  
  bool closedistance=false;  
  double cluX; double cluY; double cluR; double cluPhi; double cluEr; double cluEphi;
  double cluEx; double cluEy;

  //unsigned int clusterPlaneSymb=Map_ClustuniqueId_PlaneSymb[cluuniqueid];
  //unsigned int loopSymb=0;
  
  double xDist=0.;
  double yDist=0.;

  //GetCluster.
  std::pair <uint32_t,unsigned int> T2Det_CluIndex1=Map_ClustuniqueId_CluIndexPos[cluuniqueid];
  T2DetId T2Det1; T2Det1=T2DetId((T2Det_CluIndex1.first));
  T2PadClusterCollection::const_iterator AllClusPoint1 = T2PadClusterCollectionPtr->find(T2Det1);
  
  bool proceed1=false;
  
  if(AllClusPoint1!=T2PadClusterCollectionPtr->end()){
    if(AllClusPoint1->second.size()>(T2Det_CluIndex1.second)){       
       proceed1=true;
    }else{std::cout<<"Error in T2RoadPadFinder::CheckMatchingDistance Ai"<<std::endl;}
  }else{std::cout<<"Error in T2RoadPadFinder::CheckMatchingDistance Bi"<<std::endl;}
  


  if(proceed1){

    //clusterchecked=AllClusPoint1->second.at((T2Det_CluIndex1.second));

    if(Map_ClustuniqueId_ISBLOB[cluuniqueid]){
      //Redefine the cluster 
      T2Cluster blobclu=AllClusPoint1->second.at((T2Det_CluIndex1.second));
      
      bool everythingok=false;
      //return false;
      T2Cluster testcluster= PadExtractionFromBlob(blobclu,propagatedP,everythingok);//blobclu passed by address

      
      if(everythingok){
	clusterchecked = testcluster;
	clusterchecked.ClustId_unique=blobclu.ClustId_unique;

	if(verbosity>0)
	  std::cout<<"confirmed: "<<clusterchecked.GetClusterPhi()<<std::endl;
	
      }else{
	//Some problem wiht the blob.
	if(verbosity>0)
	  std::cout<<"WARNING: something wrong in PadExtractionFromBlob "<<std::endl;

	return false;
      }
      

    }else{

      clusterchecked=AllClusPoint1->second.at((T2Det_CluIndex1.second)); //Bug found on 11/1/2011
    }

    cluR=clusterchecked.GetClusterR();
    cluPhi=clusterchecked.GetClusterPhi()*3.14159265/180.0;
    cluEr=clusterchecked.GetClusterDR();
    cluEr=GetClusterRadialError(clusterchecked);
    cluEphi=clusterchecked.GetClusterDPhi()*3.14159265/180.0;	
    cluX=cluR*cos(cluPhi); cluY=cluR*sin(cluPhi);     
    //The error on X,Y are from propagation formula	 
    cluEx=sqrt(cluEr*cluEr*cos(cluPhi)*cos(cluPhi) + cluEphi*cluEphi*cluR*cluR*sin(cluPhi)*sin(cluPhi));
    cluEy=sqrt(cluEr*cluEr*sin(cluPhi)*sin(cluPhi) + cluEphi*cluEphi*cluR*cluR*cos(cluPhi)*cos(cluPhi));      
    cluEx+=TolleranceDX;
    cluEy+= TolleranceDY;
    
    xDist=cluX-propagatedP.at(0);
    yDist=cluY-propagatedP.at(1);
    distance=sqrt(xDist*xDist+yDist*yDist);
    //I suppose that Alignment tollerance are already 
    closedistance=false;
    if(fabs(xDist)<sqrt(NumSigma*NumSigma*propagatedP.at(2)*propagatedP.at(2)+cluEx*cluEx)) //4* is because I'm looking to 2 sigma
      if(fabs(yDist)<sqrt(NumSigma*NumSigma*propagatedP.at(3)*propagatedP.at(3)+cluEy*cluEy))
	closedistance=true;
      

    if((verbosity>=1)/*&&(closedistance)*/)
      {
	std::cout<<"Check Matching with Clu-ID: "<<cluuniqueid<<" .  Point prop X:"<<propagatedP.at(0)<<" pm "<< propagatedP.at(2)<<" Clu X:"<<cluX<<"  pm"<< cluEx<<"   X dist:"<<fabs(xDist)<<" |||  Point prop Y:"<<propagatedP.at(1)<<" pm "<< propagatedP.at(3)<<" Clu Y:"<<cluY<<"    Y dist:"<<fabs(yDist)<<std::endl;
        if( (propagatedP.at(3)>100)||(propagatedP.at(2)>100) )
            {
         	std::cout<<">>>> Too big error in propagation: "<<propagatedP.at(2)<<" "<<propagatedP.at(3)<<std::endl;
             }
    /*
    if(closedistance)
      std::cout<<"Clu Y:"<<cluY<<" CluX:"<<cluX<<" added for propagation. Point prop X:"<<propagatedP.at(0)<<" pm "<< propagatedP.at(2)<<" Clu X:"<<cluX<<"  pm"<< cluEx<<"   X dist:"<<fabs(xDist)<<" |||  Point prop Y:"<<propagatedP.at(1)<<" pm "<< propagatedP.at(3)<<" Clu Y:"<<cluY<<"    Y dist:"<<fabs(yDist)<<std::endl;
    */
      }


    //Always Allows for a quantum Jump from Pad Col-Row in tube vector to new cluster pad-row.
    bool quantumchance=false;
    if((closedistance==false)&&((distance<3*sqrt(cluEx*cluEx+cluEy*cluEy))/*||(tubeCluIds.size()<=3)*/)){
      if(clusterchecked.GetNoOfEntries()<MinCluSize_considered_asBlobs){
	//	std::cout<<"Quantum jump will be called from CheckMatchingDistance"<<std::endl;
	quantumchance= QuantumJumpCheck(clusterchecked,cluuniqueid,tubeCluIds,T2PadClusterCollectionPtr);
      }
    }
    
    if(quantumchance){
      closedistance=true;     
    }

  }

  return closedistance;
}



bool T2RoadPadFinder::QJumpCheckForTubeGen(int candidateCluIdStart, int candidateCluIdStop,const T2PadClusterCollection *T2PadClusterCollectionPtr)
{

std::pair <uint32_t,unsigned int> T2Det_CluIndex1=Map_ClustuniqueId_CluIndexPos[candidateCluIdStart];
  std::pair <uint32_t,unsigned int> T2Det_CluIndex2=Map_ClustuniqueId_CluIndexPos[candidateCluIdStop];
  T2Cluster cl1;T2Cluster cl2;
  bool proceed1=false; bool proceed2=false;

  // std::cout<<"I'm here1"<<std::endl;
  T2DetId T2Det1;
  T2Det1=T2DetId((T2Det_CluIndex1.first));
  T2DetId T2Det2;
  T2Det2=T2DetId((T2Det_CluIndex2.first));
  //std::cout<<"I'm here2 Pad index1-2:"<<(T2Det_CluIndex1.first)<<"   "<<(T2Det_CluIndex2.first)<<std::endl;  
  T2PadClusterCollection::const_iterator AllClusPoint1 = T2PadClusterCollectionPtr->find(T2Det1);
  //T2PadClusterCollection::const_iterator AllClusPoint1 = T2PadClusterCollectionPtr->
  //std::cout<<"I'm here3"<<std::endl;
  T2PadClusterCollection::const_iterator AllClusPoint2 = T2PadClusterCollectionPtr->find(T2Det2);
  if(AllClusPoint1!=T2PadClusterCollectionPtr->end()){
    if(AllClusPoint1->second.size()>(T2Det_CluIndex1.second)){
       cl1=AllClusPoint1->second.at((T2Det_CluIndex1.second));
       proceed1=true;
    }else{std::cout<<"Error in T2RoadPadFinder::New_TwoPoints_Tube Ai"<<std::endl;}
  }else{std::cout<<"Error in T2RoadPadFinder::New_TwoPoints_Tube Bi"<<std::endl;}
  
  if(AllClusPoint2!=T2PadClusterCollectionPtr->end()){
    if(AllClusPoint2->second.size()>(T2Det_CluIndex2.second)){
       cl2=AllClusPoint2->second.at((T2Det_CluIndex2.second));
      proceed2=true;
    }else{std::cout<<"Error in T2RoadPadFinder::New_TwoPoints_Tube Af"<<std::endl;}
  }else{std::cout<<"Error in T2RoadPadFinder::New_TwoPoints_Tube Bf"<<std::endl;}


  bool CloseCluster=false;


  if(cl1.GetNoOfEntries()>=MinCluSize_considered_asBlobs)
    proceed1=false;
  
  if(cl2.GetNoOfEntries()>=MinCluSize_considered_asBlobs)
    proceed2=false;

  if((proceed2)&&(proceed1))
    {
      int row1=0;int col1=0; 
      int row2=0;int col2=0;
      
      int cluplane1=Map_ClustuniqueId_PlaneSymb[candidateCluIdStart];
      int cluplane2=Map_ClustuniqueId_PlaneSymb[candidateCluIdStop];
  
  for(unsigned int i=0;i<cl1.GetNoOfEntries();i++)
    {
      row1=cl1.GetEntries()[i].rad_coord;
      col1=cl1.GetEntries()[i].ang_coord;

     
      for(unsigned int j=0;j<cl2.GetNoOfEntries();j++)
	{
	  
	  row2=cl2.GetEntries()[j].rad_coord;
	  col2=cl2.GetEntries()[j].ang_coord;

	  
	  int odd_evencheck=cluplane1-cluplane2;//cluplane
	  if(abs(odd_evencheck)%2==0){
	    //It means that planes for tubeCluIds.at(j) and clusterchecked are both Even or Odd so you don't need convertion
	    if((((abs(row2-row1))<=1)&&((abs(col2-col1))<=1))/*||(((abs(col2-col1))==0)&&((abs(row2-row1))==0))*/)
	      CloseCluster=true;
	  }
	  else
	    {
	      col2=64-col2;
	      if((((abs(row2-row1))<=1)&&((abs(col2-col1))<=1))/*||(((abs(col2-col1))==0)&&((abs(row2-row1))==0))*/)
		CloseCluster=true;
	    }
	}

    }
    }

  return CloseCluster;

}

//Review this: we want to allow quantum jump !only! if the prediction is parallel and you fall in digi domains.
bool T2RoadPadFinder::QuantumJumpCheck(T2Cluster &clusterchecked, int cluuniqueid, std::vector<int> &tubeCluIds,const T2PadClusterCollection *T2PadClusterCollectionPtr)
{
  
  bool CloseCluster=false;

  unsigned int cluplane=Map_ClustuniqueId_PlaneSymb[cluuniqueid];
  unsigned int cluintubeid=0;

  //find the closest cluster in tubeCluIds to cluplane
  int distance=50;// int tubecluIdtocompare=-10;
  int Actdistance=50;
  int selcluId_tube=tubeCluIds.at(0);

  
  int totalcluintube=tubeCluIds.size();
  int counterclu=0;

  bool onegoodclufound=false;
  T2Cluster closestClu;
  // I don't want to use blobs as the closerplane. So I'm doing a while loop until I get the Z-closer and with a small clusize
  while((onegoodclufound==false)&&(counterclu<totalcluintube)){
    
  
    for(unsigned int j=counterclu;j<tubeCluIds.size();j++)
      {
	cluintubeid=Map_ClustuniqueId_PlaneSymb[tubeCluIds.at(j)];
	Actdistance=(cluintubeid-cluplane);
	
	if(abs(Actdistance)<distance)
	  {
	    distance=Actdistance; 
	    selcluId_tube=tubeCluIds.at(j);
	  }
      }
    
    // unsigned int closeClutubePlane=Map_ClustuniqueId_PlaneSymb[selcluId_tube];

    //GetCluster.
    T2Cluster closestCluTest;
    std::pair <uint32_t,unsigned int> T2Det_CluIndex1=Map_ClustuniqueId_CluIndexPos[selcluId_tube];
    T2DetId T2Det1; T2Det1=T2DetId((T2Det_CluIndex1.first));
    T2PadClusterCollection::const_iterator AllClusPoint1 = T2PadClusterCollectionPtr->find(T2Det1);
    bool proceed1=false;
    
    if(AllClusPoint1!=T2PadClusterCollectionPtr->end()){
      if(AllClusPoint1->second.size()>(T2Det_CluIndex1.second)){       
	proceed1=true;
	closestCluTest= AllClusPoint1->second.at((T2Det_CluIndex1.second));
      }else{std::cout<<"Error in QuantumJumpCheck::CheckMatchingDistance Ai"<<std::endl;}
    }else{std::cout<<"Error in QuantumJumpCheck::CheckMatchingDistance Bi"<<std::endl;}
    
    counterclu++;

    if(proceed1 && (closestClu.GetNoOfEntries()<MinCluSize_considered_asBlobs)){
      onegoodclufound=true;
      closestClu=closestCluTest;
    }
  
  }//end while
  
  int row1=0;int col1=0; 
  int row2=0;int col2=0;

  if(onegoodclufound){
    for(unsigned int i=0;i<clusterchecked.GetNoOfEntries();i++)
      {
	row1=clusterchecked.GetEntries()[i].rad_coord;
	col1=clusterchecked.GetEntries()[i].ang_coord;
	
     
	for(unsigned int j=0;j<closestClu.GetNoOfEntries();j++)
	  {
	    
	    row2=closestClu.GetEntries()[j].rad_coord;//bug found on 3 jan 010. i instead of j
	    col2=closestClu.GetEntries()[j].ang_coord;

	  
	    int odd_evencheck=cluintubeid-cluplane;//cluplane
	    if(abs(odd_evencheck)%2==0){
	      //It means that planes for tubeCluIds.at(j) and clusterchecked are both Even or Odd so you don't need convertion
	      if((((abs(row2-row1))<=1)&&((abs(col2-col1))<=1))/*||(((abs(col2-col1))==0)&&((abs(row2-row1))==0))*/)
		CloseCluster=true;
	    }
	    else
	      {
		col2=64-col2;
		if((((abs(row2-row1))<=1)&&((abs(col2-col1))<=1))/*||(((abs(col2-col1))==0)&&((abs(row2-row1))==0))*/)
		  CloseCluster=true;
	      }
	  }
	
      }
    //@@@@ Phi:62.0368@@@@ Phi:62.0334@@@@ Phi:280.189
    
    if(CloseCluster){
      // std::cout<<"Cluster PhiClose-PhiBeg:"<<closestClu.GetClusterPhi()<<" "<<clusterchecked.GetClusterPhi()<<" merged for quantum jump"<<std::endl;
      if(fabs(closestClu.GetClusterPhi()-clusterchecked.GetClusterPhi())>25.)
	{
	  std::cout<<"WARNING!!: Qjump with too big dphi. Planes: "<<cluintubeid<<" "<<cluplane<<"  rows:"<<row1<<" "<<row2<<" cols:"<<col1<<" "<<col2<<std::endl;
	  if(fabs(closestClu.GetClusterPhi())<0.01)
	    std::cout<<"||| Size ClosestCluFound:"<<closestClu.GetNoOfEntries()<<std::endl;
	  if(fabs(clusterchecked.GetClusterPhi())<0.01)
	    std::cout<<"||| Size BeginningClu:"<<clusterchecked.GetNoOfEntries()<<std::endl;
	}
    }
    
    if(verbosity>=3){
      if(CloseCluster)  
	std::cout<<"Quantum-jump satisfied"<<std::endl;
      else
	std::cout<<"Quantum-jump fail since:r1-r2 | c1-c2 are: "<<row1<<" "<<row2<<" | "<<col1<<" "<<col2<<std::endl;
    }
  }
  else
    CloseCluster=false;
  
  return CloseCluster;
}




T2Cluster T2RoadPadFinder::PadExtractionFromBlob(T2Cluster &blobclu, std::vector<double> &propagatedP, bool &allOK)
{

  //Set Geometry of the proper plane and found the best pad.
  int row=0; int col=0; double myphimin=0.;double myRmin=0.;  double myRmax=0.; 
  double x=0.; double y=0.; double myphi=0.;double myR=0.;
  double mind=300.; double d=0.;
  int minrow=0;int mincol=0;
  T2ROGeometry t2rogeo(blobclu.GetDetID());
  //std::cout<<"blob size:"<<blobclu.GetNoOfEntries()<<std::endl;
  for(unsigned int i=0;i<blobclu.GetNoOfEntries();i++)
    {
      row=blobclu.GetEntries()[i].rad_coord;
      col=blobclu.GetEntries()[i].ang_coord;
      
	 
      myphi=t2rogeo.GetPadPhiMin(row,col);
      // double myphimax=t2rogeo.GetPadPhiMin(row,col);
      myR=t2rogeo.GetPadRMin(row,col);   //bug found 11/1/2011
      x=myR*cos(myphi*3.14159265/180.0);
      y=myR*sin(myphi*3.14159265/180.0);
      //std::cout<<"blob extraction x-y:"<<x<<" "<<y<<" exp xy:"<<propagatedP.at(0)<<propagatedP.at(1)<<std::endl;
      d=sqrt((x-propagatedP.at(0))*(x-propagatedP.at(0))+(y-propagatedP.at(1))*(y-propagatedP.at(1)));
      if(d<mind)
	{
	  
	  myphimin=myphi; 
	  myRmin=myR;  myRmax=t2rogeo.GetPadRMax(row,col);
	  mind=d;
	  minrow=row; mincol=col;
	}
    }
  
  //Build a cluster with minrow, mincol.

  if(verbosity>0)
    std::cout<<"PadExtractionFromBlob: minD:"<<mind<<" Phi:"<<myphimin<<" size:"<<blobclu.GetNoOfEntries()<<std::endl;

  if(mind>=300.)//300 is 150*2
    allOK=false;
  else
    allOK=true;

  T2Cluster SinglePadCluster;
  cluster_entry onePad; 
  onePad= cluster_entry(minrow,mincol);
  SinglePadCluster.AddEntry(onePad); SinglePadCluster.SetNoOfEntries(1);
  SinglePadCluster.SetRadialSpread(1); 
  SinglePadCluster.SetRadialCentrePos(minrow); 
  SinglePadCluster.SetAngularSpread(1); 
  SinglePadCluster.SetAngularCentrePos(mincol); 
  SinglePadCluster.SetClusterR(myRmin); 
  SinglePadCluster.SetClusterDR(fabs(myRmin-myRmax)/2.0); 
  SinglePadCluster.SetClusterPhi(myphimin);
  SinglePadCluster.SetClusterDPhi(fabs(3.5)); 
  SinglePadCluster.SetDetID(blobclu.GetDetID()); 


  return SinglePadCluster;
}



void T2RoadPadFinder::MakethreeWithThisPlaneSeed(unsigned int seedplane,std::vector<Tube> &AllInitialTubes, std::vector<Tube> &AllFinalTubes,const T2PadClusterCollection*  T2PadClusterCollectionPtr)
{
  
  //Basic steps for the "seedplane" n: 
  //a)propagate the tubes from this plane to the forward direction (if plane is 1,2...).
  //b)check criteria on tube: if the tubes has the closer pad with a distance bigger than 2 planes from this plane and does not
  //  match in the current plane, and if its size has smaller than the minimum
  //  it has to be killed and removed frome the tubes. The status of the composing pads will be updated as 3 (killed)
  //c)New Tubes will be created from this plane to the next one, using as seeds only pads not already utilized in some 
  //  alive tube. There will be also a backward propagation to include also dead or unassociated seeds.
  //d) Make an unique tube collection for tubes performed propagation with the newtubes created from the plane seed  



  int cluuniqueid=0;
  unsigned int seedplane_localnumber=seedplane%10;
  unsigned int firstplane_number=(seedplane/10)*10;
  
  //Warning propagation step or matching backward can happens only from plane 2 
if(verbosity>=1)
    std::cout<<"-------------------------------------------------------------------------PLANE LEVEL"<<seedplane<<std::endl;

  
  if((seedplane_localnumber>=1)&&(seedplane_localnumber<9)){
    
if(verbosity>=1)
    std::cout<<"Tubes propagation from plane"<<seedplane<<std::endl;


  // a)propagate the tubes from this plane to the forward(seedplane is called from 0 to 10)
     std::map<int,std::vector<unsigned int> > tubeindex_padidV;
      std::map<int,unsigned int> padid_tubeindex; //map to add after the matching conditions and the loop, the cluster to the Tube
      std::map<unsigned int, vector<int> >::iterator itm;
      itm=Map_PlaneSymb_ClustersuniqueId.find(seedplane+1);
      std::vector<unsigned int> tubeIndexToEventuallykill;
      unsigned int ConcurrentTubeCounter=0;

      std::vector<Tube> addedConcurrentTubes;addedConcurrentTubes.clear();
      
      for(unsigned int tubind=0;tubind<AllInitialTubes.size();tubind++)
	//if(tubind==13)
	{
	  if(verbosity>=1)
	    std::cout<<"Propagation of tube"<<tubind<<" with size "<<AllInitialTubes.at(tubind).uniqueCluIdV.size()<<" whit final plane "<<AllInitialTubes.at(tubind).planef<<" (Initial R-Phi | X Y= "<<AllInitialTubes.at(tubind).Ri<<"-"<<AllInitialTubes.at(tubind).Phii*180.0/3.14159<<"  | Xi:"<<AllInitialTubes.at(tubind).Xi<<" Yi:"<<AllInitialTubes.at(tubind).Yi<<")"<<" (final R-Phi | X Y= "<<AllInitialTubes.at(tubind).Rf<<"-"<<AllInitialTubes.at(tubind).Phif*180.0/3.14159<<"  | X:"<<AllInitialTubes.at(tubind).Xf<<" Y:"<<AllInitialTubes.at(tubind).Yf<<")"<<"       to plane "<<(seedplane+1)<<std::endl;
	  bool tubematchcondition=false;//bool thismatchingcondition=false;
	  
	  if(AllInitialTubes.at(tubind).planef>seedplane)
	    std::cout<<"CODE ERROR IN PROPAGATON"<<std::endl;
	  
	  ConcurrentTubeCounter=0;
	  if(seedplane-AllInitialTubes.at(tubind).planef<=InefficiencyMaxJump){

	    // unsigned int Numcluinplane=itm->second.size();
	    //	    std::cout<<(seedplane+1)<<" has "<<Numcluinplane<<" cluster "<<std::endl;
	    double alldist=0.;/*int allcuunique=0;*/ 
	    int bestcluunique=0; double bestdist=1000.;
	    
	    if(itm!=Map_PlaneSymb_ClustersuniqueId.end()){

	      for(unsigned int cl=0;cl<itm->second.size();cl++)
		if(((Map_ClustuniqueId_CluStatus[itm->second.at(cl)]==-1)&&(AllowsPadReAssociation==false))||(AllowsPadReAssociation==true))
		  {
		    bool thismatchingcondition=false;
		    cluuniqueid=itm->second.at(cl);	
		    //cluplaneid=Map_ClustuniqueId_PlaneSymb[cluuniqueid];

		    Tube atube=AllInitialTubes.at(tubind);
		    Tube newTube=MatchConditonAndTubePropagation(atube,cluuniqueid, T2PadClusterCollectionPtr,alldist,thismatchingcondition); //If It match tube will be automatically updated.		
		  
		    if(thismatchingcondition)
		      {
			
			if(alldist<bestdist)
			  {
			    tubematchcondition=true;
			    bestdist=alldist;
			    bestcluunique=cluuniqueid;
			    //std::cout<<" "<<bestdist<<" "<<bestcluunique<<std::endl;
			  }
		      
		      //  if(verbosity>=1)
		      //std::cout<<"Tube match. Tube size:"<<atube.uniqueCluIdV.size()<<std::endl;	    	
		      
		      //padid_tubeindex[cluuniqueid]=tubind;
		      //tubeindex_padidV[tubind].push_back(cluuniqueid);
		      
		      //changestatus of the clust
		      //-- Map_ClustuniqueId_CluStatus[cluuniqueid]=1;//Associated to some road
		      //AllInitialTubes.at(tubind)=atube;
		      //continue;//could not be 2 match per one tube
		      }

		    if((AllowsConcurrentBranches==true)&&(thismatchingcondition==true))
		      {//Allow concurrent tubes
			if(verbosity>=1)
			  std::cout<<"@@@Tube match. Tube size:"<<newTube.uniqueCluIdV.size()<<"->";
			if(ConcurrentTubeCounter==0){
			  addedConcurrentTubes.push_back(newTube);
			  //AllInitialTubes.at(tubind)=newTube;
			  ConcurrentTubeCounter++; //bug found on 3/1/011
			}
			else
			  {
			    //AllInitialTubes.push_back(newTube); //same Id of the tube.
			    addedConcurrentTubes.push_back(newTube);
			    ConcurrentTubeCounter++;
			    Map_ConcurrentTubeId[newTube.uniquetubeId]=ConcurrentTubeCounter;			    
			  }
			Map_ClustuniqueId_CluStatus[cluuniqueid]=1;//Associated to some road
		      }
		  
		  }//end if-for on target plane clusters.
	    
	      // std::cout<<"B"<<std::endl;
	      //WARNING: save somehow your tube or it will be lost!
	      if((tubematchcondition)&&(AllowsConcurrentBranches==false))
		{
		  //std::cout<<"C"<<std::endl;
		  padid_tubeindex[bestcluunique]=tubind;
		  tubeindex_padidV[tubind].push_back(bestcluunique);
		  //std::cout<<"D"<<std::endl;
		  double gdist=0.;
		  Tube ChosedPropagatedTube=AllInitialTubes.at(tubind);
		  //std::cout<<"E"<<std::endl;
		  if(verbosity>=1){
		    std::cout<<"@@@Tube match. Tube size:"<<ChosedPropagatedTube.uniqueCluIdV.size()<<"->";
		  }
		  //std::cout<<"F -> id clu to add:"<<bestcluunique<<std::endl;
		  Tube newTube=MatchConditonAndTubePropagation(ChosedPropagatedTube,bestcluunique, T2PadClusterCollectionPtr,gdist,tubematchcondition);
		  /*
		  for(unsigned int ig=0;ig<newTube.uniqueCluIdV.size();ig++){
		    std::cout<<"@@@@ Phi:"<<newTube.ClustV.at(ig).GetClusterPhi();
		  }
		  std::cout<<std::endl;
		  */

		  AllInitialTubes.at(tubind)=newTube; Map_ClustuniqueId_CluStatus[bestcluunique]=1;//Associated to some road
		  //  std::cout<<ChosedPropagatedTube.uniqueCluIdV.size()<<std::endl;
		}
	      // else
		//if(AllowsPadReAssociation==false)
		//if(verbosity>=1)
		//  std::cout<<"@@@Tube Not match. Tube size:"<<ChosedPropagatedTube.uniqueCluIdV.size()<<"->";
	      
	    }//if there is something in the target plane

	  }else{ 

	    if(verbosity>=1)
	      std::cout<<"Tube "<<tubind<<" can be killed or win in plane"<<seedplane<<std::endl;

	    // This tube should be killed If it hasn't the final requirement 
	    //since does not match with any plane pads of the last 2 planes.
	    tubeIndexToEventuallykill.push_back(tubind);
	  }
if(verbosity>=1)
  std::cout<<"At the end of the propagation the size of the tube "<<tubind<<" is: "<<AllInitialTubes.at(tubind).uniqueCluIdV.size()<<std::endl;

	} //End of for AllInitialTubes
      
      for(unsigned int jj=0;jj<addedConcurrentTubes.size();jj++)
	AllInitialTubes.push_back(addedConcurrentTubes.at(jj));

      //  b)check criteria on tube: if the tubes has the closer pad with a distance bigger than 2 planes and does not
      //  match in the current plane, and if its size has smaller than the minimum
      //  it has to be killed and removed frome the tubes. The status of the composing pads will be updated as 3 (killed).
      //  if the tube as the quality criteria it can be stored as a FinalTube since the algorithm criteria don't 
      //  gives further possibility to change its status.
      unsigned int tubcluid=0;  

    for(unsigned int u=0;u<tubeIndexToEventuallykill.size();u++)
	{
	  unsigned int ti=tubeIndexToEventuallykill.at(u);
	  //  std::cout<<"tube "<<ti<<"kill or win"<<std::endl;
	  if((int)AllInitialTubes.at(ti).uniqueCluIdV.size()<Nmin_padsFinal){ // Tube killed
	    // Put its cluster as dead
	    for(unsigned int clu=0;clu<AllInitialTubes.at(u).uniqueCluIdV.size();clu++)
	      {
		tubcluid=AllInitialTubes.at(u).uniqueCluIdV.at(clu);
		Map_ClustuniqueId_CluStatus[tubcluid] =3;
		//	if(verbosity>=1)std::cout<<"...tube "<<ti<<" dead"<<std::endl;
	      }
	    //Note: It is not important if it is already shared with another road. In future eventually one can associate to each cluster
	    // a vector of status.
	    //AllInitialTubes.erase(AllInitialTubes.begin()+ti);
	  }
	  else{// Tube stored as final and removed from AllInitialTubes

	    WinningTubes.push_back(AllInitialTubes.at(ti));
	    if(verbosity>=1)std::cout<<"...tube "<<ti<<" STORED AS GOOD, with size"<<AllInitialTubes.at(ti).uniqueCluIdV.size()<<std::endl;
	    //AllInitialTubes.erase(ti);
	  }
	}

if(verbosity>=1)
      std::cout<<"Removing dead tubes...: "<<std::endl;
    for(unsigned int u=0;u<AllInitialTubes.size();u++)
	{
	  if(std::find(tubeIndexToEventuallykill.begin(),tubeIndexToEventuallykill.end(),u)==tubeIndexToEventuallykill.end())
	    {
	      AllFinalTubes.push_back(AllInitialTubes.at(u));
	    }
	}
if(verbosity>=1)
  {
      std::cout<<"Tube collection size at the end of the propagation: "<<AllFinalTubes.size()<<std::endl;
     
  }

  }//if(seedplane_localnumber>=1 and <9 )


  //... but if the plane==9 copy AllinitialTubes in AllFinalTubes
  if(seedplane_localnumber==9)
   for(unsigned int u=0;u<AllInitialTubes.size();u++)
     {
       AllFinalTubes.push_back(AllInitialTubes.at(u));
     }






  //STEP c)New Tubes will be created from this plane to the next one, using as seeds only pads not already utilized in some 

  //Match ALL Seeds With Next Plane (or with backword zombie.)
  std::vector<Tube> NewForwTubes;    
  if(seedplane_localnumber<9)    //Bug found: was (seedplane_localnumber<=9)
    {

if(verbosity>=1)
      std::cout<<"Tubes creation from plane"<<seedplane<<"..."<<std::endl;


      int candidateCluIdStart=0; int candidateCluIdStop=0;
      unsigned int targetplane=0;
      std::vector<Tube> ZombieTubes; 
      
      std::map<unsigned int, std::vector<int> >::iterator Map_PlaneSymb_ClustersuniqueId_IT;
      Map_PlaneSymb_ClustersuniqueId_IT=Map_PlaneSymb_ClustersuniqueId.find(seedplane);
      
      if(Map_PlaneSymb_ClustersuniqueId_IT!=Map_PlaneSymb_ClustersuniqueId.end())
	{
	  //loop over the seedplane plane
	  for(unsigned int i=0;i<Map_PlaneSymb_ClustersuniqueId_IT->second.size();i++)
	    {
	       
	      candidateCluIdStart=Map_PlaneSymb_ClustersuniqueId_IT->second.at(i);
	      //Here You can have some efficiency dependence from the angular tube
	      if((Map_ClustuniqueId_CluStatus[candidateCluIdStart]==-1)||(Map_ClustuniqueId_ISBLOB[candidateCluIdStart]))
		{
		  //Forward generation
		  //Look Forward (jump Only one) and collegate everything.
		  targetplane=seedplane+1;
		  if(verbosity>=3)
		    std::cout<<"starting Cluster id: "<<candidateCluIdStart<<std::endl;


		  // Map_ClustuniqueId_CluStatus[candidateCluIdStart]=1;
		  std::map<unsigned int, std::vector<int> >::iterator Map_PlaneSymb_ClustersuniqueId_IT_targ;
		  Map_PlaneSymb_ClustersuniqueId_IT_targ=Map_PlaneSymb_ClustersuniqueId.find(targetplane);
		  if(Map_PlaneSymb_ClustersuniqueId_IT_targ!=Map_PlaneSymb_ClustersuniqueId.end()){	  
		    for(unsigned int i=0;i<Map_PlaneSymb_ClustersuniqueId_IT_targ->second.size();i++){

		     
		      candidateCluIdStop=Map_PlaneSymb_ClustersuniqueId_IT_targ->second.at(i);
		       

		       if(SmallAngleTube(candidateCluIdStart,candidateCluIdStop,T2PadClusterCollectionPtr)){

			 if(verbosity>=3)
			   std::cout<<" ... to target plane "<<targetplane<<" with clu-id "<<candidateCluIdStop<<"    ";
			 
			 Tube atube; //std::cout<<candidateCluIdStop<<"||||"<<std::endl;
			 atube.seedId=candidateCluIdStart;
			 atube=New_TwoPoints_Tube(candidateCluIdStart,candidateCluIdStop,false,T2PadClusterCollectionPtr);
			 NewForwTubes.push_back(atube);
			 //->update seed status 
			 Map_ClustuniqueId_CluStatus[candidateCluIdStop]=1;//Associated to some road
			 Map_ClustuniqueId_CluStatus[candidateCluIdStart]=1; //Bug found in 21/12/2010
		       }else{
			 
			 bool testqjump=QJumpCheckForTubeGen(candidateCluIdStart,candidateCluIdStop,T2PadClusterCollectionPtr);
			 if(testqjump){
			   Tube atube; //std::cout<<candidateCluIdStop<<"||||"<<std::endl;
			   atube.seedId=candidateCluIdStart;
			   atube=New_TwoPoints_Tube(candidateCluIdStart,candidateCluIdStop,false,T2PadClusterCollectionPtr);
			   NewForwTubes.push_back(atube);
			   //->update seed status 
			   Map_ClustuniqueId_CluStatus[candidateCluIdStop]=1;//Associated to some road
			   Map_ClustuniqueId_CluStatus[candidateCluIdStart]=1; //Bug found in 21/12/2010
			
			 }
			    
		       }


		       //else
		       //if(seedplane==10)
		       //  std::cout<<"Too large angle for clu-id: "<<candidateCluIdStop<<" "<<candidateCluIdStart<<std::endl;
		    }
		  }
		  
		  //Look Backword (steps MaxJump.) to check possible Zombie seeds (id=-1)to re-associate.
		  //Backword generation: in order to be fully efficient and without making double-counts tube one Should:
		  //a) make a new tube containig the zombie.
		  //b) see if the zombie-tube are compatible with one the new tubes. If so the 3 points tube should be done.
		  targetplane=seedplane-2; 
		  unsigned int NumPlanesdeepBackw=1;//How many backw plane to recheck.
		  unsigned int backwcount=0;
		 
		  while((targetplane>=firstplane_number)&&(backwcount<NumPlanesdeepBackw))//
		   {
		     backwcount++; //std::cout<<" Zombie Recover test0"<<std::endl;
		     std::map<unsigned int, std::vector<int> >::iterator Map_PlaneSymb_ClustersuniqueId_IT_targ;
		     Map_PlaneSymb_ClustersuniqueId_IT_targ=Map_PlaneSymb_ClustersuniqueId.find(targetplane);
		      if(Map_PlaneSymb_ClustersuniqueId_IT_targ!=Map_PlaneSymb_ClustersuniqueId.end()){	  
			for(unsigned int i=0;i<Map_PlaneSymb_ClustersuniqueId_IT_targ->second.size();i++){
			  candidateCluIdStop=Map_PlaneSymb_ClustersuniqueId_IT_targ->second.at(i);
			   //Backword seed should be a death (killed=3) or never used (-1) 
			  //std::cout<<" Zombie Recover test1 status zomb-hit:"<<Map_ClustuniqueId_CluStatus[candidateCluIdStop]<<std::endl;
			  if((Map_ClustuniqueId_CluStatus[candidateCluIdStop]==3)||(Map_ClustuniqueId_CluStatus[candidateCluIdStop]==-1))
			    {


			       if(SmallAngleTube(candidateCluIdStart,candidateCluIdStop,T2PadClusterCollectionPtr)){

				 Tube atube;
				 if(verbosity>=1)
				   {
				     std::cout<<" Looking if I can Recover a Zombie hit from plane"<<targetplane<<std::endl;
				     std::cout<<"Target Clu id: "<<candidateCluIdStop<<" has status "<<Map_ClustuniqueId_CluStatus[candidateCluIdStop]<<std::endl;
				   }
				 atube=New_TwoPoints_Tube(candidateCluIdStop,candidateCluIdStart,true,T2PadClusterCollectionPtr);//First argument always the small |Z| clust
				 ZombieTubes.push_back(atube);
				 //->update seed status
				 Map_ClustuniqueId_CluStatus[candidateCluIdStop]=1;
			       }
			    }
			}
		      }
		      targetplane--;
		    }
		}
	      //else
	      //if(seedplane==10)
	      //  std::cout<<"This cluster is not free"<<std::endl;
	    }
	  //end loop on the seeds in the seed plane
	  if(verbosity>=1)
	    std::cout<<" all seed review: #New tubes="<<NewForwTubes.size()<<" Zombie:"<<ZombieTubes.size()<<std::endl;
	  //Eventually merge Zombie inside Forward; First of all find the tubes with the same seed
	  std::map<int,int> tubeidFw_ZombtoMerge; std::map<int,int>::iterator mapmergeIT ;
	  for(unsigned int z=0;z<ZombieTubes.size();z++)
	  for(unsigned int f=0;f<NewForwTubes.size();f++)	    
	      {
		
		if(ZombieTubes.at(z).seedId==NewForwTubes.at(f).seedId)		  
		  {
		    //  if(verbosity>=1)
		    //std::cout<<"Seed Forw-Zombie Ok... check Collinearity:"<<std::endl;
		    //but you have to check if the zombie is collinear!
		    
		    if(ZombieForwCompatibility(NewForwTubes.at(f),ZombieTubes.at(z),T2PadClusterCollectionPtr))
		      {
			//if(verbosity>=1)
			//std::cout<<" Found a Zombie and a Forward to merge"<<std::endl;
			tubeidFw_ZombtoMerge[z]=f;		      
		      }
		  }
	      }
	  
	  for(unsigned int z=0;z<ZombieTubes.size();z++)
	    {
	      mapmergeIT=tubeidFw_ZombtoMerge.find(z);
	      
	      if(mapmergeIT!=tubeidFw_ZombtoMerge.end())//Merge Zombie and Forward
		{
		  Tube MergedTube=MergeSmallTubes(NewForwTubes.at((*mapmergeIT).second),ZombieTubes.at(z));
		  //std::cout<<"merged tube has size:"<<MergedTube.uniqueCluIdV.size()<<std::endl;
		  NewForwTubes.at((*mapmergeIT).second)=MergedTube;
		}
	      //Erase now The zombie but does not matter
	      //....
	    }

	  for(unsigned int z=0;z<ZombieTubes.size();z++)
	    {
	      mapmergeIT=tubeidFw_ZombtoMerge.find(z);
	      
	      if(mapmergeIT==tubeidFw_ZombtoMerge.end())//Classify Zombie as new forward
		{
		  NewForwTubes.push_back(ZombieTubes.at(z));
		}
	      
	    }

	 
	}
      else
if(verbosity>=1)
	std::cout<<"... This plane does not have cluster..."<<std::endl;
      
    }

  //------------------------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------------------------

  //d) Make an unique tube collection for tubes performed propagation with the newtubes created from the plane seed  
  for(unsigned int NewForwTubesInd=0;NewForwTubesInd<NewForwTubes.size();NewForwTubesInd++) 
    {
      AllFinalTubes.push_back(NewForwTubes.at(NewForwTubesInd));
    }
  if(verbosity>=1)
    {
      std::cout<<"At the end of the plane level there are "<<AllFinalTubes.size()<<" candidate tubes, with sizes:"<<std::endl;
      for(unsigned int NewForwTubesInd=0;NewForwTubesInd<AllFinalTubes.size();NewForwTubesInd++) 
	std::cout<<" "<<AllFinalTubes.at(NewForwTubesInd).uniqueCluIdV.size()<<"  ";
      std::cout<<std::endl;
    }
  //Next call will work with AllFinalTubes
  
  
}



void T2RoadPadFinder::StraightTubeFinder(int quarter, std::vector<Tube> &WinningTubes,const T2PadClusterCollection*  T2PadClusterCollectionPtr){
    

    unsigned int StraightThr=4;
    /*
    for(Int_t j = 0; j < 65; j++) {
    for(Int_t i = 0; i < 24; i++) {
    */
    std::pair<int,int> apair(0,0);
    std::vector<std::pair<int,int> > RowColsAboveThr; RowColsAboveThr.clear();
    if(quarter==0)
      for(unsigned int i=0;i<24;i++)
	for(unsigned int j=0;j<65;j++){
	  if(H0CluPadIdsMatrix[i][j].size()>=StraightThr)
	    {
	      apair.first=i; apair.second=j;
	      if(H0CluPadIdsMatrix[i][j].size()>10)
		{
		  
		  std::cout<<"ERRORRR H0CluPadIdsMatrix["<<i<<"]["<<j<<"].size()>10"<<std::endl;
		  for(unsigned int Kj=0;Kj<H0CluPadIdsMatrix[i][j].size();Kj++){
		    std::cout<<"     |"<<i<<" "<<j<<": ID:"<<H0CluPadIdsMatrix[i][j].at(Kj);
		  }
		  std::cout<<""<<std::endl;
		  
		}

	      RowColsAboveThr.push_back(apair);
	    }
	}
    if(quarter==1)
      for(unsigned int i=0;i<24;i++)
	for(unsigned int j=0;j<65;j++){
	  if(H1CluPadIdsMatrix[i][j].size()>StraightThr)
	    {
	      apair.first=i; apair.second=j;
	      RowColsAboveThr.push_back(apair);
	    }
	}

    if(quarter==2)
      for(unsigned int i=0;i<24;i++)
	for(unsigned int j=0;j<65;j++){
	  if(H2CluPadIdsMatrix[i][j].size()>StraightThr)
	    {
	      apair.first=i; apair.second=j;
	      RowColsAboveThr.push_back(apair);
	    }
	}

    if(quarter==3)
      for(unsigned int i=0;i<24;i++)
	for(unsigned int j=0;j<65;j++){
	  if(H3CluPadIdsMatrix[i][j].size()>StraightThr)
	    {
	      apair.first=i; apair.second=j;
	      RowColsAboveThr.push_back(apair);
	    }
	}
    
    

    
    //choose only one pair: depending on the cls, a  sinngle trk can produce two or more tower of pad above thr. 
    int row1=0;int row2=0; int rowmin=0;// int rowmax=0;
    int col1=0;int col2=0; int colmin=0;// int colmax=0;
    std::vector<std::pair<int,int> > RowColsAboveThrUnique; RowColsAboveThrUnique.clear();

    for(unsigned int i=0;i<RowColsAboveThr.size();i++){

      bool singlecountGuarantee=true;
      
       for(unsigned int j=0;j<RowColsAboveThrUnique.size();j++)
	 {
	   //In order to avoid double count you should also "clusterize" tower above THR

	   if(abs(RowColsAboveThr.at(i).first-RowColsAboveThrUnique.at(j).first)==0)
	     if(abs(RowColsAboveThr.at(i).second-RowColsAboveThrUnique.at(j).second)==1)
	       {
		 singlecountGuarantee=false;
	       }
	   
	   if(abs(RowColsAboveThr.at(i).first-RowColsAboveThrUnique.at(j).first)==1)
	     if(abs(RowColsAboveThr.at(i).second-RowColsAboveThrUnique.at(j).second)==0)
	       {
		 singlecountGuarantee=false;
	       }

	   
	   //Try to clusterize "Triplet or Quadruplet Square Towers", avoiding double count of the same clusters.
	   //Warning: 4 pad clusters with non square shape or blobs can be assigned to more than one towers.
	   if(abs(RowColsAboveThr.at(i).first-RowColsAboveThrUnique.at(j).first)==1)
	     if(abs(RowColsAboveThr.at(i).second-RowColsAboveThrUnique.at(j).second)==1)	       
	       {
		 row1=RowColsAboveThr.at(i).first; col1=RowColsAboveThr.at(i).second;
		 row2=RowColsAboveThrUnique.at(j).first; col2=RowColsAboveThrUnique.at(j).second;
		 rowmin=min(row1,row2);
		 // rowmax=max(row1,row2);
		 colmin=min(col1,col2);
		 //colmax=max(col1,col2);
		 if(quarter==0){
		 if(H0CluPadIdsMatrix[rowmin+1][colmin+1].size()>StraightThr)
		   if(((rowmin+1)!=row1)||((colmin+1)!=col1))
		     if(((rowmin+1)!=row2)||((colmin+1)!=col2))
		       singlecountGuarantee=false;

		 if(H0CluPadIdsMatrix[rowmin+1][colmin].size()>StraightThr)
		   if(((rowmin+1)!=row1)||((colmin)!=col1))
		     if(((rowmin+1)!=row2)||((colmin)!=col2))
		       singlecountGuarantee=false;

		 if(H0CluPadIdsMatrix[rowmin][colmin+1].size()>StraightThr)
		   if(((rowmin)!=row1)||((colmin+1)!=col1))
		     if(((rowmin)!=row2)||((colmin+1)!=col2))
		       singlecountGuarantee=false;
		 
		 if(H0CluPadIdsMatrix[rowmin][colmin].size()>StraightThr)
		   if(((rowmin)!=row1)||((colmin)!=col1))
		     if(((rowmin)!=row2)||((colmin)!=col2))
		       singlecountGuarantee=false;
		 }
		 
		 if(quarter==1){
		 if(H1CluPadIdsMatrix[rowmin+1][colmin+1].size()>StraightThr)
		   if(((rowmin+1)!=row1)||((colmin+1)!=col1))
		     if(((rowmin+1)!=row2)||((colmin+1)!=col2))
		       singlecountGuarantee=false;

		 if(H1CluPadIdsMatrix[rowmin+1][colmin].size()>StraightThr)
		   if(((rowmin+1)!=row1)||((colmin)!=col1))
		     if(((rowmin+1)!=row2)||((colmin)!=col2))
		       singlecountGuarantee=false;

		 if(H1CluPadIdsMatrix[rowmin][colmin+1].size()>StraightThr)
		   if(((rowmin)!=row1)||((colmin+1)!=col1))
		     if(((rowmin)!=row2)||((colmin+1)!=col2))
		       singlecountGuarantee=false;
		 
		 if(H1CluPadIdsMatrix[rowmin][colmin].size()>StraightThr)
		   if(((rowmin)!=row1)||((colmin)!=col1))
		     if(((rowmin)!=row2)||((colmin)!=col2))
		       singlecountGuarantee=false;
		 }

		 if(quarter==2){
		 if(H2CluPadIdsMatrix[rowmin+1][colmin+1].size()>StraightThr)
		   if(((rowmin+1)!=row1)||((colmin+1)!=col1))
		     if(((rowmin+1)!=row2)||((colmin+1)!=col2))
		       singlecountGuarantee=false;

		 if(H2CluPadIdsMatrix[rowmin+1][colmin].size()>StraightThr)
		   if(((rowmin+1)!=row1)||((colmin)!=col1))
		     if(((rowmin+1)!=row2)||((colmin)!=col2))
		       singlecountGuarantee=false;

		 if(H2CluPadIdsMatrix[rowmin][colmin+1].size()>StraightThr)
		   if(((rowmin)!=row1)||((colmin+1)!=col1))
		     if(((rowmin)!=row2)||((colmin+1)!=col2))
		       singlecountGuarantee=false;
		 
		 if(H2CluPadIdsMatrix[rowmin][colmin].size()>StraightThr)
		   if(((rowmin)!=row1)||((colmin)!=col1))
		     if(((rowmin)!=row2)||((colmin)!=col2))
		       singlecountGuarantee=false;
		 }
		 
		 if(quarter==3){
		 if(H3CluPadIdsMatrix[rowmin+1][colmin+1].size()>StraightThr)
		   if(((rowmin+1)!=row1)||((colmin+1)!=col1))
		     if(((rowmin+1)!=row2)||((colmin+1)!=col2))
		       singlecountGuarantee=false;

		 if(H3CluPadIdsMatrix[rowmin+1][colmin].size()>StraightThr)
		   if(((rowmin+1)!=row1)||((colmin)!=col1))
		     if(((rowmin+1)!=row2)||((colmin)!=col2))
		       singlecountGuarantee=false;

		 if(H3CluPadIdsMatrix[rowmin][colmin+1].size()>StraightThr)
		   if(((rowmin)!=row1)||((colmin+1)!=col1))
		     if(((rowmin)!=row2)||((colmin+1)!=col2))
		       singlecountGuarantee=false;
		 
		 if(H3CluPadIdsMatrix[rowmin][colmin].size()>StraightThr)
		   if(((rowmin)!=row1)||((colmin)!=col1))
		     if(((rowmin)!=row2)||((colmin)!=col2))
		       singlecountGuarantee=false;
		 }

		 
	       }
	   
	   //Avoid double count
	   /*
	   if(abs(RowColsAboveThr.at(i).first-RowColsAboveThrUnique.at(j).first)==0)
	     if(abs(RowColsAboveThr.at(i).second-RowColsAboveThrUnique.at(j).second)==0)
	       {
		 singlecountGuarantee=false;
	       }
	   */
	 }

       if(singlecountGuarantee)
	 RowColsAboveThrUnique.push_back(RowColsAboveThr.at(i));
       
       
    }

    

     row2=0; col2=0;std::vector<int> cluIds; 
    for(unsigned int h=0;h<RowColsAboveThrUnique.size();h++)
      {
	row2=RowColsAboveThrUnique.at(h).first;
	col2=RowColsAboveThrUnique.at(h).second;
	cluIds.clear();
	if(quarter==0){
	  cluIds=H0CluPadIdsMatrix[row2][col2];
	}
	if(quarter==1){
	  cluIds=H1CluPadIdsMatrix[row2][col2];
	}
	if(quarter==2){
	  cluIds=H2CluPadIdsMatrix[row2][col2];
	}
	if(quarter==3){
	  cluIds=H3CluPadIdsMatrix[row2][col2];
	}
	//std::cout<<"Calling MakeStraightTubeForThisclusters "<<cluIds.size()<<std::endl;
	Tube atube=MakeStraightTubeForThisclusters(cluIds,T2PadClusterCollectionPtr);	
	
	WinningTubes.push_back(atube);
      }


    
  }


Tube T2RoadPadFinder::MakeStraightTubeForThisclusters(std::vector<int> &cluIds,const T2PadClusterCollection*  T2PadClusterCollectionPtr){
  
  Tube atube;
  // //std::cout<<"-----E"<<std::endl;
  /*
  int minplaneind;
  int maxplaneind;
  */
  double cluX=0.; double cluY=0.; double cluR=0.; double cluPhi=0.; double cluEr=0.; double cluEphi=0.;
  double cluEx=0.; double cluEy=0.;
  double cluZ=0.;

  //unsigned int loopSymb=0;
  //double xDist=0.;
  //double yDist=0.;
  T2Cluster clusterchecked;
  //GetCluster.
  
  int cluuniqueid=0;
  T2DetId T2Det1; 
  bool proceed1=false;
  std::pair <uint32_t,unsigned int> T2Det_CluIndex1;
  unsigned int clusterPlaneSymb=0;
  int minplane=41; int minplanepos=0;
  int maxplane=-1; int maxplanepos=0;

  // std::cout<<"Here"<<std::endl;
  if(cluIds.size()<2)
    std::cout<<"ERROR IN MakeStraightTubeForThisclusters, cluIds.size()<2"<<std::endl;

  for(unsigned int ll=0;ll<cluIds.size();ll++){

    

    cluuniqueid=cluIds.at(ll);

    
    clusterPlaneSymb=Map_ClustuniqueId_PlaneSymb[cluuniqueid];
    // std::cout<<"Index "<<ll<<"Symb: "<<cluuniqueid<<" Plane:"<<clusterPlaneSymb<<std::endl;
    if((int)clusterPlaneSymb>maxplane){
      maxplanepos=ll;
      maxplane=clusterPlaneSymb;
    }
    if((int)clusterPlaneSymb<minplane){
      minplanepos=ll;
      minplane=clusterPlaneSymb;
    }

    
    T2Det_CluIndex1=Map_ClustuniqueId_CluIndexPos[cluuniqueid];
    
    T2Det1=T2DetId((T2Det_CluIndex1.first));
    //std::cout<<"BB"<<std::endl;
    T2PadClusterCollection::const_iterator AllClusPoint1 = T2PadClusterCollectionPtr->find(T2Det1);
    //std::cout<<"---Bo   "<<std::endl;
    proceed1=false;
    
    if(AllClusPoint1!=T2PadClusterCollectionPtr->end()){
      if(AllClusPoint1->second.size()>(T2Det_CluIndex1.second)){       
	proceed1=true;
      }else{
	std::cout<<"Error in T2RoadPadFinder::MakeStraightTubeForThisclusters Ai"<<std::endl;
	return atube;
      }
    }else{
      std::cout<<"Error in T2RoadPadFinder::MakeStraightTubeForThisclusters Bi"<<std::endl;
      return atube;
    }
    //std::cout<<"---Bo2   "<<std::endl;

    if(proceed1){
     
      clusterchecked=AllClusPoint1->second.at((T2Det_CluIndex1.second));
   
      cluR=clusterchecked.GetClusterR();
      cluPhi=clusterchecked.GetClusterPhi()*3.14159265/180.0;
      cluEr=clusterchecked.GetClusterDR();
      cluEphi=clusterchecked.GetClusterDPhi()*3.14159265/180.0;	
      cluX=cluR*cos(cluPhi); cluY=cluR*sin(cluPhi);     
      cluEx=sqrt(cluEr*cluEr*cos(cluPhi)*cos(cluPhi) + cluEphi*cluEphi*cluR*cluR*sin(cluPhi)*sin(cluPhi));
      cluEy=sqrt(cluEr*cluEr*sin(cluPhi)*sin(cluPhi) + cluEphi*cluEphi*cluR*cluR*cos(cluPhi)*cos(cluPhi));  
     
      cluZ=Map_ClustuniqueId_ZPos[cluuniqueid];

      atube.ClustV.push_back(clusterchecked);
      atube.uniqueCluIdV.push_back(cluuniqueid);

      Map_ClustuniqueId_CluStatus[cluuniqueid]=1;

      atube.VectR.push_back(cluR);
      atube.VectPhi.push_back(cluPhi);
      atube.VectX.push_back(cluX);
      atube.VectY.push_back(cluY);
      atube.VectZ.push_back(cluZ);
      atube.VectEX.push_back(cluEx);
      atube.VectEY.push_back(cluEy);  
     //std::cout<<"---Bo3   "<<std::endl;
    }
    //std::cout<<"--ll"<<std::endl;
  }//for cluIds end 
  //   std::cout<<"Filling tube.."<<minplanepos<<" "<<maxplanepos<<std::endl;
  atube.Xi=atube.VectX.at(minplanepos);
  atube.Yi=atube.VectY.at(minplanepos);
  atube.Zi=atube.VectZ.at(minplanepos);
  atube.Ri=atube.VectR.at(minplanepos);
  atube.Phii=atube.VectPhi.at(minplanepos);
  atube.planei=minplane;
  
  atube.Xf=atube.VectX.at(maxplanepos);
  atube.Yf=atube.VectY.at(maxplanepos);
  atube.Zf=atube.VectZ.at(maxplanepos);
  atube.Rf=atube.VectR.at(maxplanepos);
  atube.Phif=atube.VectPhi.at(maxplanepos);
  atube.planef=maxplane;
  
  TubeGeneratorId++;
  atube.uniquetubeId=TubeGeneratorId; 
  if((maxplane>39)||(minplane<0))
    std::cout<<"Error in MakeStraightTubeForThisclusters. Wrong minplane-maxplane"<<std::endl;

  //double Ax,Ay,Xi,Yi,Xf,Yf,Zi,Zf,Ri,Rf,Phii,Phif,planei,planef;
  atube.Ax=(atube.Xi-atube.Xf)/(atube.Zi-atube.Zf);
  atube.Ay=(atube.Yi-atube.Yf)/(atube.Zi-atube.Zf);  
  atube.seedId=cluIds.at(minplanepos);
  

  return atube;
}







void T2RoadPadFinder::produce(edm::Event& event, const edm::EventSetup& setup) {


  //  std::cout<<"Roadproduce start"<<std::endl;
  
  
  if(verbosity>=1)
    std::cout<<"*****************************************************Produce Begin*****************************************"<<std::endl;


  auto_ptr<T2RoadCollection> theT2Roads (new T2RoadCollection());  //Container for all the event PAD-roads
  
  auto_ptr<T2RoadCollection> theT2Roads_withStrips(new T2RoadCollection());  //Container for all the event Standard-roads
  
  auto_ptr<T2RoadCollection> theT2Roads_withStripsOverlapFree(new T2RoadCollection());  //Container for all the event Standard-roads
 
  auto_ptr<T2Hit_to_Track_Map> theT2Hit_to_Track_Map(new T2Hit_to_Track_Map());
  
  //  std::cout<<"A"<<std::endl;
  //T2PadClusterCollectionPtr= new T2PadClusterCollection();
  //event.getByLabel("T2MCl","T2PadClusters",t2padclcoll);
  event.getByLabel(/*"T2MCl"*/CluLabel,"T2PadClusters",t2padclcoll);
  

  //Static Matrix gave me a lot of problems so I will use this std::vector.
  // std::vector< std::vector< std::vector<int> > > H0CluPadIdsMatrix;
  H0CluPadIdsMatrix.clear();
  H1CluPadIdsMatrix.clear();
  H2CluPadIdsMatrix.clear();
  H3CluPadIdsMatrix.clear();
  for(unsigned int i=0;i<24;i++)
    {
      std::vector< std::vector<int> > initvectColsforARow; 
      for(unsigned int j=0;j<65;j++){
	std::vector<int> initvectInner; 
	initvectInner.clear();
	initvectColsforARow.push_back(initvectInner);
      }
      H0CluPadIdsMatrix.push_back(initvectColsforARow);
      H1CluPadIdsMatrix.push_back(initvectColsforARow);
      H2CluPadIdsMatrix.push_back(initvectColsforARow);
      H3CluPadIdsMatrix.push_back(initvectColsforARow);
      
    }

  //std::cout<<"Matrix size: "<<H0CluPadIdsMatrix.size()<<"-"<<H0CluPadIdsMatrix.at(0).size()<<std::endl;
  
  // std::cout<<"A. Verbosity is "<<verbosity<<std::endl;
  InitializeMapClusterAssociators(event);
  const T2PadClusterCollection *T2PadClusterCollectionPtr = t2padclcoll.product();

 
  // event.getByLabel("T2MCl","T2StripClusters",t2strclcoll);
  event.getByLabel(CluLabel,"T2StripClusters",t2strclcoll);
  t2strclcoll.product();


  event.getByLabel(HitLabel,"T2Hits",t2hitcoll);
  const T2HitCollection *T2HitCollectionPtr = t2hitcoll.product();
  
  event.getByLabel(HitLabel,"T2HitsMap",t2hitcollmap);
  const T2HitCollectionMapping *T2HitCollectionMapPtr = t2hitcollmap.product();

  event.getByLabel(HitLabel,"T2HitsMapPadonStrip",t2hitpadonstripMap);
  const T2PadStripAssociator *T2HitsMapPadonStripPtr = t2hitpadonstripMap.product();

  //Load a maps containing <ClustuniqueId->CluRoadId> and <ClustuniqueId->CluStatus>  
  //The map will be filled for the quarters eventually included in QuartedSelected
  //In the meantime also a single map of <ClustuniqueId->index position inside the t2padclcoll will filled>
  
  if(verbosity>=1)
    std::cout<<"Maps Initialized with Size: "<<Map_ClustuniqueId_CluStatus.size()<<" clusters"<<std::endl;


  std::vector<Tube> AllInitialTubes; 
  std::vector<Tube> AllFinalTubes;

  WinningTubes.clear();
  AllInitialTubes.clear();
  AllFinalTubes.clear();

  //std::cout<<"---C "<<std::endl;
  //In High multiplicity events you can try to look to straight towers of pads before the standard finding.
  if(useStraightPadTowers){
    if((occupadcluH0/10.)>NumPadCluOccupancyAlert)
      StraightTubeFinder(0,WinningTubes,T2PadClusterCollectionPtr);

    if((occupadcluH1/10.)>NumPadCluOccupancyAlert)
      StraightTubeFinder(1,WinningTubes,T2PadClusterCollectionPtr);

    if((occupadcluH2/10.)>NumPadCluOccupancyAlert)
      StraightTubeFinder(2,WinningTubes,T2PadClusterCollectionPtr);

    if((occupadcluH3/10.)>NumPadCluOccupancyAlert)
      StraightTubeFinder(3,WinningTubes,T2PadClusterCollectionPtr);
  }
  //std::cout<<"A"<<std::endl;

  ////std::cout<<"---D "<<std::endl;
  if(verbosity>=3.)
    for(unsigned int a=0;a<WinningTubes.size();a++)
      {
	if(verbosity>=3.)
	  std::cout<<"Straight tube:"<<WinningTubes.at(a).Xi<<" Yi:"<<WinningTubes.at(a).Yi<<" size:"<<WinningTubes.at(a).uniqueCluIdV.size()<<" ID:"<<WinningTubes.at(a).uniquetubeId<<" r0-c0:  ";
	for(unsigned int as=0;as<WinningTubes.at(a).ClustV.size();as++)
	  std::cout<<WinningTubes.at(a).ClustV.at(as).GetEntries()[0].rad_coord<<"-"<<WinningTubes.at(a).ClustV.at(as).GetEntries()[0].ang_coord<<"  R:"<<WinningTubes.at(a).ClustV.at(as).GetClusterR()<<" Phi:"<<WinningTubes.at(a).ClustV.at(as).GetClusterPhi()<<" CLSize:"<<WinningTubes.at(a).ClustV.at(as).GetNoOfEntries()<<" | ";
	std::cout<<std::endl;     
      }
  


  //  std::vector<Tube> FinalTubesFromAllQuarters;
   BiggestTubeAngleConsidered=BiggestTubeAngleConsideredH0;
   
   if(verbosity>=1)
     {
       std::cout<<"*******Quarter 0 processing********* with BiggestTubeAngleConsidered="<<BiggestTubeAngleConsidered<<std::endl;    
     }
 
  
  for(unsigned int seedplanelevel=0; seedplanelevel<10; seedplanelevel++){ //Piano da dove parte la costr del seed 
   
    MakethreeWithThisPlaneSeed(seedplanelevel,AllInitialTubes,AllFinalTubes,T2PadClusterCollectionPtr); //fromseedplane to forward.
    AllInitialTubes.swap(AllFinalTubes);
    AllFinalTubes.clear();  //The updated final tube will stay in AllInitialTubes
  }
 
  //Saving also the tube arriving at the end (some winning was already stored)
  for(unsigned int ft=0; ft<AllInitialTubes.size(); ft++){  
    if((int)AllInitialTubes.at(ft).uniqueCluIdV.size()>=Nmin_padsFinal)
      {
	WinningTubes.push_back(AllInitialTubes.at(ft));
      }
  }

  /*
std::cout<<"second part...."<<std::endl;
if(true){
    std::cout<<"second part on"<<std::endl;
  */



  AllInitialTubes.clear();
  BiggestTubeAngleConsidered=BiggestTubeAngleConsideredH1;
 
  if(verbosity>=1)
    {
      std::cout<<"*******Quarter 1 processing********* with BiggestTubeAngleConsidered="<<BiggestTubeAngleConsidered<<std::endl;    
    }

  
  for(unsigned int seedplanelevel=10; seedplanelevel<20; seedplanelevel++){ //Piano da dove parte la costr del seed    
    MakethreeWithThisPlaneSeed(seedplanelevel,AllInitialTubes,AllFinalTubes,T2PadClusterCollectionPtr); //fromseedplane to forward.
    AllInitialTubes.swap(AllFinalTubes);
    AllFinalTubes.clear();
  }
  
  //Saving also the tube arriving at the end (some winning was already stored)
  for(unsigned int ft=0; ft<AllInitialTubes.size(); ft++){  
    if((int)AllInitialTubes.at(ft).uniqueCluIdV.size()>=Nmin_padsFinal)
      {
	WinningTubes.push_back(AllInitialTubes.at(ft));
      }
  }
  //std::cout<<"Q1SIZE-"<<WinningTubes.size()<<std::endl;
  

  
 
 BiggestTubeAngleConsidered=  BiggestTubeAngleConsideredH2; 
 if(verbosity>=1)
    {
      std::cout<<"*******Quarter 2 processing********* with BiggestTubeAngleConsidered="<<BiggestTubeAngleConsidered<<std::endl;    
    }

  AllInitialTubes.clear();
  for(unsigned int seedplanelevel=20; seedplanelevel<30; seedplanelevel++){ //Piano da dove parte la costr del seed 
   
    MakethreeWithThisPlaneSeed(seedplanelevel,AllInitialTubes,AllFinalTubes,T2PadClusterCollectionPtr); //fromseedplane to forward.
    AllInitialTubes.swap(AllFinalTubes);
    AllFinalTubes.clear();
  }
 
  //Saving also the tube arriving at the end (some winning was already stored)
  for(unsigned int ft=0; ft<AllInitialTubes.size(); ft++){ 
    
    if((int)AllInitialTubes.at(ft).uniqueCluIdV.size()>=Nmin_padsFinal)
      {
	WinningTubes.push_back(AllInitialTubes.at(ft));
      }
  }

  //std::cout<<"Q2SIZE->"<<WinningTubes.size()<<std::endl;
  

  
  

  BiggestTubeAngleConsidered=  BiggestTubeAngleConsideredH3;
  if(verbosity>=1)
    {
      std::cout<<"*******Quarter 3 processing********* with BiggestTubeAngleConsidered="<<BiggestTubeAngleConsidered<<std::endl;    
    }

  AllInitialTubes.clear();
  for(unsigned int seedplanelevel=30; seedplanelevel<40; seedplanelevel++){ //Piano da dove parte la costr del seed 
   
    MakethreeWithThisPlaneSeed(seedplanelevel,AllInitialTubes,AllFinalTubes,T2PadClusterCollectionPtr); //fromseedplane to forward.
    AllInitialTubes.swap(AllFinalTubes);
    AllFinalTubes.clear();
  }

  //Saving also the tube arriving at the end (some winning was already stored)
  for(unsigned int ft=0; ft<AllInitialTubes.size(); ft++){  
    if((int)AllInitialTubes.at(ft).uniqueCluIdV.size()>=Nmin_padsFinal)
      {
	WinningTubes.push_back(AllInitialTubes.at(ft));
      }
  }
  //std::cout<<"Q3SIZE->"<<WinningTubes.size()<<std::endl;
  

if(verbosity>=1)
  std::cout<<"Tube finding complete!!! Tubes founds: "<<WinningTubes.size()<<"writing a first T2Road(Pad) collection"<<std::endl;



  //#ifndef DEBUG_T2Road  
  //Final Check on AllInitialTubes to put in WinningTubes;
  for(unsigned int ft=0; ft<WinningTubes.size(); ft++){  
    
    T2Road theRoad; 
    unsigned int sizeroad=WinningTubes.at(ft).ClustV.size();
    // std::cout<<" saving a road with size "<<sizeroad<<std::endl;
    theRoad.RoadID=WinningTubes.at(ft).uniquetubeId;
    
     for(unsigned int k=0; k<sizeroad;k++)
       {
	 //Eventual Info for smearing studies.
	 Map_ClustuniqueId_VectorCluRoadId[WinningTubes.at(ft).uniqueCluIdV.at(k)].push_back(WinningTubes.at(ft).uniquetubeId);
	 
	 //	 std::cout<<"         | "<<WinningTubes.at(ft).VectX.at(k)<<" "<<WinningTubes.at(ft).VectY.at(k)<<" "<<WinningTubes.at(ft).ClustV.at(k).GetClusterR()<<"-"<<WinningTubes.at(ft).ClustV.at(k).GetClusterPhi()<<" id:"<<WinningTubes.at(ft).uniquetubeId;
	 theRoad.thisPadRoad.push_back(WinningTubes.at(ft).ClustV.at(k));

	 
       }
     //std::cout<<std::endl;
     //std::cout<<"Road created boh "<<theRoad.thisPadRoad.size()<<"..."<<std::endl; 
   // swap(theRoad.thisPadRoad,WinningTubes.at(ft).ClustV);
     /*theRoad.SetRoadID(AllInitialTubes.at(ft).uniquetubeId);*/
     //std::cout<<"Road created forse "<<theRoad.thisPadRoad.size()<<std::endl;
   theT2Roads->push_back(theRoad);	
  
   
  }
  if(verbosity>1)
    std::cout<<"All events Pad-Road saved. #Roads: "<<theT2Roads->size()<<std::endl; 

  //
  //
  //Redefine PadTubes-Roadnumbering using winning tubes and redefine the map pad association.
  //Now combine the strip. A criteria to eliminate possible repetion could be to see how many track 
  //
  // cout<<"End Road2 producer"<<endl;
  // std::cout<<"Road C"<<std::endl;
  //Maybe otput should be two road collection named T2RoadCollPad and T2RoadColl, the last one containing the strip->Hit association



 
  
  //------------------------------------------------------------------------------------------------------------------------------
  //START Now To CREATE CL1-HIT ROADS, combining the strips:
  //Assumed that close roads are already resolved without strip helps. Eventually overlapping road should be threated 
  //before using already the strip information. Now I suppose that eventually I have already used that information to keep only
  //unambiguous roads set of Pad road that can have one or more tracks. 
  //------------------------------------------------------------------------------------------------------------------------------
  
  int padsymbidusedbefore=0; double planeZ=0.; long int realHardwCluId=0;
  /*
  std::vector<std::vector <long> > allROADcombinationsID;
  std::vector<std::vector <long> >  allROADcombinationsID_pad;
  std::vector<std::vector <long> >  allROADcombinationsID_strip;
  std::vector<int> vectStripInd;

  
  std::map<double, std::vector<long> > mapZvsUniqueHitId;

  std::map<double, std::vector<long> > mapZvsUniquePadCluId;
  std::map<double, std::vector<long> > mapZvsUniqueStripCluId;
  */
  

  std::vector<long int> vectStripInd;
  std::vector<std::vector <unsigned long int> >  allROADcombinationsID_pad;
  std::vector<std::vector <unsigned long int> >  allROADcombinationsID_strip;

  std::map<double, std::vector<unsigned long int> > mapZvsUniquePadCluId;
  std::map<double, std::vector<unsigned long int> > mapZvsUniqueStripCluId;


  

  std::vector<T2Hit> CL1RoadHits;
  std::vector<std::vector<T2Hit> > ALLCL1RoadHits;  //Since one Pad road can have more cl1 Hit roads.
  std::vector<std::vector<T2Hit> > SuccessFull_CL1RoadHits;  //Since one Pad road can have more cl1 Hit roads.
  std::vector<std::vector<T2Hit> > RoadInOverlapRegion;
  std::vector<std::vector<T2Hit> > OverlapRoadMerged;
  OverlapRoadMerged.clear(); RoadInOverlapRegion.clear(); 
  std::vector<int> RoadInOverlapRegion_index; RoadInOverlapRegion_index.clear();

  std::vector<std::vector<T2Hit> > CollectionOFSuccessFull_CL1RoadHits; //All successfull CL1Hit road *of the event*;
  // BxByAxAy_AllRoads
  std::vector<std::vector<double> > BxByAxAy_AllRoads;
  std::vector<std::vector<double> > SuccessFull_BxByAxAy_Roads;

  std::vector<std::vector<double> > StripJoined_BxByAxAy;
  std::vector<std::vector<T2Hit> > StripJoinedRoadHits;
  std::vector<std::vector<T2Hit> > ParallelBranches;

  CollectionOFSuccessFull_CL1RoadHits.clear();
  
  std::map<int, unsigned int>::const_iterator ConcurrentmapIt; 

  int aroadID=0;
for(T2RoadCollection::iterator itroad = theT2Roads->begin(); itroad != theT2Roads->end(); itroad++)
  // if(itroad->RoadID==4)
  {

    aroadID++;
    if(verbosity>=1)
      std::cout<<" Working on road id "<<itroad->RoadID<<" and size: "<<itroad->thisPadRoad.size()<<std::endl;
    
    //This data structure contain information about candidate trk-road extracted from *ONE ROAD PAD*
    CL1RoadHits.clear(); 
    /* mapZvsUniqueHitId.clear();*/  
    mapZvsUniquePadCluId.clear();   mapZvsUniqueStripCluId.clear();
    StripJoined_BxByAxAy.clear(); 
    StripJoinedRoadHits.clear();
    
    
    
    ALLCL1RoadHits.clear(); BxByAxAy_AllRoads.clear();    
    SuccessFull_CL1RoadHits.clear(); SuccessFull_BxByAxAy_Roads.clear();
    
    
     //a)Match the pad of the road with the class1 Hits. Need a function that takes the Padroad and come back with a matrix of 
      //CL1 Hit-ids associaed where same row mean same Z Hits. But Hit-ids are of the form <int,int>. The pair should be compacted. 
      //Since each hit index is less than 150000 one can take the rule=UniqueHitIndx= padIndx+150000*stripindex. 
      //The long data type should be enogh to contain everything. 
    // PAD HARD ID IS CALCULATED AS uniquePadIDinEvt=symbplane*1560+col_beg*24+row_beg;

    std::vector<unsigned int> PadonlyHit_vctInd; 
    PadonlyHit_vctInd.clear();
    std::vector<T2Hit>  PadonlyHit;
    PadonlyHit.clear();
    
    // thisPadRoad is a std::vector<T2Cluster>
    for(unsigned int c=0;c<itroad->thisPadRoad.size();c++)
      {
	//takeClupadId
	realHardwCluId=itroad->thisPadRoad.at(c).ClustId_unique;
	padsymbidusedbefore=Map_RealClusterID_ToSymbUsed[realHardwCluId];
	//if(verbosity>=1){
	//     std::cout<<"CLU Pad HrdwID: "<<realHardwCluId<<"corresponding to symb"<<padsymbidusedbefore<<std::endl;
	//     std::cout<<"HRD-plane:"<<realHardwCluId/1560<<" col="<<(realHardwCluId%1560)/24<<std::endl;
	//}
	planeZ=Map_ClustuniqueId_ZPos[padsymbidusedbefore];	
	vectStripInd.clear();
	

	 T2PadStripAssociator::const_iterator padonstripIDSptr = T2HitsMapPadonStripPtr->find(realHardwCluId);
	 if(padonstripIDSptr!=T2HitsMapPadonStripPtr->end()){
	   
	   //	   std::cout<<"Pad id:"<<padsymbidusedbefore<<" is a class 1 hit"<<std::endl;
	   
	   vectStripInd = padonstripIDSptr->second;

	 }else{
	   //std::cout<<"Pad id:"<<padsymbidusedbefore<<" is a pad only"<<std::endl;
	   //This is a pad only Hit  
	   PadonlyHit_vctInd.push_back(c);
	    
	   //Lets also check if is a class 1 Hit. In this case there is an error otherwise not. I'm looking if this pad is a Only-ClusterHit.
	   //if(verbosity>=1)
	   std::pair<long int, long int> hitidtest (realHardwCluId,-1);
	   //  std::pair<int,int> hitidtest (realHardwCluId,-1);
	       T2HitCollectionMapping::const_iterator iterT2Map = T2HitCollectionMapPtr->find(hitidtest);
	       if(iterT2Map!=T2HitCollectionMapPtr->end())
		 {
		   int position=iterT2Map->second;
		   if(verbosity>=1) 
		     std::cout<<padsymbidusedbefore<<"Hit is a pad Only with position in the collection="<<position<<std::endl;
		   
		   T2HitCollection::const_iterator iterCOLL = T2HitCollectionPtr->begin()+position;

		   //iterCOLL->Set_RoadID(aroadID);

		   if(position<(T2HitCollectionPtr->end() - T2HitCollectionPtr->begin()))
		     {
		       PadonlyHit.push_back((*iterCOLL));
		       if(verbosity>=1)
			 std::cout<<"Valid index: Num Pad-Strip="<<iterCOLL->GetHitNumPad()<<"-"<<iterCOLL->GetHitNumStrip()<<std::endl;
		     }
		   else
		     {
		       
		       std::cout<<"Error in T2RoadPadFinder:: T2HitCollection index out of range"<<std::endl;
		     }
		   
		 }
	       else
		 {
		   std::cout<<"Error in T2RoadPadFinder::Cannot find the Hit Id associated to the cluster pad:"<<std::endl;
		   std::cout<<"Hardware-CluID:"<<realHardwCluId<<" R-Phi:"<<itroad->thisPadRoad.at(c).GetClusterR()<<"-"<<itroad->thisPadRoad.at(c).GetClusterPhi()<<" Symb of the cluster in the tracker:"<<padsymbidusedbefore<<" Size-Type: "<<itroad->thisPadRoad.at(c).GetNoOfEntries()<<" "/*<<itroad->thisPadRoad.at(c).GetClusterType*/<<std::endl;
		 }
	 }
	 
	 //PROCEDURE WRONG!!!!!
	 //	 std::cout<<"Strip size associated to pad "<<realHardwCluId<<": "<<vectStripInd.size()<<std::endl;
	 /*
	 for(unsigned int stripcounter=0; stripcounter<vectStripInd.size();stripcounter++)
	   {
	     int stripHrdwId=vectStripInd.at(stripcounter);
	     unsigned long long HitUniqueId=realHardwCluId+(unsigned long long)(150000*stripHrdwId);
	     
	     unsigned long long TEST_HitUniqueId_PAD=(unsigned long long)(HitUniqueId%150000);
	     unsigned long long TEST_HitUniqueId_STRIP=(unsigned long long)(HitUniqueId/150000);
	    
	     mapZvsUniqueHitId[planeZ].push_back(HitUniqueId);
	     mapZvsUniqueStripCluId[planeZ].push_back(stripHrdwId);
	     mapZvsUniquePadCluId[planeZ].push_back(realHardwCluId);
	   }
	 */
	  for(unsigned int stripcounter=0; stripcounter<vectStripInd.size();stripcounter++)
	   {
	     
	     unsigned long int stripHrdwId=static_cast<unsigned long int>(vectStripInd.at(stripcounter));	     	   
	     unsigned long int padHrdwId=static_cast<unsigned long int>(realHardwCluId);

	     mapZvsUniqueStripCluId[planeZ].push_back(stripHrdwId);
	     mapZvsUniquePadCluId[planeZ].push_back(padHrdwId);
	   }
    
       }

     
    if(verbosity>=1)
      std::cout<<"Number of *different* Z CL1 Hits in combinations for this road (id:"<< itroad->RoadID<<"):"<<mapZvsUniqueStripCluId.size()<<std::endl;

    // if(mapZvsUniqueHitId.size() >= MinimumNumCl1Hit)
    if((int)mapZvsUniqueStripCluId.size() >= MinimumNumCl1Hit)
      {

	//b)Once you have constructed the matrix you can use the T2trackCollection methods for the combiantion of different
	//number representing the CL1 Hit.

	//allROADcombinationsID.clear(); 
	//allROADcombinationsID= RoadPadCL1HitsCombination(mapZvsUniqueHitId);//Each vector is a candidate

	allROADcombinationsID_pad.clear(); 
	allROADcombinationsID_pad= RoadPadCL1HitsCombination(mapZvsUniquePadCluId);//Each vector is a candidate
	allROADcombinationsID_strip.clear(); 
	allROADcombinationsID_strip= RoadPadCL1HitsCombination(mapZvsUniqueStripCluId);//Each vector is a candidate
    /*
     //TEST COMBINATION
      mapZvsUniqueHitId.clear();
      std::vector<long> row; row.push_back(1);row.push_back(2);row.push_back(3);
      mapZvsUniqueHitId[1.0]=row;
      row.clear();
      row.push_back(4);
      mapZvsUniqueHitId[2.0]=row;
      row.clear();
      row.push_back(5); row.push_back(6);
      mapZvsUniqueHitId[3.0]=row;
      */

    //Construct Now the real vectorS of Cl1 Hits candidateS
    // For each combination decompose hituniqueId in a pair <padid,stripId>. Then found the corresponding Cl1 T2Hit in the T2CollectionMapping and saved it in a ALLCL1RoadHits.
	
	/*
	int HitUniqueId=0;
	int HitUniqueId_STRIP=0;
	int HitUniqueId_PAD=0;
	*/

	long int HitUniqueId_STRIP=0;
	long int HitUniqueId_PAD=0;
	//=realHardwCluId+150000*stripHrdwId; 
	if(verbosity>0){
	  std::cout<<"Num CL1 Hits in A combinations for this road "<< itroad->RoadID<<":"<<mapZvsUniqueStripCluId.size()<<std::endl;
	  std::cout<<"# Combinazioni: (idroad:"<< itroad->RoadID<<") combination:"<<allROADcombinationsID_pad.size()<<"-"<<allROADcombinationsID_strip.size()<<std::endl;
	}

	for(unsigned int m1=0;m1<allROADcombinationsID_strip.size();m1++){
	  //std::cout<<"    Comb #"<<m1<<":  "<<std::endl;
	  CL1RoadHits.clear();
	  for(unsigned int m2=0;m2<allROADcombinationsID_strip.at(m1).size();m2++){
	    
	    /*
	    HitUniqueId=allROADcombinationsID.at(m1).at(m2);
	    HitUniqueId_PAD=allROADcombinationsID_pad.at(m1).at(m2);// (int)(HitUniqueId%150000);
	    HitUniqueId_STRIP=allROADcombinationsID_strip.at(m1).at(m2);//(int)(HitUniqueId/150000);	    
	    std::pair<int,int> hitidtest (HitUniqueId_PAD,HitUniqueId_STRIP);
	    */

	    HitUniqueId_PAD=static_cast<long int>(allROADcombinationsID_pad.at(m1).at(m2));// (int)(HitUniqueId%150000);
	    HitUniqueId_STRIP=static_cast<long int>(allROADcombinationsID_strip.at(m1).at(m2));//(int)(HitUniqueId/150000);
	    
	    std::pair<long int,long int> hitidtest (HitUniqueId_PAD,HitUniqueId_STRIP);


	    T2HitCollectionMapping::const_iterator iterT2Map = T2HitCollectionMapPtr->find(hitidtest);
	    if(iterT2Map!=T2HitCollectionMapPtr->end())
	      {
		int position=iterT2Map->second;	   
		T2HitCollection::const_iterator iterCOLL = T2HitCollectionPtr->begin()+position;	    
		if(position<(T2HitCollectionPtr->end() - T2HitCollectionPtr->begin()))
		  {
		    
		    CL1RoadHits.push_back((*iterCOLL));	     
		    
		    //std::cout<<"  |Class:"<<(*iterCOLL).GetHitClass()<<" R:"<<(*iterCOLL).GetHitR()<<" Phi:"<<(*iterCOLL).GetHitPhi()<<" NumStrip:"<<(*iterCOLL).GetHitNumStrip()<<" Z:"<<(*iterCOLL).GetHitZ()<<std::endl;
		  }
		else
		  std::cout<<"ERROR on position"<<std::endl;
	      }
	    else
	      std::cout<<"Hit"<<HitUniqueId_PAD<<"-"<<HitUniqueId_STRIP<<"not found"<<std::endl;
	    
	  }
	  
	  //std::cout<<"-Other comb-"<<std::endl;
	  ALLCL1RoadHits.push_back(CL1RoadHits);
	}
	

	
	if(verbosity>=1){

	  std::cout<<"Num Uncleaned Roads from all the cl1Hit combiantions of the Pad ROAD:"<<ALLCL1RoadHits.size()<<std::endl;

	  for(unsigned int y=0;y<ALLCL1RoadHits.size();y++)
	    {
	      std::cout<<std::endl;
	      std::cout<<"Cl1 raw road R-Z:"<<std::endl;
	      for(unsigned int yy=0;yy<ALLCL1RoadHits.at(y).size();yy++)
		std::cout<<ALLCL1RoadHits.at(y).at(yy).GetHitR()<<"  "<<ALLCL1RoadHits.at(y).at(yy).GetHitZ()<<" | ";
	      
	    }
	}

	//c) perform track selection from one PadRoad+CL1Hit to several Pad+ strip roads.
	// I want to change strategy respect to TrackProducer2.
	//---> Takes best combination under threshold -> Remove it. If nothing is under chi threshold? Drop the worst outlier from each combination and reevaluate the best allows removing until 3 CL1 hits are there. Than drop the combiantions.
	//idea per il futuro: uno puo' anche cercare di minimizzare un chi2 globale, che massimizzi il mulitracking.
	
	bool allCandidateExtracted=false;
	double chi2XProb=0;
	double chi2YProb=0;

	
	bool Atleast_aHitRoad_Created=false;

	while(allCandidateExtracted==false)
	  {
	    
	    
	    // for(unsigned int jj=0;jj<ALLCL1RoadHits.size();jj++)
	    //{
	    //std::cout<<"RESTARTinwhile: Road with size "<<ALLCL1RoadHits.at(jj).size()<<" Left"<<std::endl;		
	    //}
	    //SAFEcounter++;
	    
	    /*
	      for(unsigned int jj=0;jj<ALLCL1RoadHits.size();jj++)
	      {
		std::cout<<"RESTARTinwhile: Road with size "<<ALLCL1RoadHits.at(jj).size()<<" Left"<<std::endl;		
	      }
	    	std::vector<unsigned int> IdtoremoveAgain2; IdtoremoveAgain2.clear();
		unsigned int alive=0;
		for(unsigned int jj=0;jj<ALLCL1RoadHits.size();jj++)
		  {
		    if(ALLCL1RoadHits.at(jj).size()>=MinimumNumCl1Hit)
		      {alive++;}
		    else
		      IdtoremoveAgain2.push_back(jj); 
		  }

		for(unsigned int m=0;m<IdtoremoveAgain2.size();m++)
		  ALLCL1RoadHits.erase(ALLCL1RoadHits.begin()+IdtoremoveAgain2.at(m));
	    */
		
	   
	    //Note: GiveMeBestCandidate clear and recalculate: BxByAxAy_AllRoads
	    
	    
	    int bestroadIndex=GiveMeBestCandidate(ALLCL1RoadHits, chi2XProb, chi2YProb,BxByAxAy_AllRoads);
	    if(verbosity>=1){

	      std::cout<<"BestRoadIndex: "<<bestroadIndex<<" chi2XProb:"<<chi2XProb<<"   best road size:"<<ALLCL1RoadHits.at(bestroadIndex).size()<<"X-Y-Z:"<<std::endl;
	      
	      for(unsigned int iu=0;iu<ALLCL1RoadHits.at(bestroadIndex).size();iu++)
		std::cout<<" X:"<<ALLCL1RoadHits.at(bestroadIndex).at(iu).GetHitX()<<" | Y:"<<ALLCL1RoadHits.at(bestroadIndex).at(iu).GetHitY()<<" | Z:"<<ALLCL1RoadHits.at(bestroadIndex).at(iu).GetHitZ()<<"            ";
	      std::cout<<std::endl;
	    }


	    if((chi2XProb>=chi2XYProb_Thr))//actually a global chi-prob is calculated
	      {
	
		//Promoved the cl1 road as final and removed from ALLCL1RoadHits
		SuccessFull_CL1RoadHits.push_back(ALLCL1RoadHits.at(bestroadIndex));
		SuccessFull_BxByAxAy_Roads.push_back(BxByAxAy_AllRoads.at(bestroadIndex));
		
		//std::cout<<"OK BestRoadInd: "<<bestroadIndex<<" ProbX:"<<chi2XProb<<" CollSize:"<<SuccessFull_CL1RoadHits.size()<<" Trk size:"<<ALLCL1RoadHits.at(bestroadIndex).size()<<std::endl;
		
		std::vector<T2Hit> winningRoad=ALLCL1RoadHits.at(bestroadIndex);
		ALLCL1RoadHits.erase(ALLCL1RoadHits.begin()+bestroadIndex);

		RemoveHitsFromLoserCombinations(ALLCL1RoadHits,winningRoad); 


		BxByAxAy_AllRoads.erase(BxByAxAy_AllRoads.begin()+bestroadIndex);
		Atleast_aHitRoad_Created=true;
	
	      }
	    else//No tracks is good enough. Start to remove outliers until a best candidate is found
	      {
		if(verbosity>=1)
		  std::cout<<"Trk Is not good: P="<<chi2XProb<<" Trk size:"<<ALLCL1RoadHits.at(bestroadIndex).size()<<std::endl;

		bool goodtrkobtained=false;
		bool MinimumNumCl1HitReached=false;
		//Remove outlaiers to each candidate until a good trk is found or all trks are under thresholds
		int bestroadIndex2=bestroadIndex;

		//I'm going to remove one Hit, usually  MinimumNumCl1Hit=3
		if(((int)ALLCL1RoadHits.at(bestroadIndex2).size()<= (MinimumNumCl1Hit+1))) 
		  MinimumNumCl1HitReached=true;

		while((goodtrkobtained==false)&&(MinimumNumCl1HitReached==false))
		  {		    	    		    
		    if((int)ALLCL1RoadHits.at(bestroadIndex2).size()>= (MinimumNumCl1Hit+1)) //+1 because I'm going to remove one Hit
		      {
			//std::cout<<"Remove outliers"<<std::endl;
			RemoveOutliers(ALLCL1RoadHits,BxByAxAy_AllRoads);
			//std::cout<<"Trk now have "<<ALLCL1RoadHits.at(bestroadIndex2).size()<<" Hits"<<std::endl;
			bestroadIndex2=GiveMeBestCandidate(ALLCL1RoadHits, chi2XProb, chi2YProb,BxByAxAy_AllRoads);
			//std::cout<<"check after GiveBestCandidate: "<<chi2XProb<<" "<<ALLCL1RoadHits.at(bestroadIndex2).size()<<std::endl;
			if((chi2XProb>=chi2XYProb_Thr))
			  {
			    goodtrkobtained=true;
			    //Promoved as final and removed from *allcandidateRoads_cleaned*
			    SuccessFull_CL1RoadHits.push_back(ALLCL1RoadHits.at(bestroadIndex2));
			    SuccessFull_BxByAxAy_Roads.push_back(BxByAxAy_AllRoads.at(bestroadIndex2));
			    //  std::cout<<" @@@@@@ road successfully cleaned->New good CollSize:"<<SuccessFull_CL1RoadHits.size()<<std::endl;
			    std::vector<T2Hit> winningRoad=ALLCL1RoadHits.at(bestroadIndex2);
			    ALLCL1RoadHits.erase(ALLCL1RoadHits.begin()+bestroadIndex2);

			    RemoveHitsFromLoserCombinations(ALLCL1RoadHits,winningRoad);
			    
			    BxByAxAy_AllRoads.erase(BxByAxAy_AllRoads.begin()+bestroadIndex2);
			    Atleast_aHitRoad_Created=true;
			  }
			//	else
			//std::cout<<"Trk Is Still not good: P="<<chi2XProb<<" Trk size:"<<ALLCL1RoadHits.at(bestroadIndex2).size()<<std::endl;
		
		      }
		    else
		      {MinimumNumCl1HitReached=true;}	
		  }

		//Bug discovered in  23/1/2013. We can loose cl1hit road in case of nothing successfull. This is true
		// especially in misaligned cases.
		////////////////////////////////////////Begin of the bug-patch
		if(goodtrkobtained==false)//It means that you have stripped all hits but saved nothing in SuccessFull_CL1RoadHits
		  if(Atleast_aHitRoad_Created==false)//It means that this is not a garbage of some other good road (avoid also ghosts)
		    if(MinimumNumCl1HitReached){//It means that after this there is no more chance to have a good track.
		      SuccessFull_CL1RoadHits.push_back(ALLCL1RoadHits.at(bestroadIndex2));
		      SuccessFull_BxByAxAy_Roads.push_back(BxByAxAy_AllRoads.at(bestroadIndex2));		     
		      std::vector<T2Hit> winningRoad=ALLCL1RoadHits.at(bestroadIndex2);
		      ALLCL1RoadHits.erase(ALLCL1RoadHits.begin()+bestroadIndex2);		     
		      RemoveHitsFromLoserCombinations(ALLCL1RoadHits,winningRoad);		     
		      BxByAxAy_AllRoads.erase(BxByAxAy_AllRoads.begin()+bestroadIndex2);
		      Atleast_aHitRoad_Created=true;
		      allCandidateExtracted=true;
		    }
		////////////////////////////////////////End of the bug-patch

		if(MinimumNumCl1HitReached)
		  allCandidateExtracted=true;
		
	      }//if(Trk chi2 is bad)
	    
	    //Check how many Roads has enough CL1-hits to continue.
	    //unsigned int candidateleft=0;
	    unsigned int thesize=ALLCL1RoadHits.size();
	    if(verbosity>=1){
	      std::cout<<"Successfull CL1 roads after Cl1-selection-A:"<<ALLCL1RoadHits.size()<<std::endl;
	      if(ALLCL1RoadHits.size()>0)
		std::cout<<"1st road size:"<<ALLCL1RoadHits.at(0).size()<<std::endl;
	    }
	    std::vector<unsigned int> IdtoremoveAgain; IdtoremoveAgain.clear();
	    std::vector<std::vector<T2Hit> >ALLCL1RoadHits_copy; ALLCL1RoadHits_copy.clear();
	    for(unsigned int jj=0;jj<thesize;jj++)
	      {
		

		if((int)ALLCL1RoadHits.at(jj).size()>=MinimumNumCl1Hit)
		  {
		    ALLCL1RoadHits_copy.push_back(ALLCL1RoadHits.at(jj));
		   
		  }
	      }
	    
	    if(ALLCL1RoadHits_copy.size()==0)
	      allCandidateExtracted=true;

	    ALLCL1RoadHits.clear(); 
	    ALLCL1RoadHits=ALLCL1RoadHits_copy;
	    if(verbosity>=1)
	      std::cout<<"Successfull CL1 roads after Cl1-selection-A2:"<<ALLCL1RoadHits.size()<<std::endl;

	  }//while(allCandidateExtracted==true)
	
      }//if(mapZvsUniqueHitId.size() >= MinimumNumCl1Hit)
    
    
    

	
    unsigned int planesimbol=0; //used to understan where is the quarter.
    unsigned int startingIndex_QuartertoLook=0;
    
    //ReJoin the winning CL1 tracks with eventuals PadONLY Class0 Hits of thestarting pad Road.
    unsigned int numcl1Hit=0;
   
    if(verbosity>=1)
      std::cout<<"Successfull CL1 roads after Cl1-selection-B:"<<SuccessFull_CL1RoadHits.size()<<std::endl;
    
    for(unsigned int t=0;t<SuccessFull_CL1RoadHits.size();t++)
      {
	
	for(unsigned int uu=0;uu<PadonlyHit.size();uu++)
	  {
	    //Add at each cl1 road the padcluster (class 0 Hits)
	    // std::cout<<"Pad-only added"<<std::endl;
	    SuccessFull_CL1RoadHits.at(t).push_back(PadonlyHit.at(uu));	    	    
	  }

	//e)Now search for compatible strips of the pad +CL1Road 
	//  using SuccessFull_CL1RoadHits SuccessFull_BxByAxAy_Roads: one problem (strategy to be improved:)
	//  How to use strip information at this stage during prediction from lower plane to higher plane.
	//  For this moment strip are collected only using initial information, without updating the cl1 road fits,
	//  so with small prediction power. It is usefull only to lower match criteria.One method could be. Calculate the best Phi at 
	// the best matching R !!
	




	std::vector<unsigned int> planewhitAlreadySomething; planewhitAlreadySomething.clear();

	for(unsigned int j=0;j<SuccessFull_CL1RoadHits.at(t).size();j++){

	  if((SuccessFull_CL1RoadHits.at(t).at(j).GetHitClass()==1)||(SuccessFull_CL1RoadHits.at(t).at(j).GetHitNumPad()>0))
	    numcl1Hit++;

	  planeinfo=convGeo.GetT2Info(SuccessFull_CL1RoadHits.at(t).at(j).GetHitDetRawId());
	  planesimbol=planeinfo.symb;
	  planewhitAlreadySomething.push_back(planesimbol);
	}
	
	startingIndex_QuartertoLook=(planewhitAlreadySomething.at(0)/10)*10;
	//std::cout<<"1)Size road: "<<SuccessFull_CL1RoadHits.at(t).size()<<" passed to fitTrackspecialXY "<<std::endl;
	T1T2Track thisTrk=fitTrackspecialXY(SuccessFull_CL1RoadHits.at(t)); //Should be Improved it taking into account also strip
	
	//scan all planes in the quarter.
	
	for(unsigned int stripplane=startingIndex_QuartertoLook; stripplane<(startingIndex_QuartertoLook+10);stripplane++){	  	    
	  /*double thisR=0.;*/ double closerDR=100.; T2Hit CloserStripHit; double thisDR=0.;
	  
	    if(std::find(planewhitAlreadySomething.begin(), planewhitAlreadySomething.end(), stripplane)==planewhitAlreadySomething.end()){
	      //This is a plane where you can recover a trip. I look in the T2HitCollection for the strip with best R-Phi respect
	      //to the Trk Projection.
	     
	      uint32_t rawiddet=0;
	      unsigned int  symb=0;
	      
	      double projX=0.;double projY=0.;
	      for(T2HitCollection::const_iterator ithit = T2HitCollectionPtr->begin();ithit!=T2HitCollectionPtr->end();ithit++){
		rawiddet=(*ithit).GetHitDetRawId();
		planeinfo=convGeo.GetT2Info(rawiddet);//RawtoSymb(rawiddet);
		symb=planeinfo.symb;
		if((symb==stripplane)&&((*ithit).GetHitNumPad()==0))
		  {
		    //std::cout<<"Candidate StripR"<<<<" Strip-Plane:"<<stripplane<<std::endl;
		    projX=thisTrk.X0()+thisTrk.GetTx()*(*ithit).GetHitZ();
		    projY=thisTrk.Y0()+thisTrk.GetTy()*(*ithit).GetHitZ();
		    thisDR=fabs(sqrt(projX*projX+projY*projY)-(*ithit).GetHitR());
		    
		    if(thisDR<closerDR){
		      CloserStripHit=(*ithit);
		      closerDR=thisDR;
		    }
		    
		  }

		}	      
	    }
	    
	    if(closerDR<1.5){ //To review better
	      // std::cout<<"Strip-only added"<<std::endl;
	      SuccessFull_CL1RoadHits.at(t).push_back(CloserStripHit);
	    }
	    //else
	    //std::cout<<"Plane:"<<stripplane<<" Max close DR="<<closerDR<<std::endl;
	    
	      
	  }


	
	 //f) copy SuccessFull_CL1RoadHits in theT2Roads_withStrips if they are NOT CLONES candidate
	 if(numcl1Hit>=3)
	   {
	     //g) only who match criteria will survive.
	     // SuccessFull_CL1RoadHits.push_back(SuccessFull_CL1RoadHits.at(t));
	     ConcurrentmapIt = Map_ConcurrentTubeId.find(itroad->RoadID);


	     if(ConcurrentmapIt==Map_ConcurrentTubeId.end())
	       {
		 T2Road theRoadwithStrip;
		 theRoadwithStrip.thisPadRoad =itroad->thisPadRoad;
		 theRoadwithStrip.RoadID =  10000*t + itroad->RoadID; //So you can reassociate the initial pad-road.
		 theRoadwithStrip.thisRoad=(SuccessFull_CL1RoadHits.at(t));

		 if(verbosity>=1)
		   std::cout<<">>>>>>>>>>>>>>>> final Road with first hit R-Phi-Z: "<<SuccessFull_CL1RoadHits.at(t).at(0).GetHitR()<<" "<<SuccessFull_CL1RoadHits.at(t).at(0).GetHitPhi()<<" "<<SuccessFull_CL1RoadHits.at(t).at(0).GetHitZ()<<" saved"<<std::endl;

		 theT2Roads_withStrips->push_back(theRoadwithStrip);
		 

		 //--------------------------------------------------------------------//
		 //------------   Part for save only overlap free tracks --------------//
		 if((ResolveOverlapDoubleCount==true)&&(TrkInOverlapRegion(theRoadwithStrip)==false)){
		   theT2Roads_withStripsOverlapFree->push_back(theRoadwithStrip);
		   //std::cout<<"Sav:NotOV"<<std::endl;
		 } else if((ResolveOverlapDoubleCount==true)&&(TrkInOverlapRegion(theRoadwithStrip)==true)){
		   //Put it in a container
		   //std::cout<<"NotSav-OV"<<std::endl;
		   RoadInOverlapRegion.push_back(SuccessFull_CL1RoadHits.at(t));
		   RoadInOverlapRegion_index.push_back(theRoadwithStrip.RoadID);
		 }else if(ResolveOverlapDoubleCount==false){
		   theT2Roads_withStripsOverlapFree->push_back(theRoadwithStrip);

		 
		   // std::cout<<"Sav:Dontcare"<<std::endl;
		 }		   
		 //---------------------------------------------------------------------//
	       
	       }
	     else
	       Map_cloneRoadID_ToRoadPosIndex[itroad->RoadID].push_back(t);
	   }


      
      }//for(unsigned int t=0;t<SuccessFull_CL1RoadHits.size();t++)
   
    
    
    
  }//T2PadRoadCollection Loop End
 

//Clone Finder and writer. Road belonging to same beninning tube are in vector
//SuccessFull_CL1RoadHits with index in Road_posindxForCloneKiller
//Map_cloneRoadID_ToRoadPosIndex

 std::map<int, std::vector<unsigned int> >::iterator Map_cloneRoadID_ToRoadPosIndex_Iter;
 for(Map_cloneRoadID_ToRoadPosIndex_Iter = Map_cloneRoadID_ToRoadPosIndex.begin();Map_cloneRoadID_ToRoadPosIndex_Iter!=Map_cloneRoadID_ToRoadPosIndex.end();Map_cloneRoadID_ToRoadPosIndex_Iter++)
   {
     if(Map_cloneRoadID_ToRoadPosIndex_Iter->second.size()==1){
       T2Road theRoadwithStrip;
       //theRoadwithStrip.thisPadRoad = SuccessFull_CL1RoadHits.at(0); //itroad->thisPadRoad;
       theRoadwithStrip.RoadID =  10000*Map_cloneRoadID_ToRoadPosIndex_Iter->second.at(0) + Map_cloneRoadID_ToRoadPosIndex_Iter->first; 
       //So you can reassociate the initial pad-road.
       theRoadwithStrip.thisRoad=(SuccessFull_CL1RoadHits.at(Map_cloneRoadID_ToRoadPosIndex_Iter->second.at(0)));
       theT2Roads_withStrips->push_back(theRoadwithStrip);


       //--------------------------------------------------------------------//
       //------------   Part for save only overlap free tracks --------------//
       if((ResolveOverlapDoubleCount==true)&&(TrkInOverlapRegion(theRoadwithStrip)==false)){
	 theT2Roads_withStripsOverlapFree->push_back(theRoadwithStrip);
	 //std::cout<<"B-Sav:NotOv"<<std::endl;
       } else if((ResolveOverlapDoubleCount==true)&&(TrkInOverlapRegion(theRoadwithStrip)==true)){
	 //Put it in a container
	 // std::cout<<"B-NotSav:Ov"<<std::endl;
	 RoadInOverlapRegion.push_back((SuccessFull_CL1RoadHits.at(Map_cloneRoadID_ToRoadPosIndex_Iter->second.at(0))));
	 RoadInOverlapRegion_index.push_back(theRoadwithStrip.RoadID);
       }else if(ResolveOverlapDoubleCount==false){
	 theT2Roads_withStripsOverlapFree->push_back(theRoadwithStrip);

	 //Save here roadID-hit association
	 for(unsigned int as=0;as< theRoadwithStrip.thisRoad.size(); as++){
	   long int posincoll= (long int) theRoadwithStrip.thisRoad.at(as).Hit_PosInCollection;	   
	   theT2Hit_to_Track_Map->insert ( pair<long int,int>(posincoll,theRoadwithStrip.RoadID) );
	 }
	 // std::cout<<"B-dontcare"<<std::endl;
       }
       //--------------------------------------------------------------------//



     }//save directly the road

     
     if(Map_cloneRoadID_ToRoadPosIndex_Iter->second.size()==0){
       std::cout<<"CloneError-Message0"<<std::endl;
     }

   
     if(Map_cloneRoadID_ToRoadPosIndex_Iter->second.size()>1){ //Activate ConcurrentBranchResolver
      
       ParallelBranches.clear();
       for(unsigned int conc=0;conc<Map_cloneRoadID_ToRoadPosIndex_Iter->second.size();conc++)
	 {
	   ParallelBranches.push_back(SuccessFull_CL1RoadHits.at((Map_cloneRoadID_ToRoadPosIndex_Iter->second.at(conc))));
	   
	 }
       
       std::vector<std::vector<T2Hit> > ResolvedRoads=ConcurrentBranchResolver(ParallelBranches);
       for(unsigned int i=0;i<ResolvedRoads.size();i++)
	 {
	    T2Road theRoadwithStrip;
	    //I will use a negative id for this roads!	    
	    theRoadwithStrip.RoadID =  -10000*i + Map_cloneRoadID_ToRoadPosIndex_Iter->first; 
	    theRoadwithStrip.thisRoad=ResolvedRoads.at(i);
	    theT2Roads_withStrips->push_back(theRoadwithStrip);
	    
	    //--------------------------------------------------------------------//
	    //------------   Part for save only overlap free tracks --------------//
	    if((ResolveOverlapDoubleCount==true)&&(TrkInOverlapRegion(theRoadwithStrip)==false)){
	      theT2Roads_withStripsOverlapFree->push_back(theRoadwithStrip);
	      // std::cout<<"C-Savenotov"<<std::endl;
	    } else if((ResolveOverlapDoubleCount==true)&&(TrkInOverlapRegion(theRoadwithStrip)==true)){
	      //Put it in a container
	      //std::cout<<"C-dontsaveOV"<<std::endl;
	      RoadInOverlapRegion.push_back(ResolvedRoads.at(i));
	      RoadInOverlapRegion_index.push_back(theRoadwithStrip.RoadID);
	    }else if(ResolveOverlapDoubleCount==false){
	      theT2Roads_withStripsOverlapFree->push_back(theRoadwithStrip);
	      //std::cout<<"C-dontcare"<<std::endl;
	    }
	    //--------------------------------------------------------------------//
       


	 }
       
       
     }
    
   }

 //Last step: Resolve RoadInOverlapRegion, merge it and put in theT2Roads_withStripsOverlapFree 
 OverlapRoadMerged.clear();
 //std::cout<<" #ReadyToMerge:"<<RoadInOverlapRegion.size()<<" #Already saved:"<<theT2Roads_withStripsOverlapFree->size()<<std::endl;
 double Deta=0.; double dphi=0.; double DR=0.;
 if(ResolveOverlapDoubleCount==true){
   for(unsigned int ov1=0;ov1<RoadInOverlapRegion.size();ov1++){
     
     bool ov1merged=false;
     int bestov1_tomerge=0; int bestov2_tomerge=0;
     double bestDeta=100.; double bestdphi=100.; double bestDR=100.;
 
      for(unsigned int ov2=ov1;ov2<RoadInOverlapRegion.size();ov2++){
	
	
	if(ov1!=ov2){
	  if(SameTrksInOverlap(RoadInOverlapRegion.at(ov1),RoadInOverlapRegion.at(ov2),Deta,dphi,DR)){
	    
	    if((fabs(Deta)<bestDeta)&&(fabs(dphi)<bestdphi)&&(fabs(DR)<bestDR)){
	      ov1merged=true;
	      bestDeta=fabs(Deta); bestdphi=fabs(dphi); bestDR=fabs(DR);
	      bestov1_tomerge=ov1;  
	      bestov2_tomerge=ov2;
	    }
	    
	    
	  }
	  
	}
      
      }
      
      if(ov1merged==false){//No counterpart found for the merging. Saved as an non-ov track.
	//OverlapRoadMerged.push_back(RoadInOverlapRegion.at(ov1));
	 T2Road theRoadwithStrip;
	    //I will use a negative id for this roads!	    
	 theRoadwithStrip.RoadID = RoadInOverlapRegion_index.at(ov1);
	 theRoadwithStrip.thisRoad=RoadInOverlapRegion.at(ov1);
	 theT2Roads_withStripsOverlapFree->push_back(theRoadwithStrip);	
	 RoadInOverlapRegion.erase(RoadInOverlapRegion.begin()+ov1);
	 //erase RoadInOverlapRegion.at(ov1)
      }
      else{ 
	//merge the roads choosing the closer distances.
	
	T2Road theRoadwithStrip;
		    
	theRoadwithStrip.RoadID = RoadInOverlapRegion_index.at(bestov1_tomerge);
	theRoadwithStrip.thisRoad=RoadInOverlapRegion.at(bestov1_tomerge);
	for(unsigned int j=0;j<RoadInOverlapRegion.at(bestov2_tomerge).size();j++){
	  theRoadwithStrip.thisRoad.push_back(RoadInOverlapRegion.at(bestov2_tomerge).at(j)); 
	}
	//std::cout<<"Saved oV. road with size:"<<theRoadwithStrip.thisRoad.size()<<std::endl;

	theT2Roads_withStripsOverlapFree->push_back(theRoadwithStrip);
	RoadInOverlapRegion.erase(RoadInOverlapRegion.begin()+bestov2_tomerge); 
	//Since the index bestov1_tomerge will be not reconsidered anymore	
	
      }
      
   
   }
 } 
  

 /* }*/


/*
 for(unsigned int t=0;t<theT2Roads_withStrips->size();t++)
   {
     //std::cout<<"-----"<<theT2Roads_withStrips->at(t).RoadID<<"-----"<<std::endl;
     for(unsigned int u=0;u<theT2Roads_withStrips->at(t).thisRoad.size();u++)
       std::cout<<theT2Roads_withStrips->at(t).thisRoad.at(u).GetHitZ()<<" "<<theT2Roads_withStrips->at(t).thisRoad.at(u).GetHitR()<<" ||";

     std::cout<<std::endl;
   }
*/

 event.put(theT2Roads_withStripsOverlapFree, T2RoadCollProdName);
 event.put(theT2Hit_to_Track_Map,"HitToRoadAssociator");

} 

// END produce




bool T2RoadPadFinder::SameTrksInOverlap(std::vector<T2Hit> &RoadInOverlapRegion1,std::vector<T2Hit> &RoadInOverlapRegion2,double &Deta, double &Dphi, double &DR){
  bool arethesame=false;
  //Strategy: Compute the minimum Phi distance.
  //Compute the average Eta.
  //Compute the minimum R.
  //If ALL This are smaller than the window the track are the same.
  double phi1=0.;
  double phi2=0.;

  double avgcount1=0.;
  double avgcount2=0.;
  double etaavg1=0.;
  double etaavg2=0.;
  double Rdist=100.; double Phidist=100.;
 
  uint32_t rawiddet=0;
  unsigned int  symb1=0;unsigned int  symb2=0;	      
  if(RoadInOverlapRegion1.at(0).GetHitClass()==1)    
    phi1=RoadInOverlapRegion1.at(0).GetHitPhi();

  rawiddet=RoadInOverlapRegion1.at(0).GetHitDetRawId();
  planeinfo=convGeo.GetT2Info(rawiddet);//RawtoSymb(rawiddet);
  symb1=planeinfo.symb;
 
  
  if(RoadInOverlapRegion2.at(0).GetHitClass()==1)    
    phi2=RoadInOverlapRegion2.at(0).GetHitPhi();
  
  rawiddet=RoadInOverlapRegion2.at(0).GetHitDetRawId();
  planeinfo=convGeo.GetT2Info(rawiddet);//RawtoSymb(rawiddet);
  symb2=planeinfo.symb;
  

  if(((symb1/20)==(symb2/20))&&((symb1/10)!=(symb2/10))&&(fabs(phi1-phi2)<20.))    
    for(unsigned int r1=0;r1<RoadInOverlapRegion1.size();r1++){
      
      if(RoadInOverlapRegion1.at(r1).GetHitClass()==1){	
	etaavg1=etaavg1+RoadInOverlapRegion1.at(r1).GetHitR()/RoadInOverlapRegion1.at(r1).GetHitZ();
	avgcount1=avgcount1+1.0;
      }

      for(unsigned int r2=0;r2<RoadInOverlapRegion2.size();r2++){
	if(r1==0)
	  {
	    //compute eta2Avg;
	    if(RoadInOverlapRegion2.at(r2).GetHitClass()==1){
		etaavg2=etaavg2+RoadInOverlapRegion2.at(r2).GetHitR()/RoadInOverlapRegion2.at(r2).GetHitZ();
		avgcount2=avgcount2+1.0;
	      }
	  }

	if(fabs(RoadInOverlapRegion1.at(r1).GetHitR()-RoadInOverlapRegion2.at(r2).GetHitR())<Rdist)
	  Rdist=fabs(RoadInOverlapRegion1.at(r1).GetHitR()-RoadInOverlapRegion2.at(r2).GetHitR());
	
	if(fabs(RoadInOverlapRegion1.at(r1).GetHitPhi()-RoadInOverlapRegion2.at(r2).GetHitPhi())<Phidist)
	  Phidist=fabs(RoadInOverlapRegion1.at(r1).GetHitPhi()-RoadInOverlapRegion2.at(r2).GetHitPhi());  

      }

    }
    
  etaavg1=etaavg1/avgcount1;
  etaavg2=etaavg2/avgcount2;
  
  double etadist=fabs(etaavg1-etaavg2);
  Deta=etadist;
  Dphi=Phidist;
  DR=Rdist;

  //std::cout<<"SameTrksInOverlap out:"<<arethesame<<"   DTheta-DR-Dphi= "<<Deta<<"-"<<DR<<"-"<<Dphi<<std::endl;
  //std::cout<<"Window-cut:"<<OverlapDoubleCountDTheta<<"-"<<OverlapDoubleCountDR<<"-"<<OverlapDoubleCountDPhi<<std::endl;


  if(Deta<OverlapDoubleCountDTheta){
    //std::cout<<"Deta OK"<<std::endl;
    if(DR<OverlapDoubleCountDR){
      //std::cout<<"DR OK"<<std::endl;
      if(Dphi<OverlapDoubleCountDPhi){
	//std::cout<<"DPhi OK"<<std::endl;
	arethesame=true;
      }
    }
  }
  
  //if(arethesame)
  //std::cout<<"SameTrksInOverlap out: TO MERGE!"<<std::endl;
  return arethesame;
}


bool T2RoadPadFinder::TrkInOverlapRegion(T2Road &aroad){

  //Strategy: true if one of the 2 cl1 hit phi checked are in the overlapping region. 

  bool trkinoverlap=false;
  
  unsigned int hitsize=aroad.thisRoad.size();
  unsigned int howmanycount=0;
  double phi1=0.; double phi2=0.;
  
  unsigned int h=0;
  while((h<hitsize)&&(howmanycount<=2)){
    if(aroad.thisRoad.at(h).GetHitClass()==1){
      
      if(howmanycount==0)
	phi1=aroad.thisRoad.at(h).GetHitPhi();
      else
	phi2=aroad.thisRoad.at(h).GetHitPhi();


      howmanycount++;
    }

    h++;
  } 

  if((phi1<=96.)&&(phi1>=84.))
    trkinoverlap=true;
  
  if((phi1<=276.)&&(phi1>=264.))
    trkinoverlap=true;
  
  if((phi2<=96.)&&(phi2>=84.))
    trkinoverlap=true;
  
  if((phi2<=276.)&&(phi2>=264.))
    trkinoverlap=true;
  
  //  std::cout<<"Trk in ov?"<<trkinoverlap<<"  :"<<phi1<<"-"<<phi2<<std::endl;

  return trkinoverlap;
}






std::vector<std::vector<T2Hit> >  T2RoadPadFinder::ConcurrentBranchResolver(std::vector<std::vector<T2Hit> > &ParallelBranches)
{
  // std::vector<std::vector<T2Hit> > distinctRoads; distinctRoads.clear();
  //Strategy:
  
  //At the moment: takes NON shared Hits from all the road.
  //Start from larger unique-Hit Road -> try to associate to them the shared one.Remove outliers if 
  std::map<pair<int,int>, std::vector<unsigned int> > Map_HitID_To_PosindexofRoads;
  // std::map<unsigned int,unsigned int> Map_RoadPosID_To_VectorUniqueID;
  std::map<pair<int,int>, std::vector<unsigned int> >::const_iterator MapIter;

  std::map<pair<int,int>, T2Hit> mapForTakebackHit;
  mapForTakebackHit.clear();

  for(unsigned int p=0;p<ParallelBranches.size();p++)
    {
      for(unsigned int p2=0;p2<ParallelBranches.size();p2++)
	{
	  Map_HitID_To_PosindexofRoads[ParallelBranches.at(p).at(p2).HitUniqueId].push_back(p);
	  mapForTakebackHit[ParallelBranches.at(p).at(p2).HitUniqueId]=ParallelBranches.at(p).at(p2);
	}
      
    }
  
  std::vector<std::vector<T2Hit> > uniqueHitsRoad; uniqueHitsRoad.clear();
  std::vector<T2Hit> aRoad;
  std::vector<unsigned int> cl1HitCounter;cl1HitCounter.clear();
  for(unsigned int p=0;p<ParallelBranches.size();p++)
    {
      aRoad.clear();
      cl1HitCounter.push_back(0);
      for(unsigned int p2=0;p2<ParallelBranches.size();p2++)
	{
	  MapIter=Map_HitID_To_PosindexofRoads.find(ParallelBranches.at(p).at(p2).HitUniqueId);
	  if(MapIter!=Map_HitID_To_PosindexofRoads.end())
	    if(MapIter->second.size()==1)
	      {
		aRoad.push_back(ParallelBranches.at(p).at(p2));
		if(ParallelBranches.at(p).at(p2).GetHitClass()==1)
		  cl1HitCounter.at(p)=cl1HitCounter.at(p)+1;
	      }
	    
	}

      if(aRoad.size()>=3)
	uniqueHitsRoad.push_back(aRoad);
    }
  
  std::vector<std::vector <T2Hit> > ArtifactRoads;
  std::vector<std::vector <double> > ArtifactFitParam;
  
  for(MapIter=Map_HitID_To_PosindexofRoads.begin();MapIter!=Map_HitID_To_PosindexofRoads.end();MapIter++)
    {
      if(MapIter->second.size()>1)
	{ 
	  //decide Who assign the double Hit. I will use association according to Final LessChi2Prob 
	  //Only Cl1Hit will be reassociated.
	  if(mapForTakebackHit[MapIter->first].GetHitClass()==1)
	    {
	      std::vector<double> allFinalChi2Prob; allFinalChi2Prob.clear();
	      double maxprob=0.; unsigned int maxprobind=0;
	      allFinalChi2Prob.clear();
	      for(unsigned int i=0;i<uniqueHitsRoad.size();i++)
		{
		  double finalchiX=0;double finalchiY=0;
		  ArtifactRoads.clear();
		  ArtifactFitParam.clear();
		  ArtifactRoads.push_back(uniqueHitsRoad.at(i));
		  GiveMeBestCandidate(ArtifactRoads,finalchiX,finalchiY,ArtifactFitParam);
		  allFinalChi2Prob.push_back(finalchiX);
		  if(finalchiX>maxprob)
		    {
		      maxprob=finalchiX;
		      maxprobind=i;
		    }
		}

	      //Assign the Hit to the max-prob chi2x
	      
	      if(maxprob>=chi2XYProb_Thr)
		uniqueHitsRoad.at(maxprobind).push_back(mapForTakebackHit[MapIter->first]);
	    }
	  
	}

    }
  
  return uniqueHitsRoad;
}


void T2RoadPadFinder::RemoveHitsFromLoserCombinations(std::vector<std::vector<T2Hit> > &allcandidateRoadsCL1,std::vector<T2Hit> &winnigRoad)
{

  std::vector<T2Hit> tobecopied; tobecopied.clear();

  std::map<pair<long int, long int>, T2Hit> surviveHits; surviveHits.clear();
  
  for(unsigned int i=0;i<allcandidateRoadsCL1.size();i++)   
    {
      tobecopied.clear();
      for(unsigned int j=0;j<allcandidateRoadsCL1.at(i).size();j++)
	{
	  bool toberased=false;
	  for(unsigned int w=0;w<winnigRoad.size();w++)
	    if(allcandidateRoadsCL1.at(i).at(j).HitUniqueId==winnigRoad.at(w).HitUniqueId)
	      {
		//allcandidateRoadsCL1.at(i).erase(allcandidateRoadsCL1.at(i).begin()+j);
		toberased=true;
	      }
	  
	  if(toberased==false)
	    {
	      
	      tobecopied.push_back(allcandidateRoadsCL1.at(i).at(j));
	      surviveHits[allcandidateRoadsCL1.at(i).at(j).HitUniqueId]=allcandidateRoadsCL1.at(i).at(j);
	      //int padHrdw=realHardwCluId
	      //int padsymbidusedbefore=Map_RealClusterID_ToSymbUsed[realHardwCluId];
	      //double planeZ=Map_ClustuniqueId_ZPos[padsymbidusedbefore];

	    }
	}
      allcandidateRoadsCL1.at(i)=tobecopied;
    }

  //Added on 23/12/2010
  //Here I want to recompute the combination and save in allcandidateRoadsCL1 basicly this is because I want to restart looking at 
  //the bigger size combination.

  
  
  std::map<double, std::vector<T2Hit> > MapZ_to_surviveHits; MapZ_to_surviveHits.clear();
  std::map<double, std::vector<unsigned long int> > MapZ_to_surviveHits_supportComb; MapZ_to_surviveHits_supportComb.clear();

  std::map<pair<long int,long int>, T2Hit>::iterator surviveHits_Iter;
    
  unsigned long int localHitIndex=0; 


  std::vector<T2Hit> survivedHitsVect;survivedHitsVect.clear();

  for(surviveHits_Iter=surviveHits.begin();surviveHits_Iter!=surviveHits.end();surviveHits_Iter++)
    {
      T2Hit survHit= surviveHits_Iter->second;

      double theZ=survHit.GetHitZ();
      survivedHitsVect.push_back(survHit);

      MapZ_to_surviveHits[theZ].push_back(survHit); 
      
      localHitIndex=survivedHitsVect.size()-1;
      MapZ_to_surviveHits_supportComb[theZ].push_back(localHitIndex); 
      
      
    }
  
   std::vector<std::vector <unsigned long int> > allSymbcombination= RoadPadCL1HitsCombination(MapZ_to_surviveHits_supportComb);

  allcandidateRoadsCL1.clear();
  for(unsigned int u1=0;u1<allSymbcombination.size();u1++){
    std::vector <T2Hit> onecombination;
    onecombination.clear();
    for(unsigned int u2=0;u2<allSymbcombination.at(u1).size();u2++){
      unsigned long int index=(unsigned int) allSymbcombination.at(u1).at(u2);
      onecombination.push_back(survivedHitsVect.at(index));
    }
    allcandidateRoadsCL1.push_back(onecombination);
  }
  
  


}

// ----------------------------------------------------------------------
//  SUPPORTING FUNCTIONS FOR THE STRIP-MATCHING PART of the ALGORITHM
// ----------------------------------------------------------------------




int T2RoadPadFinder::GiveMeBestCandidate(std::vector<std::vector<T2Hit> > & allcandidateRoads, double &chi2XProb, double &chi2YProb,std::vector<std::vector<double> > &BxByAxAy_AllRoads)
{
  //Find best candidate index       
  BxByAxAy_AllRoads.clear();
  int bestcandidatexy=-1;
  
  //double chir,chiphi,chirprob,chiphiprob;
  //double bestchir,bestchiphi,bestchirprob,bestchiphiprob;

 double chi=0;double chiprob=0;
  double bestchiprob=0;



  std::vector<double> oneBxByAxAy;

  std::vector<T2Hit> thiscandidate;
  
  for(unsigned int r=0;r<allcandidateRoads.size();r++)
    {
   
      T1T2Track thetrkstoredXY(2);//T1T2Track thetrack(2);
      if(verbosity>=3)
	std::cout<<"GiveMeBestCandidate: This candidate size: ( "<< allcandidateRoads.at(r).size()<<" ) "<<std::endl;

      allcandidateRoads.at(r).size();
      
      	  
      //The strategy is to take the best between chixprob and chiyprob since I have seen that the two fit can 
      //have very different value of prob (0.00001 vs 0.5 for example) and I want a conservative track finding.

      //std::cout<<"2)Size road: "<<allcandidateRoads.at(r).size()<<"passed to fitTrackspecialXY "<<std::endl;
      
      //thetrkstoredXY=fitTrackspecialXY(allcandidateRoads.at(r));

      thetrkstoredXY=TrackerFitter(allcandidateRoads.at(r));//CHANGED ON 26/1/2011.

      oneBxByAxAy.clear();
      oneBxByAxAy.push_back(thetrkstoredXY.X0());oneBxByAxAy.push_back(thetrkstoredXY.Y0());oneBxByAxAy.push_back(thetrkstoredXY.GetTx());oneBxByAxAy.push_back(thetrkstoredXY.GetTy());
      BxByAxAy_AllRoads.push_back(oneBxByAxAy);

      thetrkstoredXY.ChiSquaredX();
      thetrkstoredXY.ChiSquaredY();
      (TMath::Prob(thetrkstoredXY.ChiSquaredX(),(thetrkstoredXY.GetHitEntries()-2)));
      (TMath::Prob(thetrkstoredXY.ChiSquaredY(),(thetrkstoredXY.GetHitEntries()-2)));
      chi=thetrkstoredXY.ChiSquared();
      chiprob=(TMath::Prob(chi,(2*thetrkstoredXY.GetHitEntries()-4)));
      
      
      /*//Commented  ON 26/1/2011.
      if(chixprob>chiyprob)
	chiprob=chixprob;
      else
	chiprob=chiyprob;
      */

      if(r==0)
	{
	  bestchiprob=chiprob;
	  bestcandidatexy=0;

	}
      else
	{
	  
	  //if(((chix<bestchix)&&(chiy<bestchiy))||((chixprob<bestchixprob)&&(chiyprob<bestchiyprob)))	
	  //if((chiprob<bestchiprob))
	  if((chiprob>bestchiprob))
	    {
	      bestchiprob=chiprob;

	      bestcandidatexy=r;


	    }
	  
	}      
     
      //  std::cout<<r<<" "<<chix<<" "<<chiy<<" "<<chixprob<<" "<<chiyprob<<std::endl;
    }
  
  chi2XProb=bestchiprob;//bestchixprob;
  chi2YProb=bestchiprob;

  //std::cout<<"GiveMeBestCandidate: Index:"<<bestcandidatexy<<"P: "<<chiprob<<std::endl;//<<" Eta:"<<thetrkstoredXY.Eta()<<" Phi:"<<thetrkstoredXY.Phi()*180./3.14159<<std::endl;
  

  return bestcandidatexy;
}




/*std::vector<std::vector<T2Hit> >*/void  T2RoadPadFinder::RemoveOutliers(std::vector<std::vector<T2Hit> > &analyzedRoadS,std::vector<std::vector<double> > &BxByAxAy_AllRoads)
{
  //  std::vector<std::vector<T2Hit> >   cleanedroadS; cleanedroadS.clean();
  
  for(unsigned int u=0;u<analyzedRoadS.size();u++)
    {
      int worstpoint=findWorstPoint(analyzedRoadS.at(u),BxByAxAy_AllRoads.at(u));
      if(worstpoint>=0)
	{
	  // erase the 6th element
	  // myvector.erase (myvector.begin()+5);
	  analyzedRoadS.at(u).erase(analyzedRoadS.at(u).begin()+worstpoint);
	  //BxByAxAy_AllRoads.at(u).erase( BxByAxAy_AllRoads.at(u).begin()+worstpoint);

	}
    }

}


 
 

std::vector<std::vector <unsigned long int> >  T2RoadPadFinder::RoadPadCL1HitsCombination(std::map<double, std::vector<unsigned long int> > &ZvsUniqueHitId_inroad)
{

  // std::vector<std::vector <long> > cl1Hitcombiation;
  //cl1Hitcombiation.clear();

  //This is the same strategy as in trackproducer2. Just changing T2Hits with its id (unsigned long int)
  // std::vector<unsigned long int> sameZhits;  
  std::vector<std::vector<unsigned long int> > hitmatrix;  hitmatrix.clear();
  
  
  std::vector<unsigned long int> candidatepiece; candidatepiece.clear();
  
  std::map<double, std::vector<unsigned long int> >::const_iterator ZvsUniqueHitId_inroad_iter;
 
  //Estimation of the total combination:
  
  

  std::vector<unsigned int> numcol_perrow;numcol_perrow.clear();
  std::vector<unsigned int> rowindex; rowindex.clear();
  unsigned int rowcounter=0;

  // Hit Matrix construction with in same row->same Z hit
  for (ZvsUniqueHitId_inroad_iter=ZvsUniqueHitId_inroad.begin();ZvsUniqueHitId_inroad_iter!=ZvsUniqueHitId_inroad.end();ZvsUniqueHitId_inroad_iter++)  
    {  
      numcol_perrow.push_back((*ZvsUniqueHitId_inroad_iter).second.size());
      rowindex.push_back(rowcounter);
      rowcounter++;
      //hitmatrix.push_back((*ZvsUniqueHitId_inroad_iter).second);           
    }
  
  
  bool reasonablenumberofcomb=false;
  
  while((reasonablenumberofcomb==false)&&(numcol_perrow.size()>0)){
   
    unsigned long int totalcomb=1;
    int maxrowentry=-1;
    unsigned int indexmaxrow=0;
    
    for(unsigned int i=0;i<numcol_perrow.size();i++){
      totalcomb*=numcol_perrow.at(i);
      if((int)numcol_perrow.at(i)>maxrowentry)
	{
	  maxrowentry=numcol_perrow.at(i);
	  indexmaxrow=i;
	}
    }

    //std::cout<<"HereB Total comb expected before suppression:"<<totalcomb<<std::endl;

    if(totalcomb<20000)
      reasonablenumberofcomb=true;
    else
      {
	std::vector<unsigned int> numcol_perrowLOWCOMB;
	std::vector<unsigned int> rowindexLOWCOMB;
      
	for(unsigned int i=0;i<numcol_perrow.size();i++){
	  if(i!=indexmaxrow){
	    rowindexLOWCOMB.push_back(rowindex.at(i));
	    numcol_perrowLOWCOMB.push_back(numcol_perrow.at(i));
	  }
	}
	numcol_perrow.clear();
	numcol_perrow=numcol_perrowLOWCOMB;
	rowindex.clear();
	rowindex=rowindexLOWCOMB;
      }

  }

  unsigned int rowcounterLOWCOMB=0;
  // Hit Matrix construction with in same row->same Z hit (below a reasonable number of combinations)
  for (ZvsUniqueHitId_inroad_iter=ZvsUniqueHitId_inroad.begin();ZvsUniqueHitId_inroad_iter!=ZvsUniqueHitId_inroad.end();ZvsUniqueHitId_inroad_iter++)  
    {  

      if(std::find(rowindex.begin(),rowindex.end(), rowcounterLOWCOMB)!= rowindex.end())      
	hitmatrix.push_back((*ZvsUniqueHitId_inroad_iter).second);   
      
      rowcounterLOWCOMB++;        
    }
  
  
   std::vector<std::vector<unsigned long int> > allcandidate;  allcandidate.clear();
                //allcandidate is a std::vector<std::vector<T2Hit> > ; vector of track-hits_ID
  for(unsigned int k=0;k<hitmatrix.size();k++)
    {
      //std::vector<unsigned long int> matrixraw=hitmatrix.at(k);
      
      if(allcandidate.size()==0){ // fill allcandidate with the first row of hitmatrix
	for(unsigned int l=0;l<hitmatrix.at(k).size();l++)
	  {
	    candidatepiece.push_back(hitmatrix.at(k).at(l));
	    allcandidate.push_back(candidatepiece);
	    
	    FreeVectorMemory(candidatepiece);candidatepiece.clear();
	  }
      }
      else
	{
	  unsigned int currentsz=allcandidate.size();
	  unsigned int m=0;
	  std::vector<std::vector<unsigned long int> > newcandidatestoput;
	  
	  currentsz=allcandidate.size();
	  
	  while(m<currentsz)
	    {
	      //join to all the candidates each element of the full next-Z row 
	      //so if the input was 2 candidate and the next Z have 3 elements at the end you have 2*3 candidates, each one differs for a next Z value.
	      candidatepiece=allcandidate.at(m);			 		       			 		     		     
	      for(unsigned int l=0;l<hitmatrix.at(k).size();l++)
		{
		  std::vector<unsigned long int> newcandidatepiece=candidatepiece;
		  newcandidatepiece.push_back(hitmatrix.at(k).at(l));
		  //allcandidate.push_back(newcandidatepiece);
		  newcandidatestoput.push_back(newcandidatepiece);
		  newcandidatepiece.clear();
		  
		}
	      FreeVectorMemory(candidatepiece);
	      candidatepiece.clear();
	      m++;
	    }
	  
	  
	  allcandidate=newcandidatestoput;
	  newcandidatestoput.clear();		     
	  
	}	               
      
    }

  

  return allcandidate;
}


/*




std::vector<std::vector <long> >  T2RoadPadFinder::RoadPadCL1HitsCombination(std::map<double, std::vector<long> > &ZvsUniqueHitId_inroad)
{

  // std::vector<std::vector <long> > cl1Hitcombiation;
  //cl1Hitcombiation.clear();

  //This is the same strategy as in trackproducer2. Just changing T2Hits with its id (long)
  // std::vector<long> sameZhits;  
  std::vector<std::vector<long> > hitmatrix;  hitmatrix.clear();
  std::vector<std::vector<long> > allcandidate;  allcandidate.clear();
  long lasthit;
  
  std::vector<long> candidatepiece; candidatepiece.clear();
  
  std::map<double, std::vector<long> >::const_iterator ZvsUniqueHitId_inroad_iter;

  // Hit Matrix construction with in same row->same Z hit
  for (ZvsUniqueHitId_inroad_iter=ZvsUniqueHitId_inroad.begin();ZvsUniqueHitId_inroad_iter!=ZvsUniqueHitId_inroad.end();ZvsUniqueHitId_inroad_iter++)  
    {      
      hitmatrix.push_back((*ZvsUniqueHitId_inroad_iter).second);      
    }
  
  //if(verbosity>=1)
    //for(unsigned int u=0;u<hitmatrix.size();u++)
      //{
	//std::cout<<" "<<std::endl;
//	std::vector<long> uxmatrixraw=hitmatrix.at(u);
	//for (unsigned int ii=0;ii<uxmatrixraw.size();ii++) 
	  //{
	    //std::cout<<"uxmatrixraw.at("<<ii<<") ="<<uxmatrixraw.at(ii)<<"  "<<std::endl;	 
	 // }
      //}
  

  allcandidate.clear();                 //allcandidate is a std::vector<std::vector<T2Hit> > ; vector of track-hits_ID
  for(unsigned int k=0;k<hitmatrix.size();k++)
    {
      //std::vector<long> matrixraw=hitmatrix.at(k);
      
      if(allcandidate.size()==0){ // fill allcandidate with the first row of hitmatrix
	for(unsigned int l=0;l<hitmatrix.at(k).size();l++)
	  {
	    candidatepiece.push_back(hitmatrix.at(k).at(l));
	    allcandidate.push_back(candidatepiece);
	    candidatepiece.clear();
	  }
      }
      else
	{
	  unsigned int currentsz=allcandidate.size();
	  unsigned int m=0;
	  std::vector<std::vector<long> > newcandidatestoput;
	  
	  currentsz=allcandidate.size();
	  
	  while(m<currentsz)
	    {
	      //join to all the candidates each element of the full next-Z row 
	      //so if the input was 2 candidate and the next Z have 3 elements at the end you have 2*3 candidates, each one differs for a next Z value.
	      candidatepiece=allcandidate.at(m);			 		       			 		     		     
	      for(unsigned int l=0;l<hitmatrix.at(k).size();l++)
		{
		  std::vector<long> newcandidatepiece=candidatepiece;
		  newcandidatepiece.push_back(hitmatrix.at(k).at(l));
		  //allcandidate.push_back(newcandidatepiece);
		  newcandidatestoput.push_back(newcandidatepiece);
		  newcandidatepiece.clear();
		  
		}
	      
	      candidatepiece.clear();
	      m++;
	    }
	  
	  
	  allcandidate=newcandidatestoput;
	  newcandidatestoput.clear();		     
	  
	}	               
      
    }

  
//if(verbosity>=1)
//  {
//    std::cout<<"Found "<<allcandidate.size()<<" combinations: "<<std::endl;
//    std::cout<<"All track candidates ( "<<allcandidate.size()<<"trk candidates ) in the road founds"<<std::endl;
//    for(unsigned int jj=0;jj<allcandidate.size();jj++)
//{
//  std::vector<long> onecandidate=allcandidate.at(jj);
//  for(unsigned int mm=0;mm<onecandidate.size();mm++)
//    { 
//      if(mm==0)
//	std::cout<<"Hit-Ids :  ";
//      
//      std::cout<<onecandidate.at(mm)<<"                          ";
//    }
//  std::cout<<" "<<std::endl;
//}
//    std::cout<<" |-------| "<<std::endl;
//  }    
  

  return allcandidate;
}


*/







int T2RoadPadFinder::findWorstPoint(std::vector<T2Hit> &trk,std::vector<double>  &BxByAxAy)
{

  std::vector<T2Hit> hitvec;  
  std::vector<float> r;std::vector<float> x;std::vector<float> y;
  std::vector<float> z;
  std::vector<float> er;
  std::vector<float> ez;
  //std::cout<<"A"<<std::endl;
  unsigned int sizeHitv=trk.size();
  bool biggestContributionChiRSet = false;
  double biggestContributionChiR = 0;
  unsigned int worstHitPointIndexR = 0;
  double chi2r = 0.0;  double chi2=0.;
  //Push the staff into vectors
  
  for (unsigned int jj = 0;jj<sizeHitv;jj++){
    T2Hit hit = trk.at(jj);


     TMatrixD OneVy(2,2); 
     OneVy.Zero();
     double phirad=trk.at(jj).GetHitPhi()*3.14159/180.0;
     double r= trk.at(jj).GetHitR(); 
     double sigmaR=trk.at(jj).GetHitDR();
     
      if(trk.at(jj).GetHitNumStrip()>4)
	sigmaR=(trk.at(jj).GetHitDR()+(trk.at(jj).GetHitNumStrip()-4)*0.4);
     
     double sigmaPhi=trk.at(jj).GetHitDPhi()*(2/sqrt(12.0))*3.14159/180.0;


     double x= trk.at(jj).GetHitX();  double y= trk.at(jj).GetHitY();       
     double z=trk.at(jj).GetHitZ();      
     OneVy(0,0)=sigmaR*sigmaR*cos(phirad)*cos(phirad)+r*r*sin(phirad)*sin(phirad)*sigmaPhi*sigmaPhi;
     OneVy(0,1)=cos(phirad)*sin(phirad)*(sigmaR*sigmaR-r*r*sigmaPhi*sigmaPhi);  
     OneVy(1,0)=OneVy(0,1);
     OneVy(1,1)=sigmaR*sigmaR*sin(phirad)*sin(phirad)+ r*r*cos(phirad)*cos(phirad)*sigmaPhi*sigmaPhi;
   
     //Invert Vy matrix
     TMatrixD OneVym1(2,2); 
     OneVym1=OneVy;
     Double_t deti;	
     OneVym1.Invert(&deti);

     
     TVectorD Residui(2);
     Residui(0)=(BxByAxAy.at(2)*z + BxByAxAy.at(0) - x);
     Residui(1)=(BxByAxAy.at(3)*z + BxByAxAy.at(1)  - y);
     
     TVectorD VdotRes(2);
     VdotRes=(OneVym1*Residui);
     chi2=Residui*VdotRes;
     
     chi2r += chi2;
    
    if(biggestContributionChiRSet == false || biggestContributionChiR < chi2){     
      biggestContributionChiRSet = true;
      biggestContributionChiR = chi2;
      worstHitPointIndexR = jj;
    }
     
  }
  // std::cout<<"B"<<std::endl;
  
 
  
  // if (TMath::Prob(chi2r,2*(sizeHitv-2)) >= chi2XYProb_Thr || biggestContributionChiRSet == false) {
    //Everything is OK. No worst point found.
  // return -1;
  // }
  
  /*
  if(verbosity>=1)
    std::cout<<"Worst index in trk: "<<(int) worstHitPointIndexR<<std::endl;
  */
   

  return (int) worstHitPointIndexR;
  
}




T1T2Track T2RoadPadFinder::TrackerFitter(std::vector<T2Hit> &hitvec2)
{

  TMatrixD par_covariance_matrix(4,4);
  
  
  int RoadID=0;
  T2Hit worstHit;



  std::vector<T2Hit> hitvec; hitvec.clear();
  bool Usestrip=false; 
 
  double sigmaR=(0.12+0.05);
  double sigmaPhi=0.015;
  unsigned int numphimin20=0;
  unsigned int numphimag340=0;

  double phirad=0.;
  double r=0.;


  bool inserted=false;

  unsigned int numHit_t2=0;
  unsigned int numCl1HitHit_t2=0;
  unsigned int numStripOnlyHit_t2=0;
  unsigned int numPadOnlyHit_t2=0;
 
  if(verbosity){
    std::cout<<" MyLinearfitCorrDEV start with "<<hitvec2.size()<<std::endl;
    
  }

  ////////////////////////////////////////////////////////////////////////////////
  ///           SAVE IN HITVEC ONLY WANTED HITS
  ///////////////////////////////////////////////////////////////////////////////

  for(unsigned int m=0;m<hitvec2.size();m++)
    {
      inserted=false;
      if(hitvec2.at(m).GetHitNumPad()>0)
	{
	  if((hitvec2.at(m).GetHitPhi()<20)&&(hitvec2.at(m).GetHitPhi()>0))
	    numphimin20++;
	
	  if((hitvec2.at(m).GetHitPhi()>340))
	    numphimag340++; 
	}
   
      
      if(hitvec2.at(m).GetHitNumPad()==0){
	if(Usestrip==true){
	  inserted=true ;
	  hitvec.push_back(hitvec2.at(m));
	}
      }
      else
	{
	  inserted=true ;
	  hitvec.push_back(hitvec2.at(m));
	}


      if(inserted){

	if(verbosity)
	  std::cout<<"Z-Phi: "<<hitvec2.at(m).GetHitZ()<<" "<<hitvec2.at(m).GetHitPhi()<<" Num Pad-Strip:"<<hitvec2.at(m).GetHitNumPad()<<" "<<hitvec2.at(m).GetHitNumStrip()<<std::endl;


	numHit_t2++;
	if((hitvec2.at(m).GetHitNumPad()>0)&&(hitvec2.at(m).GetHitNumStrip()==0))
	  numPadOnlyHit_t2++;

	if((hitvec2.at(m).GetHitNumPad()==0)&&(hitvec2.at(m).GetHitNumStrip()>0))
	  numStripOnlyHit_t2++;

	if((hitvec2.at(m).GetHitNumPad()>0)&&(hitvec2.at(m).GetHitNumStrip()>0))
	  {
	    numCl1HitHit_t2++;
	   hitvec2.at(m).GetHitPhi();
	    hitvec2.at(m).GetHitDPhi();
	  }

      }
      
    }

 
  /////////////////////////////////////////////////////////////////////////////////
  ///           Variable INITIALIZATION
  ///////////////////////////////////////////////////////////////////////////////


 int hemisphere = 0;

  unsigned int  sizeHitv=hitvec.size();
  if(sizeHitv<2)
   {
     std::cout<<" T2SelectionCutUtils::MyLinearfitCorr problem: Track with less than 2 Cl1 hits!! Continue the fitting.."<<std::endl;        }
  else
    hemisphere = (int)(hitvec[0].GetHitZ()/fabs(hitvec[0].GetHitZ()));

  
  TMatrixD ParCov_matr(4,4); //matrice di covarianza dei parametri trovati;
  unsigned int sizeArighe=sizeHitv*2;
  TMatrixD A(sizeArighe,4);
  TMatrixD At(4,sizeArighe);
  //A è la matrice per cui Mis=A(Param)
  TMatrixD Vy(sizeArighe,sizeArighe); //matrice di covarianza delle misure (una per ogni xy, quindi è diag a blocchi);
  TMatrixD Vym1(sizeArighe,sizeArighe); 

  TMatrixD Ls(4,4);//((A^T)(Vy^-1)A)
  // TMatrixD Cs(4,4);//(A^T)(Vy^-1)
  TMatrixD Cs(4,sizeArighe);//(A^T)(Vy^-1)

  TVectorD FittedParam(4);
 
  TVectorD Yvect(sizeArighe);


  Ls.Zero();
  TMatrixD Ls00(2,2);//un quarto della matrice LS.
  TMatrixD Ls01(2,2);
  TMatrixD Ls10(2,2);
  TMatrixD Ls11(2,2);

  TMatrixD Mia(2,2);
  TMatrixD Mib(2,2);

  TMatrixD MiaT(2,2);
  TMatrixD MibT(2,2);

  Ls00.Zero();
  Ls01.Zero();
  Ls10.Zero();
  Ls11.Zero();

  Vy.Zero();
  Vym1.Zero();

  if(verbosity)
    std::cout<<"MyLinearfitCorr Start computation  .. "<<std::endl;
  
  std::vector<std::vector<double> > All_YMeasures;
  std::vector<std::vector<double> > All_Vym1Measures;//11 12 21 22
  std::vector<std::vector<double> > All_VyMeasures;
  std::vector<double> All_ZMeasures; 





  for(unsigned int k =0; k<hitvec.size(); k++)
    {
      TMatrixD OneVy(2,2); 
      OneVy.Zero();
      
      phirad=hitvec[k].GetHitPhi()*3.14159/180.0;  
      sigmaPhi=hitvec[k].GetHitDPhi()*3.14159/180.0;
      
      r=hitvec[k].GetHitR();
      sigmaR=hitvec[k].GetHitDR();

      if(hitvec[k].GetHitNumStrip()>4)
	sigmaR=(hitvec[k].GetHitDR()+(hitvec[k].GetHitNumStrip()-4)*0.4);

 
   /////////////////////////////////////////////////////////////////////////////////
   ///         FITTING PART
   ///////////////////////////////////////////////////////////////////////////////



   OneVy(0,0)=sigmaR*sigmaR*cos(phirad)*cos(phirad)+r*r*sin(phirad)*sin(phirad)*sigmaPhi*sigmaPhi;//ex
   OneVy(0,1)=cos(phirad)*sin(phirad)*(sigmaR*sigmaR-r*r*sigmaPhi*sigmaPhi);  
   OneVy(1,0)=OneVy(0,1);
   OneVy(1,1)=sigmaR*sigmaR*sin(phirad)*sin(phirad)+ r*r*cos(phirad)*cos(phirad)*sigmaPhi*sigmaPhi;//ey
   
   //Invert Vy matrix
   TMatrixD OneVym1(2,2); 
   OneVym1=OneVy;
   Double_t deti;	
   OneVym1.Invert(&deti);
   
   //OneVy.Print();
   if(fabs(deti)<0.001)
     {
       std::cout<<"Possible Vy Zero Determinant in one point error matrix!!!"<<std::endl;
       std::cout<<"sigmaR-R-phi-SigmaPhi:"<<sigmaR<<" "<<r<<" "<<phirad*180.0/3.14159<<" "<<sigmaPhi*180.0/3.14159<<" NumStrip: "<<hitvec[k].GetHitNumStrip()<<" NumPad:"<<hitvec[k].GetHitNumPad()<<std::endl;
     }

   Vym1(2*k,2*k)=OneVym1(0,0);
   Vym1(2*k,2*k+1)= OneVym1(0,1); 
   Vym1(2*k+1,2*k)=OneVym1(1,0);
   Vym1(2*k+1,2*k+1)=OneVym1(1,1);
   

   Vy(2*k,2*k)=OneVy(0,0);
   Vy(2*k,2*k+1)= OneVy(0,1); 
   Vy(2*k+1,2*k)=OneVy(1,0);
   Vy(2*k+1,2*k+1)=OneVy(1,1);

   
   std::vector<double>recYvect; 
   recYvect.push_back(hitvec[k].GetHitX());
   recYvect.push_back(hitvec[k].GetHitY());

   std::vector<double>  recVy;
   recVy.push_back(OneVy(0,0)); recVy.push_back(OneVy(0,1)); 
   recVy.push_back(OneVy(1,0)); recVy.push_back(OneVy(1,1));

   std::vector<double>  recVym1; 
   recVym1.push_back(OneVym1(0,0)); recVym1.push_back(OneVym1(0,1)); 
   recVym1.push_back(OneVym1(1,0)); recVym1.push_back(OneVym1(1,1));
   // std::vector<double> recZvect; recZvect.push_back(hitvec[k].GetHitZ());

   All_YMeasures.push_back(recYvect);
   All_Vym1Measures.push_back(recVym1);
   All_VyMeasures.push_back(recVy);
   All_ZMeasures.push_back(hitvec[k].GetHitZ());

   
   



   Yvect(2*k)=hitvec[k].GetHitX();
   Yvect(2*k+1)=hitvec[k].GetHitY();  

   A(2*k,0)=hitvec[k].GetHitZ();
   A(2*k,1)=1.0;
   A(2*k,2)=0.0;
   A(2*k,3)=0.0;
   A(2*k+1,0)=0.0;
   A(2*k+1,1)=0.0;
   A(2*k+1,2)=hitvec[k].GetHitZ();
   A(2*k+1,3)=1.0;

     

   Mia.Zero();
   Mib.Zero();
   
   
   Mia(0,0)=hitvec[k].GetHitZ();
   Mia(0,1)=0.;
   Mia(1,0)=1.;
   Mia(1,1)=0.;
   Mib(0,0)=0.;
   Mib(0,1)=hitvec[k].GetHitZ();
   Mib(1,0)=0.;
   Mib(1,1)=1.;
   
   MiaT=Mia;
   MibT=Mib;
   MiaT.Transpose(MiaT);
   MibT.Transpose(MibT);  
      

   Ls00=Ls00+Mia*OneVym1*MiaT;
   
   Ls01=Ls01+Mia*OneVym1*MibT;
   
   Ls10=Ls10+Mib*OneVym1*MiaT;
   
   Ls11=Ls11+Mib*OneVym1*MibT;  


 } //end loop on hits.

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

 if(verbosity)
   std::cout<<"MyLinearfitCorr Parameter calculation  .. "<<std::endl;

//Nota che qui Vy è in realtà Vy^-1 

 Ls(0,0)=Ls00(0,0);  Ls(0,1)=Ls00(0,1);   Ls(0,2)=Ls01(0,0);   Ls(0,3)=Ls01(0,1);
 Ls(1,0)=Ls00(1,0);  Ls(1,1)=Ls00(1,1);   Ls(1,2)=Ls01(1,0);   Ls(1,3)=Ls01(1,1); 

 Ls(2,0)=Ls10(0,0);  Ls(2,1)=Ls10(0,1);   Ls(2,2)=Ls11(0,0);   Ls(2,3)=Ls11(0,1);
 Ls(3,0)=Ls10(1,0);  Ls(3,1)=Ls10(1,1);   Ls(3,2)=Ls11(1,0);   Ls(3,3)=Ls11(1,1);

 
 Double_t determ;	
 Ls.Invert(&determ);
 if(fabs(determ)<0.001)
     std::cout<<"WARNING: Possible Zero LS  Determinant!!!"<<std::endl;

 At.Zero();
 At.Transpose(A);
 Cs.Zero();
 Cs=At*Vym1;


 
 FittedParam= Ls*Cs*Yvect;//ax,bx.ay.by..
 
 //std::cout<<"MyLinearfitCorr FittedParam: "<<std::endl;
 //FittedParam.Print();

 ParCov_matr=At*Vym1*A;


 ParCov_matr.Invert(&determ);
 //std::cout<<"MyLinearfitCorr Error Matrix: "<<std::endl;
 //ParCov_matr.Print();

 if(fabs(determ)<0.001)
   std::cout<<"WARNING: Possible Zero ParCov_matr  Determinant!!!"<<std::endl;  

 //Calcolo chi2
 double chi2=0.;
 TVectorD Residui(sizeArighe);
 Residui=(Yvect-A*FittedParam);
 
 //TVectorD ResiduiT;
 //ResiduiT.Transpose(Residui);
 TVectorD VdotRes(sizeArighe);
 VdotRes=(Vym1*Residui);
 chi2=Residui*VdotRes;

 
 
 //////////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////// SAVE THE WORST T2-HIT , Calculate single projection chi2X,chi2Y.
 /////////////////////////////////////////////////////////////////////////////////////

 // std::vector<std::vector<double> > All_YMeasures;
 //std::vector<std::vector<double> > All_Vym1Measures;//11 12 21 22
 //std::vector<double> All_ZMeasures; 

 double chi2X=0.; double chi2Y=0.; 
 double maxchi2contr=0.;double worsthitz=hitvec[0].GetHitZ();
 for(unsigned int i=0;i<All_YMeasures.size();i++)
   {     
     double chi2contr=0.;
     TVectorD Misura(2);   Misura(0)=All_YMeasures.at(i).at(0);  Misura(1)=All_YMeasures.at(i).at(1);

     TMatrixD OneVym1(2,2);  
     OneVym1(0,0)= All_Vym1Measures.at(i).at(0); 
     OneVym1(0,1)= All_Vym1Measures.at(i).at(1);
     OneVym1(1,0)= All_Vym1Measures.at(i).at(2); 
     OneVym1(1,1)= All_Vym1Measures.at(i).at(3);
     
     TMatrixD OneVy(2,2);  
     OneVy(0,0)= All_VyMeasures.at(i).at(0); 
     OneVy(0,1)= All_VyMeasures.at(i).at(1);
     OneVy(1,0)= All_VyMeasures.at(i).at(2); 
     OneVy(1,1)= All_VyMeasures.at(i).at(3);
     

     TVectorD Residuo(2);

     Residuo(0)=FittedParam(0)*All_ZMeasures.at(i)+FittedParam(1);
     Residuo(1)=FittedParam(2)*All_ZMeasures.at(i)+FittedParam(3);
     

     Residuo=Misura-Residuo;//Bug found 25/1/2011. Subtraction missing.
     chi2X=chi2X+Residuo(0)*Residuo(0)/OneVy(0,0);
     chi2Y=chi2Y+Residuo(1)*Residuo(1)/OneVy(1,1);

     
     TVectorD Temp(2);
     Temp=OneVym1*Residuo;
     chi2contr=Residuo*Temp;

     //OneVym1(2,2)=;OneVym1(2,2)=;OneVym1(2,2)=;OneVym1(2,2)=;
     if(chi2contr>maxchi2contr){
       maxchi2contr=chi2contr;
       worsthitz=All_ZMeasures.at(i);
     }
   }

 for(unsigned int u=0;u<hitvec.size();u++)
   {
     if(fabs(hitvec.at(u).GetHitZ()-worsthitz)<0.5)
       worstHit=hitvec.at(u);
   }


 //////////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////// SAVE THE OUTPUT
 /////////////////////////////////////////////////////////////////////////////////////


 
 if(chi2<0)
   {
     std::cout<<"WARNING: MyLinearfitCorr Chi2: "<<chi2<<std::endl;
     std::cout<<"Residui"<<std::endl;
     Residui.Print();
     std::cout<<"Cov Matrix Vy"<<std::endl;
     Vy.Print();
     std::cout<<"Vym1*Residui"<<std::endl;
     VdotRes.Print();

   }
 par_covariance_matrix=ParCov_matr;
 
 //std::cout<<"MyLinearfitCorr Error Matrix: "<<std::endl;
 //par_covariance_matrix.Print();


 std::vector<double> vect;

 for (unsigned int i=0;i<10;i++)
   vect.push_back(0.);

 vect[0]=FittedParam(0);
 vect[1]=FittedParam(1);
 vect[2]=FittedParam(2);
 vect[3]=FittedParam(3);
 vect[4]=ParCov_matr(0,0);
 vect[5]=ParCov_matr(1,1);
 vect[6]=ParCov_matr(2,2);
 vect[7]=ParCov_matr(3,3);

 double correlabX=ParCov_matr(0,1);
 vect[8]=correlabX;
 double correlabY=ParCov_matr(2,3);
 vect[9]=correlabY;


 //Swap of positions in order to be compatible with the track object.
 TMatrixD ParCov_matrConv(4,4); 
 TVectorD FittedParamConv(4);
 

 for(unsigned int yy=0;yy<4;yy++)
   for(unsigned int ll=0;ll<4;ll++)
     ParCov_matrConv(yy,ll)=0.;


 ParCov_matrConv(0,0)=ParCov_matr(1,1);
 ParCov_matrConv(1,1)=ParCov_matr(3,3);
 ParCov_matrConv(2,2)=ParCov_matr(0,0);
 ParCov_matrConv(3,3)=ParCov_matr(2,2);

 FittedParamConv(0)=FittedParam(1);
 FittedParamConv(1)=FittedParam(3);
 FittedParamConv(2)=FittedParam(0);
 FittedParamConv(3)=FittedParam(2);


 // par_covariance_matrix=ParCov_matr;

 /*
   FittedParam(0)=b_xz;
   FittedParam(1)=b_yz;
   FittedParam(2)=a_xz;
   FittedParam(3)=a_yz;


 TMatrixD ParCov_matr(4,4); 
 for(unsigned int uu=0;uu<4;uu++)
   for(unsigned int ii=0;ii<4;ii++)
     ParCov_matr(uu,ii)=0.;


     ParCov_matr(0,0)=e_b_xz;
     ParCov_matr(1,1)=e_b_yz;
     ParCov_matr(2,2)=e_a_xz;
     ParCov_matr(3,3)=e_a_yz;
 */


  
  
 /*
  if((fabs(vect[0])>1.8)||(fabs(vect[2])>1.8))
    {
      std::cout<<"WARNING in T2TrackProducer3 Ax:"<<(vect[0])<<"  Ay:"<<(vect[2])<<" too big!! Trk Hit (all) print :"<<std::endl;
      for(unsigned int jj =0; jj<hitvec.size(); jj++)
	{
	  std::cout<<"x-y-z-phi:  "<<hitvec[jj].GetHitX()<<"  "<<hitvec[jj].GetHitY()<<"  "<<hitvec[jj].GetHitZ()<<"  "<<hitvec[jj].GetHitPhi()<<"     DX:"<<hitvec[jj].GetHitDX()<<"     DY:"<<hitvec[jj].GetHitDY()<<" NumPad:"<<hitvec[jj].GetHitNumPad()<<std::endl;
	}
    }

  if(verbosity){
    std::cout<<"Before save: XY:"<<std::endl;
    // std::cout<<"ax.bx.ay.by:"<<vect[0]<<" "<<vect[1]<<" "<<vect[2]<<" "<<vect[3]<<std::endl;
  
    std::cout<<"Before save: RZ:"<<std::endl;
    std::cout<<"Phi "<<phiRZ*180.0/3.14159<<" TanTheta"<<TanthetaRZ<<" Brz:"<<bRZ<<" Eta:"<<trk2rz.Eta()<<std::endl;

    std::cout<<" TotHit:"<<numHit_t2<<" cl1Hit:"<<numCl1HitHit_t2<<" StripOnly:"<<numStripOnlyHit_t2<<" PadOnly:"<<numPadOnlyHit_t2<<std::endl;
  }
 */
  double chi2Phi=0.;
  double chi2R=0.;
  
  double phiRZ=0.;               //trk2rz.Phi(); 
  double e_phiRZ=0.;
  
  double TanthetaRZ=0.;     //trk2rz.GetTy(); 
  double e_TanthetaRZ=0.; //trk2rz.GetTySigma();
  
  double bRZ=0.;          //trk2rz.GetTx(); 
  double e_bRZ=0.;      //trk2rz.GetTxSigma();
  


  T1T2Track fittedtrack(FittedParamConv,ParCov_matrConv,chi2,chi2X,chi2Y,hemisphere,2 , TanthetaRZ,  bRZ,   phiRZ,   e_TanthetaRZ,  e_bRZ,  e_phiRZ, chi2R, chi2Phi, RoadID, numHit_t2, numCl1HitHit_t2, numStripOnlyHit_t2, numPadOnlyHit_t2);
  

  if(verbosity){      
    std::cout<<"After save: :"<<std::endl;
    std::cout<<"PhiRZ:"<<fittedtrack._phiRZ*180.0/3.14159<<"  EtaRZ:"<< fittedtrack._etaRZ<<"  EtaXY:"<<fittedtrack.Eta()<<" ChirR/N:"<<(chi2R/(numCl1HitHit_t2+numPadOnlyHit_t2))<<" Chi2/N; "<<chi2/(numCl1HitHit_t2+numPadOnlyHit_t2)<<" NumHit:"<<fittedtrack._numHit_t2<<std::endl;
    
  }


  sizeHitv=hitvec.size();
  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      fittedtrack.AddHit(hitvec[jj]);
    }
  
 // std::cout<<"Trk with "<<sizeHitv<<" hit saved"<<std::endl;
  return fittedtrack;
  
}





T1T2Track T2RoadPadFinder::fitTrackspecialXYDEV(std::vector<T2Hit> &RoadHits)
{

  unsigned int realcandidatesize=RoadHits.size();

  T1T2Track trk;
  trk.SetDet(2);					

  for(unsigned int l=0;l<realcandidatesize;l++)
    {	 
      trk.AddHit(RoadHits.at(l));	
    }      
  
  
  unsigned int sizeHitv=trk.GetHitEntries();   
  std::vector<T2Hit> hitvec;
  
  for (unsigned int l=0;l<sizeHitv;l++)
    {
      hitvec.push_back(trk.GetHitT2(l));     	   
    }  
  
  
  //std::cout<<" Here 1 "<<std::endl;
  sizeHitv=hitvec.size();
  
  int hemisphere=0; 
   
   if(sizeHitv>=2)
     {
       hemisphere = (int)(hitvec[1].GetHitZ()/fabs(hitvec[1].GetHitZ()));
     }
   else
     {
        hemisphere = (int)(hitvec[0].GetHitZ()/fabs(hitvec[0].GetHitZ()));
	std::cout<<"ERROR in T2TrackProducer2::fitTrackspecialXY: using vtx position to determine the hemisphere "<<std::endl;
     }


double Sx=0.;
double Sxz=0.;
double Szz_x=0.;
double Sz_x=0.; 
double S0_x=0.; 

double Sy=0.;
double Syz=0.;
double Szz_y=0.;
double Sz_y=0.; 
double S0_y=0.; 

std::vector<double> x;   
std::vector<double> z;
std::vector<double> y;
std::vector<double> ex;
std::vector<double> ez;
std::vector<double> ey;
 double a_xz, b_xz, a_yz, b_yz;

  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      //   std::cout<<hitvec[jj].GetHitX()<<" "<<hitvec[jj].GetHitPhi()<<" "<<hitvec[jj].GetHitZ()<<std::endl;
      //r->x   phi->y
      
      x.push_back(hitvec[jj].GetHitX());
      y.push_back(hitvec[jj].GetHitY());
      z.push_back(hitvec[jj].GetHitZ());    
 
      double phirad=hitvec[jj].GetHitPhi()*3.14159/180.0;
      double sigmay;
      double sigmax;
      //er.push_back(0.1);
      if(hitvec[jj].GetHitClass()==1)
	{
	  sigmax=cos(phirad)*cos(phirad)*0.12*0.12+sin(phirad)*sin(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
	  ex.push_back(sqrt(sigmax));
      
	  //ephi.push_back(0.1);
	  sigmay=sin(phirad)*sin(phirad)*0.12*0.12+cos(phirad)*cos(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
	  ey.push_back(sqrt(sigmay));
	}
      else   //To be tuned properly !!!!!!!!!!!!!
	{
	  double dr=hitvec[jj].GetHitDR();

	  if(hitvec[jj].GetHitClass()!=9)//9=THE VTX
	    {
	      sigmax=cos(phirad)*cos(phirad)*dr*dr+sin(phirad)*sin(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
	      sigmay=sin(phirad)*sin(phirad)*dr*dr+cos(phirad)*cos(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
	    }
	  else
	    {	      
	       sigmax=hitvec[jj].GetHitDX()*hitvec[jj].GetHitDX();
	       sigmay=hitvec[jj].GetHitDY()*hitvec[jj].GetHitDY();		
	    }
	 
	  ex.push_back(sqrt(sigmax));      
	  //ephi.push_back(0.1);
	  
	  ey.push_back(sqrt(sigmay));

	}
      if((sigmax==0)||(sigmay==0))
	std::cout<<"ERROR in T2TrackProducer2::fitTrackspecialXY: sigmaX=  "<<sigmax<<" sigmaY= "<<sigmay<<std::endl;
      //ez.push_back(hitvec[jj].GetHitDZ());
            
      Sxz += x[jj]*z[jj]/ex[jj]/ex[jj];
      Szz_x += z[jj]*z[jj]/ex[jj]/ex[jj];
      Sz_x += z[jj]/ex[jj]/ex[jj];
      Sx += x[jj]/ex[jj]/ex[jj];
      S0_x += 1.0/ex[jj]/ex[jj];



      Syz += y[jj]*z[jj]/ey[jj]/ey[jj];
      Szz_y += z[jj]*z[jj]/ey[jj]/ey[jj];
      Sz_y += z[jj]/ey[jj]/ey[jj];
      Sy += y[jj]/ey[jj]/ey[jj];
      S0_y += 1.0/ey[jj]/ey[jj];

  }


  //std::cout<<" Here 2 "<<std::endl;

a_xz = (Sxz*S0_x - Sz_x*Sx) / (Szz_x*S0_x - Sz_x*Sz_x);   // angular coefficient
b_xz = (Sx*Szz_x - Sz_x*Sxz) / (Szz_x*S0_x - Sz_x*Sz_x);  // intercept   X=(a_xz)Z + b_xz

a_yz = (Syz*S0_y - Sz_y*Sy) / (Szz_y*S0_y - Sz_y*Sz_y);   // angular coefficient
b_yz = (Sy*Szz_y - Sz_y*Syz) / (Szz_y*S0_y - Sz_y*Sz_y);  // intercept  Y=(a_yz)Z + b_yz 

double e_a_xz = sqrt( S0_x / (S0_x*Szz_x - Sz_x*Sz_x) );
double e_b_xz = sqrt( Szz_x / (S0_x*Szz_x - Sz_x*Sz_x) );
double e_a_yz = sqrt( S0_y / (S0_y*Szz_y - Sz_y*Sz_y) );
double e_b_yz = sqrt( Szz_y / (S0_y*Szz_y - Sz_y*Sz_y) );

//USE T1 Convenction
TVectorD vect(4);
 for(int oo=0; oo<4; oo++)
   vect[oo]=0.;

 vect[0] = b_xz;
 vect[1] = b_yz;
 vect[2] = a_xz;
 vect[3] = a_yz;
 
//double correlx=(-1.0)*Sz_x*(1.0/(S0_x*Szz_x-(Sz_x*Sz_x)));
//double correly=(-1.0)*Sz_y*(1.0/(S0_y*Szz_y-(Sz_y*Sz_y)));

 

 TMatrixD mat(4,4);
 //Save time for this purpose.
 
 for(unsigned int oo=0; oo<4; oo++)
   for(unsigned int ooo=0;ooo<4; ooo++)
     mat[oo][ooo]=0.;


     mat[0][0] = e_b_xz*e_b_xz;
    mat[1][1] = e_b_yz*e_b_yz;
    mat[2][2] = e_a_xz*e_a_xz;
    mat[3][3] = e_a_yz*e_a_yz;
 

    double chi2 = 0;
    double chi2X = 0;
    double chi2Y = 0;

    for(unsigned int jjj =0; jjj<sizeHitv; jjj++){
      
      chi2X += (a_xz*z[jjj]+b_xz - x[jjj])*(a_xz*z[jjj]+b_xz - x[jjj])/ex[jjj]/ex[jjj];
      chi2Y += (a_yz*z[jjj]+b_yz - y[jjj])*(a_yz*z[jjj]+b_yz - y[jjj])/ey[jjj]/ey[jjj];

      TMatrixD OneVy(2,2); 
      OneVy.Zero();
      double phirad=trk.GetHitT2(jjj).GetHitPhi()*3.14159/180.0;
      double r=trk.GetHitT2(jjj).GetHitR();  
      double sigmaR=trk.GetHitT2(jjj).GetHitDR();double sigmaPhi=trk.GetHitT2(jjj).GetHitDPhi()*3.14159/180.0;
      
      OneVy(0,0)=sigmaR*sigmaR*cos(phirad)*cos(phirad)+r*r*sin(phirad)*sin(phirad)*sigmaPhi*sigmaPhi;
      OneVy(0,1)=cos(phirad)*sin(phirad)*(sigmaR*sigmaR-r*r*sigmaPhi*sigmaPhi);  
      OneVy(1,0)=OneVy(0,1);
      OneVy(1,1)=sigmaR*sigmaR*sin(phirad)*sin(phirad)+ r*r*cos(phirad)*cos(phirad)*sigmaPhi*sigmaPhi;
      
      //Invert Vy matrix
      TMatrixD OneVym1(2,2); 
      OneVym1=OneVy;
      Double_t deti;	
      OneVym1.Invert(&deti);
      
      TVectorD Residui(2);
      Residui(0)=(a_xz*z[jjj]+b_xz - x[jjj]);
      Residui(1)=(a_yz*z[jjj]+b_yz - y[jjj]);
      
      TVectorD VdotRes(2);
      VdotRes=(OneVym1*Residui);
      chi2=chi2+Residui*VdotRes;
    }


    std::cout<<chi2<<" BX-BY-AX-AY:"<< vect[0]<<" "<<vect[1]<<" "<<vect[2]<<" "<<vect[3]<<std::endl;

 T1T2Track fittedtrack(vect,mat,chi2,chi2X,chi2Y,hemisphere,2);

 sizeHitv=hitvec.size();
 for(unsigned int jj =0; jj<sizeHitv; jj++)
   {
     fittedtrack.AddHit(hitvec[jj]);
   }
 
   return fittedtrack; 
}






T1T2Track T2RoadPadFinder::fitTrackspecialXY(std::vector<T2Hit> &hitvec)
{
   

   
   unsigned int sizeHitv=hitvec.size();   
  
   

   int hemisphere=0; 




   if(sizeHitv>=3)
     {
       hemisphere = (int)(hitvec[1].GetHitZ()/fabs(hitvec[1].GetHitZ()));
     }
   else
     {
       if(sizeHitv>=1)
	 hemisphere = (int)(hitvec[0].GetHitZ()/fabs(hitvec[0].GetHitZ()));
       std::cout<<"ERROR in T2TrackProducer2::fitTrackspecialXY: size:"<<hitvec.size()<<std::endl;
     }


double Sx=0.;
double Sxz=0.;
double Szz_x=0.;
double Sz_x=0.; 
double S0_x=0.; 

double Sy=0.;
double Syz=0.;
double Szz_y=0.;
double Sz_y=0.; 
double S0_y=0.; 

std::vector<double> x;   
std::vector<double> z;
std::vector<double> y;
std::vector<double> ex;
std::vector<double> ez;
std::vector<double> ey;
 double a_xz, b_xz, a_yz, b_yz;

  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      //   std::cout<<hitvec[jj].GetHitX()<<" "<<hitvec[jj].GetHitPhi()<<" "<<hitvec[jj].GetHitZ()<<std::endl;
      //r->x   phi->y
      
      x.push_back(hitvec[jj].GetHitX());
      y.push_back(hitvec[jj].GetHitY());
      z.push_back(hitvec[jj].GetHitZ());    
 
      double phirad=hitvec[jj].GetHitPhi()*3.14159/180.0;
      double sigmay;
      double sigmax;
      //er.push_back(0.1);
      if(hitvec[jj].GetHitClass()==1)
	{
	  sigmax=cos(phirad)*cos(phirad)*0.12*0.12+sin(phirad)*sin(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
	  ex.push_back(sqrt(sigmax));
      
	  //ephi.push_back(0.1);
	  sigmay=sin(phirad)*sin(phirad)*0.12*0.12+cos(phirad)*cos(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
	  ey.push_back(sqrt(sigmay));
	}
      else   //To be tuned properly !!!!!!!!!!!!!
	{
	  double dr=hitvec[jj].GetHitDR();

	  if(hitvec[jj].GetHitClass()!=9)//9=THE VTX
	    {
	      sigmax=cos(phirad)*cos(phirad)*dr*dr+sin(phirad)*sin(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
	      sigmay=sin(phirad)*sin(phirad)*dr*dr+cos(phirad)*cos(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
	    }
	  else
	    {	      
	       sigmax=hitvec[jj].GetHitDX()*hitvec[jj].GetHitDX();
	       sigmay=hitvec[jj].GetHitDY()*hitvec[jj].GetHitDY();		
	    }
	 
	  ex.push_back(sqrt(sigmax));      
	  //ephi.push_back(0.1);
	  
	  ey.push_back(sqrt(sigmay));

	}
      if((sigmax==0)||(sigmay==0))
	std::cout<<"ERROR in T2TrackProducer2::fitTrackspecialXY: sigmaX=  "<<sigmax<<" sigmaY= "<<sigmay<<std::endl;
      //ez.push_back(hitvec[jj].GetHitDZ());
            
      Sxz += x[jj]*z[jj]/ex[jj]/ex[jj];
      Szz_x += z[jj]*z[jj]/ex[jj]/ex[jj];
      Sz_x += z[jj]/ex[jj]/ex[jj];
      Sx += x[jj]/ex[jj]/ex[jj];
      S0_x += 1.0/ex[jj]/ex[jj];



      Syz += y[jj]*z[jj]/ey[jj]/ey[jj];
      Szz_y += z[jj]*z[jj]/ey[jj]/ey[jj];
      Sz_y += z[jj]/ey[jj]/ey[jj];
      Sy += y[jj]/ey[jj]/ey[jj];
      S0_y += 1.0/ey[jj]/ey[jj];

  }


  //std::cout<<" Here 2 "<<std::endl;

a_xz = (Sxz*S0_x - Sz_x*Sx) / (Szz_x*S0_x - Sz_x*Sz_x);   // angular coefficient
b_xz = (Sx*Szz_x - Sz_x*Sxz) / (Szz_x*S0_x - Sz_x*Sz_x);  // intercept   X=(a_xz)Z + b_xz

a_yz = (Syz*S0_y - Sz_y*Sy) / (Szz_y*S0_y - Sz_y*Sz_y);   // angular coefficient
b_yz = (Sy*Szz_y - Sz_y*Syz) / (Szz_y*S0_y - Sz_y*Sz_y);  // intercept  Y=(a_yz)Z + b_yz 

double e_a_xz = sqrt( S0_x / (S0_x*Szz_x - Sz_x*Sz_x) );
double e_b_xz = sqrt( Szz_x / (S0_x*Szz_x - Sz_x*Sz_x) );
double e_a_yz = sqrt( S0_y / (S0_y*Szz_y - Sz_y*Sz_y) );
double e_b_yz = sqrt( Szz_y / (S0_y*Szz_y - Sz_y*Sz_y) );

//USE T1 Convenction
TVectorD vect(4);
 for(int oo=0; oo<4; oo++)
   vect[oo]=0.;

 vect[0] = b_xz;
 vect[1] = b_yz;
 vect[2] = a_xz;
 vect[3] = a_yz;
 
//double correlx=(-1.0)*Sz_x*(1.0/(S0_x*Szz_x-(Sz_x*Sz_x)));
//double correly=(-1.0)*Sz_y*(1.0/(S0_y*Szz_y-(Sz_y*Sz_y)));

 

 TMatrixD mat(4,4);
 for(unsigned int oo=0; oo<4; oo++)
   for(unsigned int ooo=0;ooo<4; ooo++)
     mat[oo][ooo]=0.;


 mat[0][0] = e_b_xz*e_b_xz;
    mat[1][1] = e_b_yz*e_b_yz;
    mat[2][2] = e_a_xz*e_a_xz;
    mat[3][3] = e_a_yz*e_a_yz;

    double chi2 = 0;
    double chi2X = 0;
    double chi2Y = 0;

for(unsigned int jjj =0; jjj<sizeHitv; jjj++)
    {
      chi2X += (a_xz*z[jjj]+b_xz - x[jjj])*(a_xz*z[jjj]+b_xz - x[jjj])/ex[jjj]/ex[jjj];
      chi2Y += (a_yz*z[jjj]+b_yz - y[jjj])*(a_yz*z[jjj]+b_yz - y[jjj])/ey[jjj]/ey[jjj];
    }

 chi2 = (chi2X+chi2Y)/2.0; 


 // std::cout<<"chi2X chi2Y "<<chi2X<<" "<<chi2Y<<std::endl;

T1T2Track fittedtrack(vect,mat,chi2,chi2X,chi2Y,hemisphere,2);
//std::cout<<"chi2X after: "<<fittedtrack.ChiSquaredX()<<"  chi2Y after: "<<fittedtrack.ChiSquaredY()<<std::endl;
//sizeHitv=trk.GetHitEntries();
 sizeHitv=hitvec.size();
for(unsigned int jj =0; jj<sizeHitv; jj++)
  {
    fittedtrack.AddHit(hitvec[jj]);
  }
 //fittedtrack.AddHit(trk.GetHitT2(jj));
  // fittedtrack.AddHit(hitvec[jj]);

//  std::cout<<" Here 4 "<<std::endl;

   return fittedtrack; 
}





/*

std::vector<std::vector <long> >  T2RoadPadFinder::NoFixPoint_HitsCombination(std::map<double, std::vector<long> > &ZvsUniqueHitId_inroad)
{

  // std::vector<std::vector <long> > cl1Hitcombiation;
  //cl1Hitcombiation.clear();

  //This is the same strategy as in trackproducer2. Just changing T2Hits with its id (long)
  // std::vector<long> sameZhits;  
  std::vector<std::vector<long> > hitmatrix;  hitmatrix.clear();
  std::vector<std::vector<long> > allcandidate;  allcandidate.clear();
  long lasthit;
  
  std::vector<long> candidatepiece; candidatepiece.clear();
  
  std::map<double, std::vector<long> >::const_iterator ZvsUniqueHitId_inroad_iter;

  // Hit Matrix construction with in same row->same Z hit
  for (ZvsUniqueHitId_inroad_iter=ZvsUniqueHitId_inroad.begin();ZvsUniqueHitId_inroad_iter!=ZvsUniqueHitId_inroad.end();ZvsUniqueHitId_inroad_iter++)  
    {      
      hitmatrix.push_back((*ZvsUniqueHitId_inroad_iter).second);      
    }
  

  allcandidate.clear();                 //allcandidate is a std::vector<std::vector<T2Hit> > ; vector of track-hits_ID
  for(unsigned int k=0;k<hitmatrix.size();k++)
    {
      //std::vector<long> matrixraw=hitmatrix.at(k);
      
      if(allcandidate.size()==0){ // fill allcandidate with the first row of hitmatrix
	for(unsigned int l=0;l<hitmatrix.at(k).size();l++)
	  {
	    candidatepiece.push_back(hitmatrix.at(k).at(l));
	    allcandidate.push_back(candidatepiece);
	    candidatepiece.clear();
	  }
      }
      else
	{
	  unsigned int currentsz=allcandidate.size();
	  unsigned int m=0;
	  std::vector<std::vector<long> > newcandidatestoput;
	  
	  currentsz=allcandidate.size();
	  
	  while(m<currentsz)
	    {
	      //join to all the candidates each element of the full next-Z row 
	      //so if the input was 2 candidate and the next Z have 3 elements at the end you have 2*3 candidates, each one differs for a next Z value.
	      candidatepiece=allcandidate.at(m);			 		       			 		     		     
	      for(unsigned int l=0;l<hitmatrix.at(k).size();l++)
		{
		  std::vector<long> newcandidatepiece=candidatepiece;
		  newcandidatepiece.push_back(hitmatrix.at(k).at(l));
		  //allcandidate.push_back(newcandidatepiece);
		  newcandidatestoput.push_back(newcandidatepiece);
		  newcandidatepiece.clear();
		  
		}
	      
	      candidatepiece.clear();
	      m++;
	    }
	  
	  
	  allcandidate=newcandidatestoput;
	  newcandidatestoput.clear();		     
	  
	}	               
      
    }

 

  return allcandidate;
}


*/

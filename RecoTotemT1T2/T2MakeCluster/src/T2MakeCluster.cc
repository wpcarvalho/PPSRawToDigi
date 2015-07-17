/**
 * Class T2MakeCluster
 *
 * Author: Mirko Berretti / University of Siena 
 * Email:  mirko.berretti@gmail.com
 * Date:   2007-12-08
 */
 

#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "RecoTotemT1T2/T2MakeCluster/interface/T2MakeCluster.h"
#include "RecoTotemT1T2/T2MakeCluster/interface/T2DetClustReconst.h"

#include "DataFormats/T2Digi/interface/T2StripDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2PadDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2StripDigi.h"
#include "DataFormats/T2Digi/interface/T2PadDigi.h"
#include "DataFormats/T2Cluster/interface/T2StripClusterCollection.h"
#include "DataFormats/T2Cluster/interface/T2PadClusterCollection.h"
#include <iostream>
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"
#include "FWCore/Framework/interface/ESHandle.h"

//clustering id: int ClustId_unique in the event
//All cluster are less than 1560*40=62400< 100000
//
//So CluId(PADCASE) = Takes GetDetID*1560 + unique_firstpadId_in_cluentry
//So CluId(StripCase) = 100000 + Takes GetDetID*512 + unique_firstpadId_in_cluentry



T2MakeCluster::T2MakeCluster(const edm::ParameterSet& paraSet) 
  :  theT2Clusterizer() {

  produces<T2StripClusterCollection>("T2StripClusters");
  produces<T2PadClusterCollection>("T2PadClusters");

  Projectionthreshold=paraSet.getParameter<double>("Projectionthreshold");
  BlobMinSize=paraSet.getParameter<int>("BlobMinSize");
  maskvect=paraSet.getParameter<std::vector<unsigned int> >("maskvect");
  TakeCleanEventOnly=paraSet.getParameter<bool>("TakeCleanEventOnly");
  SimuClusterEfficiency=paraSet.getParameter<bool>("SimuClusterEfficiency");
  EffiGeoRootFileName= paraSet.getParameter<std::string>("EffiGeoRootFileName");
  // verbosity=paraSet.getParameter<bool>("verbosity");

  T2StripDigiCollectionLabel = paraSet.getParameter<edm::InputTag>("T2StripDigiCollectionLabel");
  T2PadDigiCollectionLabel = paraSet.getParameter<edm::InputTag>("T2PadDigiCollectionLabel");



  edm::Service<edm::RandomNumberGenerator> rngPR;
  
  if (( ! rngPR.isAvailable())&& SimuClusterEfficiency) {
    throw cms::Exception("Configuration")
      << "This class requires the RandomNumberGeneratorService\n"
      "which is not present in the configuration file.  You must add the service\n"
      "in the configuration file or remove the modules that require it.";
  }
  else  
    rndEnginePR = &(rngPR->getEngine());


} // T2MakeCluster


void T2MakeCluster::beginJob(){
  std::cout<<"T2MakeCluster beginJob"<<std::endl;
  
  if(SimuClusterEfficiency){
    
     char *cmsswPath = getenv("CMSSW_BASE");
     if(EffiGeoRootFileName!=""){

        char sZname2[1024];

 
       EffiGeoRootFileName= string(cmsswPath) + string("/src/") + EffiGeoRootFileName;
 
       GeometryEffiFile = boost::shared_ptr<TFile> (new TFile(EffiGeoRootFileName.c_str())); 	  
       //EvtVfat_Strip_WhenCompletelyOn=(THnSparseD*)MCfileNoise->Get("VFAT_Monitoring/EvtVfat_Strip_WhenCompletelyOn");
       if ( GeometryEffiFile->IsOpen() ){
	 printf("File opened successfully\n");
	 GeometryEffiFile->ls();
	 //GeometryEffiFile->cd("EffiGeometry");

	 for(unsigned int m=0;m<40;m++){
	   sprintf(sZname2, "HGeometry_PadEfficiency %d", m); 
	   //sprintf(sZnamehist, "HGeometry_PadEfficiency %d", m);
	   //EffiPad[m]=boost::shared_ptr<TH2D>(  (static_cast<TH2D*>(GeometryEffiFile->Get(sZname2)->Clone() )) );
	   //   Data_distr_ = boost::shared_ptr<TH1>(  (static_cast<TH1*>(dataFile_->Get( DataHistName_.c_str() )->Clone() )) );
	   //MC_distr_ = boost::shared_ptr<TH1>(  (static_cast<TH1*>(generatedFile_->Get( GenHistName_.c_str() )->Clone() )) ); 
	   std::cout<<m<<" "<<sZname2<<std::endl;
	   TH2D *hname = (TH2D*)GeometryEffiFile->Get(sZname2);
	   
	   EffiPad.push_back((*hname));
	   std::cout<<m<<" OK  "<<sZname2<<std::endl;

	   sprintf(sZname2, "HGeometry_StripEfficiency %d", m); 
	   
	   TH2D *hname2 = (TH2D*)GeometryEffiFile->Get(sZname2);
	   EffiStrip.push_back((*hname2)); 
	   //EffiStrip[m]=boost::shared_ptr<TH2D>(  (static_cast<TH2D*>(GeometryEffiFile->Get(sZname2)->Clone() )) );
	   //=boost::shared_ptr<TH2D>(new TH2D(GeometryEffiFile->Get(sZname2)));
	   
	 } 
	 printf("40 Planes loaded !!!!\n");
	 EffiPad.at(0).Print();
	 printf("X bin...!\n");
	 numrbin=EffiPad.at(0).GetXaxis()->GetNbins();
	 numphibin=EffiPad.at(0).GetYaxis()->GetNbins();

	 printf("Effi work with  %d %d sectors",numrbin,numphibin);

       }else
	 printf("WARNING: File NOT OPENED !!!!\n");

     }else
        printf("WARNING: Empty Effi file Name !!!!\n");
  }
}



void T2MakeCluster::produce(edm::Event& ev, const edm::EventSetup& evSet) {

  //get the digitization information from the event

  bool eventpass=false;
  eventpass=true;
  /*
  edm::Handle<Totem::RawEvent> input;
  ev.getByType(input);
  
  if((input->timestamp==1342124437)&&(input->dataEventNumber==6342461))
    eventpass=true;

  if((input->timestamp==1342123077)&&(input->dataEventNumber==5469644))
    eventpass=true;
  */

 

  edm::Handle<T2PadDigiCollection> t2paddigicoll;
  ev.getByLabel(T2PadDigiCollectionLabel, t2paddigicoll);

  //T2PadDigiMod
  // ev.getByLabel("T2Digis","T2PadDigiMod",t2paddigicoll);

  const T2PadDigiCollection* PadDigiptr;
  PadDigiptr= t2paddigicoll.product();
  edm::Handle<T2StripDigiCollection> t2stripdigicoll;
  ev.getByLabel(T2StripDigiCollectionLabel, t2stripdigicoll);
  const T2StripDigiCollection* StripDigiptr;
  StripDigiptr= t2stripdigicoll.product(); 

  
  auto_ptr<T2StripClusterCollection> theStripClusters (new T2StripClusterCollection);
  auto_ptr<T2PadClusterCollection> thePadClusters (new T2PadClusterCollection) ;  

  if(eventpass){
  //std::auto_ptr<T2ClusterCollection> theClusters(new T2ClusterCollection());


  //theT2Clusterizer is a member of this class T2MakeCluster to be declared in .h file;  its type is T2DetClusterReconst
  //Take all the digipad and digistrip from all the detector for this event , and make clusterization inside *theStripClusters and *thePadClusters
  
  // std::cout<<"Event id:"<<ev.id().event()<<std::endl;
  unsigned int symb;
  std::vector<int> myvp;
  std::vector<int> myvs;

  bool detmasked;
  DigiContainerIterator<T2DetId, T2PadDigi> itp;   
  DigiContainerIterator<T2DetId, T2StripDigi> its;

  std::vector<T2Cluster> myclusters; 
  
 bool proceedWithReco = true;
 int numoffullstripplane=0;
 int numoffullpadplane=0;
 //Added on 15 June 2010
 if(TakeCleanEventOnly)
  {
    for(itp= PadDigiptr->begin(); itp!=PadDigiptr->end(); ++itp)
    {

      T2DetId mydet=(*itp).first;
      symb= mydet.plane()*2 + mydet.planeSide() + mydet.arm()*20 + mydet.halfTelescope()*10;
   
      int numofpadsintheplane=(*itp).second.second - (*itp).second.first ;
      if(numofpadsintheplane > 50)
	numoffullpadplane++;
     }

    for(its= StripDigiptr->begin(); its!=StripDigiptr->end(); ++its)
    {

      T2DetId mydet=(*itp).first;
      symb= mydet.plane()*2 + mydet.planeSide() + mydet.arm()*20 + mydet.halfTelescope()*10;
      int numofstripsintheplane=(*its).second.second - (*its).second.first ;
      if(numofstripsintheplane > 50)
	numoffullstripplane++;
     }
    
    if((numoffullstripplane>7)&&(numoffullpadplane>7))
	proceedWithReco=false;

    std::cout<<"Number of plane with more the 50 pad/Strip ON: "<<numoffullpadplane<<"/"<<numoffullstripplane<<std::endl;
  }


  if(proceedWithReco)	
  for(itp= PadDigiptr->begin(); itp!=PadDigiptr->end(); ++itp)
    {
          
      T2DetId mydet=(*itp).first;  
      
      //symb=mydet.arm()*20+mydet.plane()*4+mydet.halfTelescope()*2+mydet.planeSide();      
      symb= mydet.plane()*2 + mydet.planeSide() + mydet.arm()*20 + mydet.halfTelescope()*10;

      detmasked=false;
      for(unsigned int j=0;j<maskvect.size();j++)
	{
	  if(symb==maskvect.at(j))
	     detmasked=true;
	}


	for(std::vector<T2PadDigi>::const_iterator itpad =(*itp).second.first; itpad !=(*itp).second.second; itpad++)
	{
	  //mypadsdigi=(*itpad).second; 

	  if((*itpad).getPadNr()==0)//This discriminate data vs simu
	    myvp.push_back(((*itpad).getRow()+1)+24*((*itpad).getCol()) );   //Prima del dataraw
	  else
	    myvp.push_back(((*itpad).getPadNr())); 

	}
     
	theT2Clusterizer.Projectionthreshold=Projectionthreshold;
	theT2Clusterizer.BlobMinSize=BlobMinSize;
	theT2Clusterizer.SetDetId(mydet.rawId());

	theT2Clusterizer.SetPadHits(myvp);

	theT2Clusterizer.BuildClusters();

	myvp.clear();
	myclusters=theT2Clusterizer.GetPadClusters();
      

	////////////////////////////////////////////////////////////////////////
	//////////Simulate Cluster EffiBeg
	std::vector<T2Cluster> myclustersWeightedInEffi;
	if(SimuClusterEfficiency){
	  for(unsigned int gg=0;gg<myclusters.size();gg++){
	    double r=myclusters.at(gg).GetClusterR();
	    double phi=myclusters.at(gg).GetClusterPhi();
	    uint32_t i=myclusters.at(gg).GetDetID();
	    
	    int xx=0;
	    int yy=0;
	    int symbId=0;
	    GetEffiSector(r,phi,i,xx,yy,symbId,numrbin,numphibin);
	    //  unsigned int rseed=0;
	    //rseed=(unsigned int)( CLHEP::RandFlat::shoot(rndEngineStr) * std::numeric_limits<unsigned int>::max() );
	    //gRandom->SetSeed(rseed);
	    double randomValue=(double)(CLHEP::RandFlat::shoot(0.,1.));
	    double padeffi=EffiPad[symbId].GetBinContent(xx+1,yy+1)/1000.; 
	    //Now I force a 2% larger efficiency since in the simulation 2% is eated from digi
	    if(padeffi>0.05){//Dead VFAT are DEAD VFAT!
	      if((padeffi+0.02)<1.)
		padeffi=padeffi+0.02;
	      else
		padeffi=1.;
	    }

	    
	
	    
	    if(padeffi>randomValue){
	      myclustersWeightedInEffi.push_back(myclusters.at(gg));
//	      std::cout<<"Plane: "<<symbId<<"Rand: "<<randomValue<<" Pad Effi Read:"<<streffi<<"Sector:"<<xx<<" "<<yy<<std::endl;
	    }
	  }
	  
	  myclusters.clear();
	  myclusters=myclustersWeightedInEffi;
	}
      //////////Simulate Cluster EffiEnd
      /////////////////////////////////////////////////////////////////////////


	if(detmasked==false)
	  thePadClusters->insert(std::make_pair(mydet,myclusters));     
	
	myclusters.clear();
    }
  
 
  

  if(proceedWithReco)
  for(its= StripDigiptr->begin(); its!=StripDigiptr->end(); ++its)

    {
      T2DetId mydet=(*its).first;
      //symb=mydet.arm()*20+mydet.plane()*4+mydet.halfTelescope()*2+mydet.planeSide();
      symb= mydet.plane()*2 + mydet.planeSide() + mydet.arm()*20 + mydet.halfTelescope()*10;
      detmasked=false;
      for(unsigned int j=0;j<maskvect.size();j++)
	{
	  if(symb==maskvect.at(j))
	     detmasked=true;
	}

     
      for(std::vector<T2StripDigi>::const_iterator itstrip =(*its).second.first; itstrip !=(*its).second.second; itstrip++)
	{
	
	  if((*itstrip).getStripNr()==0)
	    myvs.push_back(((*itstrip).getRow()+1)+256*((*itstrip).getCol()) ); //PrimaDataraw
	  else
	    myvs.push_back(((*itstrip).getStripNr())); 
	       
	  //     cout<<" Detector pl: "<< mydet.plane() <<". Strip digitizzati di Erik: Row:  "<< (*itstrip).getRow() <<" , Col: "<< (*itstrip).getCol() <<endl;
    
	}
 
      theT2Clusterizer.Projectionthreshold=Projectionthreshold;
      theT2Clusterizer.BlobMinSize=BlobMinSize;
      theT2Clusterizer.SetDetId(mydet.rawId());
      theT2Clusterizer.SetStripHits(myvs); 
      theT2Clusterizer.BuildClusters();

      //Part added for noise studies      
      //StrCluStrRawID_=theT2Clusterizer.StrCluStrRawID;
      //StrCluStrColID_=theT2Clusterizer.StrCluStrColID;
      //PadCluPadRawID_=theT2Clusterizer.PadCluPadRawID;
      //PadCluPadColID_=theT2Clusterizer.PadCluPadColID;


      myvs.clear();
     
      myclusters=theT2Clusterizer.GetStripClusters();

      ////////////////////////////////////////////////////////////////////////
      //////////Simulate Cluster EffiBeg
      std::vector<T2Cluster> myclustersWeightedInEffi;
      if(SimuClusterEfficiency){
	for(unsigned int gg=0;gg<myclusters.size();gg++){
	  double r=myclusters.at(gg).GetClusterR();
	  double phi=myclusters.at(gg).GetClusterPhi();
	  uint32_t i=myclusters.at(gg).GetDetID();

	  int xx=0;
	  int yy=0;
	  int symbId=0;
	  GetEffiSector(r,phi,i,xx,yy,symbId,numrbin,numphibin);
	  //  unsigned int rseed=0;
	  //rseed=(unsigned int)( CLHEP::RandFlat::shoot(rndEngineStr) * std::numeric_limits<unsigned int>::max() );
	  //gRandom->SetSeed(rseed);
	  double randomValue=(double)(CLHEP::RandFlat::shoot(0.,1.));
	  double streffi=EffiStrip[symbId].GetBinContent(xx+1,yy+1)/1000.; 
	  //Now I force a 5% larger efficiency since in the simulation 5% is eated from digi
	 
	  if(streffi>0.05){//Dead VFAT are DEAD VFAT!
	    if((streffi+0.06)<1.)
	      streffi=streffi+0.06;
	    else
	      streffi=1.;
	  }
	  
	  
	 

	  if(streffi>randomValue){
	    myclustersWeightedInEffi.push_back(myclusters.at(gg));
	    //if(symbId==35)
	    //  std::cout<<"Pl35 strip:"<<streffi<<" Rand:"<<randomValue<<"saved"<<std::endl;
	    
	  }
	}
      
	myclusters.clear();
	myclusters=myclustersWeightedInEffi;
      }
      //////////Simulate Cluster EffiEnd
      /////////////////////////////////////////////////////////////////////////


      if(detmasked==false)
	theStripClusters->insert(std::make_pair(mydet,myclusters));

    }
  
  }//if(eventpass)

  ev.put(theStripClusters, "T2StripClusters");
  ev.put(thePadClusters, "T2PadClusters");

} // produce




void T2MakeCluster::GetEffiSector(double Rcl, double Phicl, uint32_t thedet, int &rCell, int &phiCell, int &symbId, int numRsectEffi, int numPhisectEffi){
  
  //  144 and 41 mm correspond to 5.3 and 6.5
  int rsect=0;int phisect=0;
  int plane0_39=0;
 
  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo=conv.GetT2Info(thedet);
  plane0_39=planeinfo.symb;
  symbId=planeinfo.symb;
  double zpiano=14000.;
  zpiano=planeinfo.Zdet;



  double etacluster=-log((Rcl/fabs(zpiano))/2.0);
  
  rsect=(int)(((6.5-etacluster)/(6.5-5.3))*numRsectEffi);
 

  if(Rcl<41.)
    rsect=0;
  if((Rcl>144.)||(rsect==numRsectEffi))
    rsect=numRsectEffi-1;
  /*
  if((rsect<0)||(rsect>numRsectEffi))
    std::cout<<"Error in sector of ConvertCluRPhiInEffiGeoCell "<<rsect<<std::endl;
  */


  bool isnear=false;
  bool isfar=false;
  
  if(plane0_39>=30)
    isfar=true;

  if(plane0_39>=10)
    if(plane0_39<20)
      isfar=true;

  if(isfar==false)
    isnear=true;
 
  //clockwise increasing sector

  if(isnear){

    bool warn=false;
    
    if(Phicl<264.)
      if(Phicl>180.){
	phisect=0;
	warn=true;
      }

    if(Phicl>96.)
      if(Phicl<110.){
      phisect=numPhisectEffi-1;
      warn=true;
    }

    if(warn==false){
      
      if(Phicl<250)//isUp
	Phicl=Phicl+360.;
      
      Phicl=Phicl-264.; //Now PhiCl goes from 0 to 192.
      phisect=(int)((Phicl/192.)*numPhisectEffi);
    }
    //std::cout<<"RSec: "<<rsect<<" -> Eta Cluster "<<etacluster<<" PhiSec: "<<phisect<<" Phicl: "<<Phiclbeg<<" -> "<<Phicl<<std::endl;
    

  }

  if(isfar){

    bool warn=false;
    
    if(Phicl<84.){
	phisect=0;
	warn=true;
      }

    if(Phicl>276.){
      phisect=numPhisectEffi-1;
      warn=true;
    }

    if(warn==false){
      
      
      Phicl=Phicl-84.; //Now PhiCl goes from 0 to 192.
      phisect=(int)((Phicl/192.)*numPhisectEffi);
    }

  }
  
  if(phisect==numPhisectEffi)
    phisect=numPhisectEffi-1;


 
  
  rCell=rsect;
  phiCell=phisect;
}

/**
 * Class T2TrackProducer3
 *
 * Author: Mirko Berretti / University of Siena 
 * Email:  mirko.berretti@gmail.com
 * Date:   2007-12-08
 * Revision: 2008-10-10 
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"



#include "RecoTotemT1T2/T2TrackProducer3/interface/T2TrackProducer3.h"
#include "DataFormats/T1T2Track/interface/T1T2Track.h"
#include "DataFormats/T1T2Track/interface/T1T2TrackCollection.h"
#include "DataFormats/T2Road/interface/T2Road.h"
#include "DataFormats/T2Road/interface/T2RoadCollection.h"
#include "DataFormats/T2Hit/interface/T2Hit.h"
#include "DataFormats/T2Hit/interface/T2HitCollection.h"
#include "DataFormats/T2Hit/interface/T2Hit_to_Track_Map.h"
#include <iostream>
#include <cmath>



T2TrackProducer3::T2TrackProducer3(const edm::ParameterSet &config): theChiRidThr(1.5){
  
  produces<T1T2TrackCollection>("T2TrackColl");  
  theChiRidThr = config.getParameter<double>("ChiRidThr");    
  RoadModuleLabel = config.getParameter<std::string>("RoadModuleLabel");
  RoadInstanceLabel = config.getParameter<std::string>("RoadInstanceLabel");

  theT2Hit_to_Track_Map0_Module= config.getParameter<std::string>("theT2Hit_to_Track_Map0_Module");
  theT2Hit_to_Track_Map0_Instance= config.getParameter<std::string>("theT2Hit_to_Track_Map0_Instance");

  MinHitFinalTrk = config.getParameter<unsigned int>("MinHitFinalTrk");
  MinHitFinalTrkCl1Hit= config.getParameter<unsigned int>("MinHitFinalTrkCl1Hit");
  MaxHitsInRoad= config.getParameter<unsigned int>("MaxHitsInRoad");
  verbosity=config.getParameter<bool>("verbosity");
  forceXYfit = config.getParameter<bool>("forceXYfit");    
  forceRZfit = config.getParameter<bool>("forceRZfit");    
  TrackInstanceLabel= config.getParameter<std::string>("TrackInstanceLabel");
  
  UseRefittingProcedure = config.getParameter<bool>("UseRefittingProcedure");	  
  RecoHitRErrorPerStripcount = config.getParameter<std::vector<double>  >("RecoHitRErrorPerStripcount");
  DropWorstRecHitChisquareRThreshold = config.getParameter<double>("DropWorstRecHitChisquareRThreshold");

  //Ideated for Alignment analysis bias
  PickUpDisplacedHit= config.getParameter<bool>("PickUpDisplacedHit");
  PickUpRadius= config.getParameter<double>("PickUpRadius");
  
  
  FitVertexPosition=config.getParameter<bool>("FitVertexPosition");
  
  VtxPositionX=config.getParameter<double>("VtxPositionX");
  VtxPositionY=config.getParameter<double>("VtxPositionY");
  VtxPositionZ=config.getParameter<double>("VtxPositionZ");
  
  VtxPositionEX=config.getParameter<double>("VtxPositionEX");
  VtxPositionEY=config.getParameter<double>("VtxPositionEY");
  VtxPositionEZ=config.getParameter<double>("VtxPositionEZ");
  StripFitting = config.getParameter<bool>("StripFitting"); 
  HitLabel= config.getParameter<string>("HitLabel");

  RemoveOutliers = config.getParameter<bool>("RemoveOutliers");
  GhostSuppression = config.getParameter<bool>("GhostSuppression");
}

T2TrackProducer3::~T2TrackProducer3(){

}



void T2TrackProducer3::produce(edm::Event& event, const edm::EventSetup& setup) {

  //std::cout<<"Producer Track2"<<std::endl;
 
  edm::Handle<T2RoadCollection> myRoadColl;
  //event.getByType(myRoadColl);
  //  event.getByLabel(RoadModuleLabel,"T2RoadColl",myRoadColl);
  event.getByLabel(RoadModuleLabel,RoadInstanceLabel,myRoadColl);
  

  edm::Handle<T2HitCollection> t2hitcoll;
  event.getByLabel(HitLabel,"T2Hits",t2hitcoll);
  const T2HitCollection *T2HitCollectionPtr = t2hitcoll.product();

  edm::Handle<T2Hit_to_Track_Map> theT2Hit_to_Track_Map0;
  event.getByLabel(theT2Hit_to_Track_Map0_Module,theT2Hit_to_Track_Map0_Instance,theT2Hit_to_Track_Map0);

  const T2Hit_to_Track_Map *theT2Hit_to_Track_Map = theT2Hit_to_Track_Map0.product();

  auto_ptr<T1T2TrackCollection> theT2Tracks (new T1T2TrackCollection());  //Container for all the event tracks
  auto_ptr<T1T2TrackCollection> theT2Tracks_NoGhost (new T1T2TrackCollection());  //Container for all the event tracks
  
  
  T2GeometryUtil convGeo;
  T2GeometryUtil::T2DetInfo planeinfo;


  if(MinHitFinalTrk<3)
    {
      MinHitFinalTrk=3; //3 Hit with at leat 2 class1
      std::cout<<"WARNING: T2TrackProducer3 minimum cl1 hit in track forced to 3"<<std::endl;
    }
  //    MinHitFinalTrk=3;// before was 3
 
  /*
  if(verbosity)
       {
	 std::cout<<"|||||||-------------------------------------------------------------"<<std::endl;
	 std::cout<<"-------------------------------------------------NEW EVENT |||||||||"<<std::endl;
	 std::cout<<"|||||||-------------------------------------------------------------"<<std::endl;
       }
  */

  std::vector<double> trkparam;
  //std::vector<double> trkparam2;
  std::vector<double> unused;
  unused.clear();




 for(T2RoadCollection::const_iterator itroad = myRoadColl->begin();
       itroad != myRoadColl->end(); itroad++)
   {
     T2Road theroad=(*itroad);
    
     
     if((theroad.thisRoad.size()>=3))//&&(theroad.thisRoad.size()>2) //Avoid to loop on pad-only pads
       {
	 //std::cout<<"---------------------------------------------"<<std::endl;
	 //std::cout<<"TrkProd3 RoadId: "<<theroad.RoadID<<std::endl;
	 
	 unsigned int totblob=0;
	 unsigned int totnumpadsmall=0;
	 int roadId=0;
	 for(unsigned int jj=0;jj<theroad.thisRoad.size();jj++)
	   {
	     if(theroad.thisRoad.at(jj).GetHitNumPad()>6){
	       totblob++;
	     }else
	       {
		 if(theroad.thisRoad.at(jj).GetHitNumPad()>0)
		   totnumpadsmall++;
	       }
	   }

	 if(totnumpadsmall>totblob){
	   //std::cout<<"TrkProd3 RoadId: "<<theroad.RoadID<<std::endl;
	 trkparam.clear();
	 //trkparam2.clear();
	 TMatrixD covmat(4,4);
	 covmat.Zero();
	 double chi2corr=0.;		
	 //T1T2Track firstTrk=MyLinearfitCorrDEV(theroad.thisRoad,covmat,chi2corr,false,unused,theroad.RoadID);
	 T2Hit worsthit;
	 roadId=theroad.RoadID;

	 T1T2Track firstTrkRaw=MyLinearfitCorrDEV(theroad.thisRoad,covmat,chi2corr,false,unused,roadId,worsthit);

	 double ProbChi2_XY=TMath::Prob(firstTrkRaw.ChiSquared(),(firstTrkRaw._numCl1HitHit_t2*2-4));

	 if(verbosity)
	 if(ProbChi2_XY<0.01){
	   std::cout<<"Warning! First Lowchiprob:"<<ProbChi2_XY<<" Chi2:"<<firstTrkRaw.ChiSquared()<<" NumCl1Hit:"<<firstTrkRaw._numCl1HitHit_t2<<std::endl;
	   for(unsigned int i=0; i<firstTrkRaw.GetHitEntries();i++)
	     {
	       // if((*TrkCit).GetHitT2(i).GetHitClass()==1){
	       std::cout<<"X-Y-Z: "<<firstTrkRaw.GetHitT2(i).GetHitX()<<" "<<firstTrkRaw.GetHitT2(i).GetHitY()<<" "<<firstTrkRaw.GetHitT2(i).GetHitZ()<<" class:"<<firstTrkRaw.GetHitT2(i).GetHitClass()<<" NumPad:"<<firstTrkRaw.GetHitT2(i).GetHitNumPad()<<" NumStrip:"<<firstTrkRaw.GetHitT2(i).GetHitNumStrip()<<" DX-DY:"<<firstTrkRaw.GetHitT2(i).GetHitDX()<<" "<<firstTrkRaw.GetHitT2(i).GetHitDY()<<" dR:"<<firstTrkRaw.GetHitT2(i).GetHitDR()<<std::endl;            
	     }
	 }


	 T1T2Track firstTrk;
	 




	 if(PickUpDisplacedHit){
	   
	   if(verbosity)
	     if(RemoveOutliers)
	       std::cout<<"Warning! remove outliers incompatible with PickUpDisplacedHit options"<<std::endl;
	   
	   std::vector<T2Hit> InitTrkhit;std::vector<uint32_t> trkhitplanes;
	   for(unsigned int i=0; i<firstTrkRaw.GetHitEntries();i++)
	     {	    
	      InitTrkhit.push_back(firstTrkRaw.GetHitT2(i));            
	      trkhitplanes.push_back(firstTrkRaw.GetHitT2(i).GetHitDetRawId());	     
	     }

	   std::vector<T2Hit> HitToAdd;
	   std::vector<uint32_t> IdPlaneADD;
	   uint32_t rawiddet=0;uint32_t rawiddetB=0;
	   unsigned int  symbh=0; unsigned int  symbt=0;
	   if(verbosity)
	     std::cout<<"InitTrk size: "<<InitTrkhit.size()<<std::endl;   
	   double projX=0.;double projY=0.;double thisDR=0.;
           rawiddetB=trkhitplanes.at(0);
	   planeinfo=convGeo.GetT2Info(rawiddetB);
	   symbt=planeinfo.symb;


	   for(T2HitCollection::const_iterator ithit = T2HitCollectionPtr->begin();ithit!=T2HitCollectionPtr->end();ithit++){
	     rawiddet=(*ithit).GetHitDetRawId();
	     planeinfo=convGeo.GetT2Info(rawiddet);//RawtoSymb(rawiddet);
	     symbh=planeinfo.symb;
	     
	     bool candidateToinclude=false;	     
	     if(std::find(trkhitplanes.begin(),trkhitplanes.end(),rawiddet)== trkhitplanes.end())    
	       if((symbt/10)==(symbh/10))
		 candidateToinclude=true;




	     if(candidateToinclude)
	     if(((*ithit).GetHitNumStrip()>0)&&((*ithit).GetHitNumPad()>0))
	       {
		 //std::cout<<"Candidate StripR"<<<<" Strip-Plane:"<<stripplane<<std::endl;
		 projX=firstTrkRaw.X0()+firstTrkRaw.GetTx()*(*ithit).GetHitZ();
		 projY=firstTrkRaw.Y0()+firstTrkRaw.GetTy()*(*ithit).GetHitZ();

		 thisDR=(projX-(*ithit).GetHitX())*(projX-(*ithit).GetHitX());
		 thisDR=thisDR+(projY-(*ithit).GetHitY())*(projY-(*ithit).GetHitY());
		 thisDR=sqrt(thisDR);
		 // thisDR=fabs(sqrt(projX*projX+projY*projY)-(*ithit).GetHitR());
		 
		 if(thisDR<PickUpRadius)
		   {
		     if(std::find(IdPlaneADD.begin(),IdPlaneADD.end(),rawiddet)== IdPlaneADD.end())
		       {

			 
			 long int posincoll= (long int) (*ithit).Hit_PosInCollection;
			 //std::map<long int,int>::iterator itera = theT2Hit_to_Track_Map->find(posincoll);
			 
			 //check if this position is registered in theT2Hit_to_Track_Map where only hits associated to
			 //some road are saved. Proceed with the association also if it is not already associated to some
			 //other tracks 

			 if(theT2Hit_to_Track_Map->find(posincoll)==theT2Hit_to_Track_Map->end()){
			   HitToAdd.push_back((*ithit));
			   IdPlaneADD.push_back(rawiddet);
			   trkhitplanes.push_back(rawiddet);
			 }
		       }
		     
		   }
		 
	       }
	     
	   }
	   for(unsigned int hh=0;hh<HitToAdd.size();hh++){
	     InitTrkhit.push_back(HitToAdd.at(hh));   
	     if(verbosity)
	       std::cout<<"Pick-up something in plane"<<(HitToAdd.at(hh).GetHitPlane()*2+HitToAdd.at(hh).GetHitPlaneSide())<<" "<<HitToAdd.at(hh).GetHitHalftele()<<" "<<HitToAdd.at(hh).GetHitArm()<<std::endl;
	   }
	   T1T2Track firstTrkPickUP=MyLinearfitCorrDEV(InitTrkhit,covmat,chi2corr,false,unused,roadId,worsthit);
	   
	   firstTrk=firstTrkPickUP;
	   firstTrkRaw=firstTrkPickUP;

	 }










	 if((ProbChi2_XY<0.01)&&(RemoveOutliers)&&(totnumpadsmall>3)){
	 
	   bool enoughhit=true; 
	   bool trackcleaned=false;	   
	   
	   std::vector<T2Hit> hitcleaned;
	   for(unsigned int i=0; i<firstTrkRaw.GetHitEntries();i++)
	     {	    
	       hitcleaned.push_back(firstTrkRaw.GetHitT2(i));            
	     }
	   
	   while((enoughhit==true)&&(trackcleaned==false)){
	     
	     T2Hit worsthit2;
	     T1T2Track currenttrack=MyLinearfitCorrDEV(hitcleaned,covmat,chi2corr,false,unused,roadId,worsthit2);
	     	    
	     //Remove an hit
	     std::vector<T2Hit> actualhitcleaned=VectOutliersRemoving(hitcleaned,worsthit2);
	      
	     if(verbosity)
	       std::cout<<"Removed Outlier X-Y:  "<<worsthit2.GetHitX()<<" "<<worsthit2.GetHitY()<<std::endl;
	     //Recalculate
	     T1T2Track Aftercurrenttrack=MyLinearfitCorrDEV(actualhitcleaned,covmat,chi2corr,false,unused,roadId,worsthit2);
	     double ProbChi2=TMath::Prob(Aftercurrenttrack.ChiSquared(),(Aftercurrenttrack._numCl1HitHit_t2*2-4));	     
	     //std::cout<<"Remove outliers: Hit survived: "<<Aftercurrenttrack._numCl1HitHit_t2<<". Prob:"<<ProbChi2<<std::endl;
	    
	     if(ProbChi2>0.01){
	       trackcleaned=true;
	       firstTrk=Aftercurrenttrack;
	       //std::cout<<"..Track recovered!! Numhit="<<firstTrk.GetHitEntries()<<std::endl;
	       //std::cout<<"..Track Recovered !!"<<std::endl;
	     }
	     
	     if(Aftercurrenttrack._numCl1HitHit_t2<=4){
	       enoughhit=false;
	       trackcleaned=false;
	     }
	     
	     hitcleaned.clear();
	     hitcleaned=actualhitcleaned;
	    	            
	   }
	   //Recopy in the final trk, the initial trk
	   if(trackcleaned==false){
	     firstTrk=firstTrkRaw;
	    
	     //std::cout<<"..Unable to remove enough outliers. Initial Trk saved with numhit="<<firstTrk.GetHitEntries()<<std::endl;
 	   }

	 }else{
	   
	   firstTrk=firstTrkRaw;
	   
	   // 
	 }
	 
	
	 
	 ProbChi2_XY=TMath::Prob(firstTrk.ChiSquared(),(firstTrk._numCl1HitHit_t2*2-4));
	 
	 if(verbosity)
	 if(ProbChi2_XY<0.01){
	   std::cout<<"Warning! Final Chiprob<0.01:"<<ProbChi2_XY<<" Chi2:"<<firstTrk.ChiSquared()<<" Chi2X:"<<firstTrk.ChiSquaredX()<<" Chi2Y:"<<firstTrk.ChiSquaredY()<<" NumCl1Hit:"<<firstTrk._numCl1HitHit_t2<<std::endl;
	   for(unsigned int i=0; i<firstTrk.GetHitEntries();i++)
	     {
	       
	       std::cout<<"X-Y-Z: "<<firstTrk.GetHitT2(i).GetHitX()<<" "<<firstTrk.GetHitT2(i).GetHitY()<<" "<<firstTrk.GetHitT2(i).GetHitZ()<<" class:"<<firstTrk.GetHitT2(i).GetHitClass()<<" NumPad:"<<firstTrk.GetHitT2(i).GetHitNumPad()<<" NumStrip:"<<firstTrk.GetHitT2(i).GetHitNumStrip()<<"| ExpX: "<<firstTrk.GetTx()*firstTrk.GetHitT2(i).GetHitZ()+firstTrk.X0()<<" ExpY: "<<firstTrk.GetTy()*firstTrk.GetHitT2(i).GetHitZ()+firstTrk.Y0()<<std::endl;            
	     }
	 }

	 
	 

	 //Removing outliers (only if track remain with at least 3cl1hits. Otherwise firstTrk=firstTrkRaw)  
	 //T1T2Track firstTrk=OutliersRemovingDEV(firstTrkRaw);
	 //OutliersRemovingDEV(firstTrkRaw);
	 
	 if(StripFitting)
	   {
	     trkparam.push_back(firstTrk.GetTx());trkparam.push_back(firstTrk.X0());	 
	     trkparam.push_back(firstTrk.GetTy());trkparam.push_back(firstTrk.Y0());	 
	     T2Hit worstHit;
	     // T1T2Track finalTrk=MyLinearfitCorrDEV(theroad.thisRoad,covmat,chi2corr,true,trkparam,theroad.RoadID);
	     T1T2Track finalTrk=MyLinearfitCorrDEV(theroad.thisRoad,covmat,chi2corr,true,trkparam,roadId,worstHit);
	     theT2Tracks->push_back(finalTrk);
	   }
	 else
	   {
	     theT2Tracks->push_back(firstTrk);
	   }



	 }
	 
       }
     else
       std::cout<<"Warning in T2TrackProducer3: found Road with less than 3 hits"<<std::endl;
     
   } 
 //end of event Road collection 
 
  if(GhostSuppression){
   //each row contain the tracks belonging to the same pad road. the key is the roadId moduled 10000
   //Note from roadProducer without concurrent branches: theRoadwithStrip.RoadID =  10000*t + itroad->RoadID; 
   std::map<unsigned int, std::vector<T1T2Track> > trackMatrixWithGhost; 
   std::vector<T1T2Track> trackVectorNoGhost; 
   
   for(T1T2TrackCollection::const_iterator recTrack = theT2Tracks->begin(); recTrack!=theT2Tracks->end();++recTrack){
     //unsigned int t2_roadID
     unsigned int trk_roadID=(*recTrack).t2_roadID;
     
     //rescaled trk_roadID to the common pad
     trk_roadID=trk_roadID%10000;
     std::map<unsigned int, std::vector<T1T2Track> >::iterator Matr_it;
     Matr_it=trackMatrixWithGhost.find(trk_roadID);

     if(Matr_it==trackMatrixWithGhost.end()){
       
       std::vector<T1T2Track> OneElementVectr;
       OneElementVectr.push_back((*recTrack));
       trackMatrixWithGhost.insert ( pair<unsigned int, std::vector<T1T2Track> >(trk_roadID,OneElementVectr) );             
     }else{
       //Add a proably-ghost Track
       trackMatrixWithGhost[trk_roadID].push_back((*recTrack));       
     }
         
   }
   

   //Now for each row of the matrix decide the better track. The choice is driven by 
   //Chi2-Prob and number of hits. 
   //Strategy:
   //All Chi2<0.1-> choose the one with more hits. Otherwise choose the best chi2.
   
   std::vector<double> chi2V;
   std::vector<unsigned int> multiplicityV;
   // std::vector<unsigned int> index_winnerTrack;

   for(std::map<unsigned int, std::vector<T1T2Track> >::iterator Map_iter=trackMatrixWithGhost.begin();Map_iter!=trackMatrixWithGhost.end();++Map_iter){
     //Take a row;
     std::vector<T1T2Track> matrixRow=Map_iter->second;
     
     if(matrixRow.size()==1)
       theT2Tracks_NoGhost->push_back(matrixRow.at(0));
     else{
       chi2V.clear();
       multiplicityV.clear();
       unsigned int indexBestTrk=0; 
   
       double MaxProbChi2_XY=0.;   
       unsigned int MaxProbChi2_XY_index=0;   
       
       unsigned int MaxMultipl=0.;   
       unsigned int MaxMultipl_index=0; 

       for(unsigned int j=0;j<matrixRow.size();j++){

	 unsigned int thisCl1mult=matrixRow.at(j)._numCl1HitHit_t2;
	 double ProbChi2_XY=TMath::Prob(matrixRow.at(j).ChiSquared(),(thisCl1mult*2-4));
	 
	 chi2V.push_back(ProbChi2_XY);
	 multiplicityV.push_back(thisCl1mult);
	 
	 if(ProbChi2_XY>MaxProbChi2_XY){
	   MaxProbChi2_XY=ProbChi2_XY;
	   MaxProbChi2_XY_index=j;
	 }
	 if(thisCl1mult>MaxMultipl){
	   MaxMultipl=thisCl1mult;
	   MaxMultipl_index=j;
	 }	 
       }
       
       if(MaxProbChi2_XY>0.1)
	 indexBestTrk=MaxProbChi2_XY_index;
       else{
	 indexBestTrk=MaxMultipl_index;
       }
	 

       theT2Tracks_NoGhost->push_back(matrixRow.at(indexBestTrk));
     }
     
   }
   
 }//end if(GhostSuppression)



 if(GhostSuppression)
   event.put(theT2Tracks_NoGhost, "T2TrackColl");
 else
   event.put(theT2Tracks, "T2TrackColl");

 //event.put(theT2Tracks, "T2TrackColl");

// std::cout<<"found "<<theT2Tracks->size()<<" tracks in event"<<std::endl;



} // produce


std::vector<T2Hit> T2TrackProducer3::VectOutliersRemoving(std::vector<T2Hit> &hitv, T2Hit &worsthit){

  std::vector<T2Hit> hitcleaned;


  unsigned int largerHitindex=0; int largersize=0;

  for(unsigned int i=0; i<hitv.size();i++)
    {
      if((hitv.at(i).GetHitNumPad()>6)||(hitv.at(i).GetHitNumStrip()>6)){
	double size=(hitv.at(i).GetHitNumPad())*(hitv.at(i).GetHitNumStrip());
	if(size>largersize){
	  largersize=size;
	  largerHitindex=i;
	}
      }

    }

  // std::cout<<"largersize Hit:"<<largersize<<std::endl;

  if(largersize>6){
    for(unsigned int i=0; i<hitv.size();i++)
      if(i!=largerHitindex){
	 hitcleaned.push_back(hitv.at(i));
      }
      else
	worsthit=hitv.at(i);
  }
  else{
    hitcleaned.clear();
    for(unsigned int i=0; i<hitv.size();i++)
      {
	if(hitv.at(i).HitUniqueId!=worsthit.HitUniqueId)
	  hitcleaned.push_back(hitv.at(i));
	
	// if((*TrkCit).GetHitT2(i).GetHitClass()==1){
	//    std::cout<<"X-Y-Z: "<<firstTrkRaw.GetHitT2(i).GetHitX()<<" "<<firstTrkRaw.GetHitT2(i).GetHitY()<<" "<<firstTrkRaw.GetHitT2(i).GetHitZ()<<" class:"<<firstTrkRaw.GetHitT2(i).GetHitClass()<<" NumPad:"<<firstTrkRaw.GetHitT2(i).GetHitNumPad()<<" NumStrip:"<<firstTrkRaw.GetHitT2(i).GetHitNumStrip()<<" DX-DY:"<<firstTrkRaw.GetHitT2(i).GetHitDX()<<" "<<firstTrkRaw.GetHitT2(i).GetHitDY()<<std::endl;
	
      }
  }

  return hitcleaned;
  
}

std::vector<T2Hit>  T2TrackProducer3::OutliersRemoving(T1T2Track &firstTrkRaw,T2Hit &worsthit){

  std::vector<T2Hit> hitcleaned;
  for(unsigned int i=0; i<firstTrkRaw.GetHitEntries();i++)
    {
      if(firstTrkRaw.GetHitT2(i).HitUniqueId!=worsthit.HitUniqueId)
	hitcleaned.push_back(firstTrkRaw.GetHitT2(i));
      // if((*TrkCit).GetHitT2(i).GetHitClass()==1){
	//    std::cout<<"X-Y-Z: "<<firstTrkRaw.GetHitT2(i).GetHitX()<<" "<<firstTrkRaw.GetHitT2(i).GetHitY()<<" "<<firstTrkRaw.GetHitT2(i).GetHitZ()<<" class:"<<firstTrkRaw.GetHitT2(i).GetHitClass()<<" NumPad:"<<firstTrkRaw.GetHitT2(i).GetHitNumPad()<<" NumStrip:"<<firstTrkRaw.GetHitT2(i).GetHitNumStrip()<<" DX-DY:"<<firstTrkRaw.GetHitT2(i).GetHitDX()<<" "<<firstTrkRaw.GetHitT2(i).GetHitDY()<<std::endl;
      
    }
  return hitcleaned;
}


 T1T2Track T2TrackProducer3::MyLinearfitCorr(std::vector<T2Hit> hitvec,TMatrixD &par_covariance_matrix,double &chi2_,bool Usestrip,std::vector<double> trkparam, int RoadID)
{
   

   unsigned int sizeHitv=hitvec.size();   
  

   int hemisphere=0; 

   if(sizeHitv>=2)
     {
       hemisphere = (int)(hitvec[1].GetHitZ()/fabs(hitvec[1].GetHitZ()));
     }
   else
     {
       if(sizeHitv>0)
        hemisphere = (int)(hitvec[0].GetHitZ()/fabs(hitvec[0].GetHitZ()));
	std::cout<<"ERROR in T2TrackProducer3::fitTrackspecialXY: using vtx position to determine the hemisphere "<<std::endl;
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
   double a_xz=0; double b_xz=0; double a_yz=0; double b_yz=0;
   x.clear(); y.clear(); z.clear(); ex.clear(); ey.clear(); ez.clear(); 
 

   std::vector<T2Hit> hitvec2;hitvec2.clear();


   unsigned int numHit_t2=0;
   unsigned int numCl1HitHit_t2=0;
   unsigned int numStripOnlyHit_t2=0;
   unsigned int numPadOnlyHit_t2=0;

   for(unsigned int u=0;u<sizeHitv;u++)
     {


       numHit_t2++;
       
       if(hitvec[u].GetHitNumPad()>0){

	 if(hitvec[u].GetHitNumPad()<5)
	   hitvec2.push_back(hitvec[u]);

	 if(hitvec[u].GetHitNumStrip()==0)
	   numPadOnlyHit_t2++;
	 
	 if(hitvec[u].GetHitClass()==1){
	   hitvec[u].GetHitPhi();
	   numCl1HitHit_t2++;	   
	 }	   	     
       }else
	 {
	   if(hitvec[u].GetHitNumStrip()>0)
	     numStripOnlyHit_t2++;
	 }
       

       /*
       else
	 {
	   if((hitvec[k].GetHitNumStrip()>0)&&(hitvec[k].GetHitNumPad()==0)&&(Usestrip))
	     if(refphi>0)
	       {
		 //Found a fake pad coordinate and create an hit using only strip info.
		 // double phirad=hitvec[k].GetHitPhi()*3.14159/180.0;
		 // sigmaPhi=hitvec[k].GetHitDPhi()*3.14159/180.0;
		 double bestxydist=100.;
		 double xydist=0.;
		 double fkX=0;	  
		 double fkY=0;  
		 double bestphirad=0;
		 double thisphi=0;
		 somethingmatch=false;
		 
		 for(unsigned int step=0;step<30;step++)  
		   {
		     thisphi=(refphi+1.0*step)*3.14159/180.0;
		     fkX=r*cos(thisphi); fkY=r*sin(thisphi);//ax,bx.ay.by..
		     xydist=(VectTrkacking.at(0)*hitvec[k].GetHitZ()+VectTrkacking.at(1)-fkX)*(VectTrkacking.at(0)*hitvec[k].GetHitZ()+VectTrkacking.at(1)-fkX);
		     xydist=xydist+(VectTrkacking.at(2)*hitvec[k].GetHitZ()+VectTrkacking.at(3)-fkX)*(VectTrkacking.at(2)*hitvec[k].GetHitZ()+VectTrkacking.at(3)-fkX);
		     xydist=sqrt(xydist);
		     if(xydist<bestxydist)
		       {
			 bestxydist=xydist;
			 bestphirad=thisphi;
			 somethingmatch=true;
		       }
		   }
		 
		 for(unsigned int step=0;step<30;step++)  
		   {
		     thisphi=(refphi-1.0*step)*3.14159/180.0;
		     fkX=r*cos(thisphi); fkY=r*sin(thisphi);//ax,bx.ay.by..
		     xydist=(VectTrkacking.at(0)*hitvec[k].GetHitZ()+VectTrkacking.at(1)-fkX)*(VectTrkacking.at(0)*hitvec[k].GetHitZ()+VectTrkacking.at(1)-fkX);
		     xydist=xydist+(VectTrkacking.at(2)*hitvec[k].GetHitZ()+VectTrkacking.at(3)-fkX)*(VectTrkacking.at(2)*hitvec[k].GetHitZ()+VectTrkacking.at(3)-fkX);
		     xydist=sqrt(xydist);
		     if(xydist<bestxydist)
		       {
			 bestxydist=xydist;
			 bestphirad=thisphi;
			 somethingmatch=true;
		       }
		   }
		 
		 phirad=bestphirad;
		 sigmaPhi=sigmaPhi*3.0;//The scale factor is to check.
	       }
	   

	   
	   
	 }//else end
       */

       /*
       if(somethingmatch==false)
	 continue;
       */
       





     }//end for
   


 
   sizeHitv=hitvec2.size();   
 
  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      //   std::cout<<hitvec2[jj].GetHitX()<<" "<<hitvec2[jj].GetHitPhi()<<" "<<hitvec2[jj].GetHitZ()<<std::endl;
      //r->x   phi->y
       if(verbosity)
	 std::cout<<"X-Y-R-Z:"<<hitvec2[jj].GetHitX()<<" "<<hitvec2[jj].GetHitY()<<" "<<hitvec2[jj].GetHitR()<<" "<<hitvec2[jj].GetHitZ()<<std::endl;

      x.push_back(hitvec2[jj].GetHitX());
      y.push_back(hitvec2[jj].GetHitY());
      z.push_back(hitvec2[jj].GetHitZ());
      /* 
      double r=hitvec2[jj].GetHitR();
      double phirad=hitvec2[jj].GetHitPhi()*3.14159/180.0;
      double dphirad=hitvec2[jj].GetHitDPhi()*3.14159/180.0;
      */

      double sigmay;
      double sigmax;
      //er.push_back(0.1);
      
	
      sigmax=hitvec2[jj].GetHitDX();
      ex.push_back((sigmax));
      
      //ephi.push_back(0.1);
      sigmay=hitvec2[jj].GetHitDY();
      ey.push_back((sigmay));
	
      
     
      if((sigmax==0)||(sigmay==0))
	std::cout<<"ERROR in T2TrackProducer2::fitTrackspecialXY: sigmaX=  "<<sigmax<<" sigmaY= "<<sigmay<<std::endl;
      //ez.push_back(hitvec2[jj].GetHitDZ());
            



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
 
 if(verbosity)
    std::cout<<"AX:"<<a_xz<<"BX:"<< b_yz<<std::endl;

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
 /*
T1T2Track fittedtrack(vect,mat,chi2,chi2X,chi2Y,hemisphere,2);
 */
//std::cout<<"chi2X after: "<<fittedtrack.ChiSquaredX()<<"  chi2Y after: "<<fittedtrack.ChiSquaredY()<<std::endl;
//sizeHitv=trk.GetHitEntries();


 TVectorD FittedParam(4); //ax bx ay by
//Constructor utilized for T2 XY AND RZ fit



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

 par_covariance_matrix=ParCov_matr;
 /*

//  T1T2Track(const TVectorD & track_params_vector, const TMatrixD &par_covariance_matrix, double chiSquared,double chiSquaredX,double chiSquaredY,  int hemi, int det, double TanthetaRZ, double bRZ, double  phiRZ, double  e_TanthetaRZ, double e_bRZ, double e_phiRZ, double chi2R, double  chi2Phi, int theT2RoadID, unsigned int numHit_t2, unsigned int numCl1HitHit_t2, unsigned int numStripOnlyHit_t2, unsigned int numPadOnlyHit_t2);

inline double X0() const {return track_params_vector_[0];}
  inline double X0Sigma() const {return sqrt(par_covariance_vector_[0]);}
  inline double Y0() const {return track_params_vector_[1];}
  inline double Y0Sigma() const {return sqrt(par_covariance_vector_[1]);}
  inline double Z0() const {return z0_;}
  inline void Z0(double z0) {z0_=z0;}
  inline double GetTx() const {return track_params_vector_[2];}
  inline double GetTxSigma() const {return sqrt(par_covariance_vector_[2]);}
  inline double GetTy() const {return track_params_vector_[3];}
  inline double GetTySigma() const {return sqrt(par_covariance_vector_[3]);}

 */







T1T2Track trk2rz=RPhiFit(hitvec);//RZFIT  
  
  double chi2Phi=trk2rz.ChiSquaredPhi();
  double chi2R=trk2rz.ChiSquaredR();
  
  double phiRZ=trk2rz._phiRZ;               //trk2rz.Phi(); 
  double e_phiRZ=trk2rz._e_phiRZ;
  
  double TanthetaRZ=trk2rz._TanthetaRZ;     //trk2rz.GetTy(); 
  double e_TanthetaRZ=trk2rz._e_TanthetaRZ; //trk2rz.GetTySigma();
  
  double bRZ=trk2rz._bRZ;          //trk2rz.GetTx(); 
  double e_bRZ=trk2rz._e_bRZ;      //trk2rz.GetTxSigma();
  
  if((fabs(vect[2])>1.8)||(fabs(vect[3])>1.8))
    {
      std::cout<<"ERROR in T2TrackProducer3 Ax:"<<(vect[2])<<"  Ay:"<<(vect[3])<<" too big!! Trk Hit (all) print :"<<std::endl;
      for(unsigned int jj =0; jj<hitvec2.size(); jj++)
	{
	  std::cout<<"x-y-z-phi:  "<<hitvec2[jj].GetHitX()<<"  "<<hitvec2[jj].GetHitY()<<"  "<<hitvec2[jj].GetHitZ()<<"  "<<hitvec2[jj].GetHitPhi()<<"     DX:"<<hitvec2[jj].GetHitDX()<<"     DY:"<<hitvec2[jj].GetHitDY()<<std::endl;
	}
    }

  if(verbosity){
    std::cout<<"Before save: XY:"<<std::endl;
    //std::cout<<"ax.bx.ay.by:"<<vect[0]<<" "<<vect[1]<<" "<<vect[2]<<" "<<vect[3]<<std::endl;
  
    std::cout<<"Before save: RZ:"<<std::endl;
    std::cout<<"Phi "<<phiRZ*180.0/3.14159<<" TanTheta"<<TanthetaRZ<<" Brz:"<<bRZ<<" Eta:"<<trk2rz.Eta()<<std::endl;

    std::cout<<" TotHit:"<<numHit_t2<<" cl1Hit:"<<numCl1HitHit_t2<<" StripOnly:"<<numStripOnlyHit_t2<<" PadOnly:"<<numPadOnlyHit_t2<<std::endl;
  }


  T1T2Track fittedtrack(FittedParam, ParCov_matr,chi2,chi2,chi2,hemisphere,2 , TanthetaRZ,  bRZ,   phiRZ,   e_TanthetaRZ,  e_bRZ,  e_phiRZ, chi2R, chi2Phi, RoadID, numHit_t2, numCl1HitHit_t2, numStripOnlyHit_t2, numPadOnlyHit_t2);
  
  //T1T2Track fittedtrack(FittedParamConv,ParCov_matrConv,chi2,chi2,chi2,hemisphere,2 , TanthetaRZ,  bRZ,   phiRZ,   e_TanthetaRZ,  e_bRZ,  e_phiRZ, chi2R, chi2Phi, RoadID, numHit_t2, numCl1HitHit_t2, numStripOnlyHit_t2, numPadOnlyHit_t2);
  
  if(verbosity){      
    std::cout<<"After save: :"<<std::endl;
    std::cout<<"PhiRZ:"<<fittedtrack._phiRZ*180.0/3.14159<<"  EtaRZ:"<< fittedtrack._etaRZ<<"  EtaXY:"<<fittedtrack.Eta()<<" ChirR/N:"<<(chi2R/(numCl1HitHit_t2+numPadOnlyHit_t2))<<" Chi2/N; "<<chi2/(numCl1HitHit_t2+numPadOnlyHit_t2)<<" NumHit:"<<fittedtrack._numHit_t2<<std::endl;
    
  }




 sizeHitv=hitvec2.size();
for(unsigned int jj =0; jj<sizeHitv; jj++)
  {
    fittedtrack.AddHit(hitvec2[jj]);
  }




 //fittedtrack.AddHit(trk.GetHitT2(jj));
  // fittedtrack.AddHit(hitvec[jj]);

//  std::cout<<" Here 4 "<<std::endl;

   return fittedtrack; 
}





T1T2Track T2TrackProducer3::MyLinearfitCorrDEV(std::vector<T2Hit> hitvec2,TMatrixD &par_covariance_matrix,double &chi2_,bool Usestrip,std::vector<double>  VectTrkacking, int RoadID,T2Hit &worstHit)
{

  std::vector<T2Hit> hitvec; hitvec.clear();
   
 
  double sigmaR=(0.12+0.05);
  double sigmaPhi=0.015;
  unsigned int numphimin20=0;
  unsigned int numphimag340=0;
  double refphi=0;
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

  /////////////////////////////////////////////////////////////////////////////////
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
	    refphi=hitvec2.at(m).GetHitPhi(); 
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
  //A ?? la matrice per cui Mis=A(Param)
  TMatrixD Vy(sizeArighe,sizeArighe); //matrice di covarianza delle misure (una per ogni xy, quindi ?? diag a blocchi);
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
   ///         INCLUDE EVENTUALLY STRIPONLY 
   ///////////////////////////////////////////////////////////////////////////////

   bool somethingmatch=true;
   if((hitvec[k].GetHitNumStrip()>0)&&(hitvec[k].GetHitNumPad()==0)&&(Usestrip))
     {
       //Found a fake pad coordinate and create an hit using only strip info.
       // double phirad=hitvec[k].GetHitPhi()*3.14159/180.0;
       // sigmaPhi=hitvec[k].GetHitDPhi()*3.14159/180.0;
	double bestxydist=100.;double xydist=0.;
	double fkX=0;	    double bestphirad=0;double thisphi=0;
	somethingmatch=false;

	  for(unsigned int step=0;step<30;step++)  
	    {
	      thisphi=(refphi+1.0*step)*3.14159/180.0;
	      fkX=r*cos(thisphi); //ax,bx.ay.by..
	      xydist=(VectTrkacking.at(0)*hitvec[k].GetHitZ()+VectTrkacking.at(1)-fkX)*(VectTrkacking.at(0)*hitvec[k].GetHitZ()+VectTrkacking.at(1)-fkX);
	      xydist=xydist+(VectTrkacking.at(2)*hitvec[k].GetHitZ()+VectTrkacking.at(3)-fkX)*(VectTrkacking.at(2)*hitvec[k].GetHitZ()+VectTrkacking.at(3)-fkX);
	      xydist=sqrt(xydist);
	      if(xydist<bestxydist)
		{
		  bestxydist=xydist;
		  bestphirad=thisphi;
		  somethingmatch=true;
		}
	    }

	  for(unsigned int step=0;step<30;step++)  
	    {
	      thisphi=(refphi-1.0*step)*3.14159/180.0;
	      fkX=r*cos(thisphi); //ax,bx.ay.by..
	      xydist=(VectTrkacking.at(0)*hitvec[k].GetHitZ()+VectTrkacking.at(1)-fkX)*(VectTrkacking.at(0)*hitvec[k].GetHitZ()+VectTrkacking.at(1)-fkX);
	      xydist=xydist+(VectTrkacking.at(2)*hitvec[k].GetHitZ()+VectTrkacking.at(3)-fkX)*(VectTrkacking.at(2)*hitvec[k].GetHitZ()+VectTrkacking.at(3)-fkX);
	      xydist=sqrt(xydist);
	      if(xydist<bestxydist)
		{
		  bestxydist=xydist;
		  bestphirad=thisphi;
		  somethingmatch=true;
		}
	    }
	  
	  phirad=bestphirad;
	  sigmaPhi=sigmaPhi*3.0;//The scale factor is to check.
	
     }
  
   // if(swapupto0)
   //phirad=phirad-360 //Is not important working with the projecion

   if(somethingmatch==false)
     continue;

   if((hitvec[k].GetHitNumStrip()>0)&&(hitvec[k].GetHitNumPad()==0)&&(Usestrip==false))
     continue;


 
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

//Nota che qui Vy ?? in realt?? Vy^-1 

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


  
  T1T2Track trk2rz=RPhiFit(hitvec);//RZFIT  
  
  double chi2Phi=trk2rz.ChiSquaredPhi();
  double chi2R=trk2rz.ChiSquaredR();
  
  double phiRZ=trk2rz._phiRZ;               //trk2rz.Phi(); 
  double e_phiRZ=trk2rz._e_phiRZ;
  
  double TanthetaRZ=trk2rz._TanthetaRZ;     //trk2rz.GetTy(); 
  double e_TanthetaRZ=trk2rz._e_TanthetaRZ; //trk2rz.GetTySigma();
  
  double bRZ=trk2rz._bRZ;          //trk2rz.GetTx(); 
  double e_bRZ=trk2rz._e_bRZ;      //trk2rz.GetTxSigma();
  /*
  if((fabs(vect[0])>1.8)||(fabs(vect[2])>1.8))
    {
      std::cout<<"WARNING in T2TrackProducer3 Ax:"<<(vect[0])<<"  Ay:"<<(vect[2])<<" too big!! Trk Hit (all) print :"<<std::endl;
      for(unsigned int jj =0; jj<hitvec.size(); jj++)
	{
	  std::cout<<"x-y-z-phi:  "<<hitvec[jj].GetHitX()<<"  "<<hitvec[jj].GetHitY()<<"  "<<hitvec[jj].GetHitZ()<<"  "<<hitvec[jj].GetHitPhi()<<"     DX:"<<hitvec[jj].GetHitDX()<<"     DY:"<<hitvec[jj].GetHitDY()<<" NumPad:"<<hitvec[jj].GetHitNumPad()<<std::endl;
	}
    }
*/
  if(verbosity){
    std::cout<<"Before save: XY:"<<std::endl;
    // std::cout<<"ax.bx.ay.by:"<<vect[0]<<" "<<vect[1]<<" "<<vect[2]<<" "<<vect[3]<<std::endl;
  
    std::cout<<"Before save: RZ:"<<std::endl;
    std::cout<<"Phi "<<phiRZ*180.0/3.14159<<" TanTheta"<<TanthetaRZ<<" Brz:"<<bRZ<<" Eta:"<<trk2rz.Eta()<<std::endl;

    std::cout<<" TotHit:"<<numHit_t2<<" cl1Hit:"<<numCl1HitHit_t2<<" StripOnly:"<<numStripOnlyHit_t2<<" PadOnly:"<<numPadOnlyHit_t2<<std::endl;
  }


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





T1T2Track T2TrackProducer3::RPhiFit(std::vector<T2Hit> &hitvec2)
{

  T1T2Track trk(2);
  trk.SetDet(2);

  //std::cout<<"T1T2Track BEGIN T2SelectionCutUtils::RPhiFit"<<std::endl;

  // std::vector<T2Hit> hitvec; //Only Real Hits hits
  std::vector<T2Hit> hitvecR;  //Real Hits and Vtx
  

  for(unsigned int jj =0; jj<hitvec2.size(); jj++)
   {         
	 hitvecR.push_back(hitvec2.at(jj));
      
   }
  //std::cout<<"Inside MyLinearfitCorr .. init"<<std::endl;

  
  //  unsigned int  sizeHitv=hitvec.size();
  unsigned int  sizeHitvR=hitvecR.size();
  if(sizeHitvR<2)
   {
     std::cout<<" T2SelectionCutUtils::RPhiFit problem: Track with less than 2 Cl1 hits!! Continue the fitting.."<<std::endl;
     
   }


  //std::vector<T2Hit> hitvecraw;
  
  //std::vector<T2Hit> doublehitz;
  //std::vector<T2Hit> hits0360;


 unsigned int num0=0;
 unsigned int num360=0;
 

for(unsigned int jj =0; jj<sizeHitvR; jj++)
  { //std::cout<<"Hitphi "<<hitvec[jj].GetHitPhi()<<std::endl;
    //if((hitvecR[jj].GetHitPhi()>0.)&&(hitvecR[jj].GetHitPhi()<10.))//<110 MODIFY HERE
    unsigned int half=hitvecR[jj].GetHitHalftele();
    if((hitvecR[jj].GetHitPhi()<110)&&(half==0))  
       num0++;

    if((hitvecR[jj].GetHitPhi()>260)&&(half==0))//((hitvecR[jj].GetHitPhi()>350.)&&(hitvecR[jj].GetHitPhi()<360.))//>260
      num360++;
  }


 if((num0>0) && (num360>0)) 
   {
    
     for(unsigned int jj =0; jj<sizeHitvR; jj++)
       {	
	 if((hitvecR[jj].GetHitPhi()>0.)&&(hitvecR[jj].GetHitPhi()<10.))
	   hitvecR[jj].SetHitPhi(hitvecR[jj].GetHitPhi()+360.);
       }
   }



 int hemisphere = (int)(hitvecR[0].GetHitZ()/fabs(hitvecR[0].GetHitZ()));

 std::vector<float> r;
 std::vector<float> phi;
 std::vector<float> z;
 std::vector<float> er;
 std::vector<float> ephi;
 std::vector<float> ez;


float Sr=0;
float Srz=0;
float Szz_r=0;
float Sz_r=0; 
float S0_r=0; 
float Swphi=0;
float Sphiwphi=0;
 
 float cl1err=(0.12+0.05);
 
 //loop on cl1hits+vtx
  for(unsigned int jj =0; jj<sizeHitvR; jj++)
    {
      
      r.push_back(hitvecR[jj].GetHitR());
      
      z.push_back(hitvecR[jj].GetHitZ());    
    //  er.push_back(hitvecR[jj].GetHitDR());
      if(hitvecR[jj].GetHitClass()==1)
	 er.push_back(cl1err);
      else
	er.push_back(hitvecR[jj].GetHitDR());

      //	er.push_back(0.12+0.05);

      //ephi.push_back(hitvecR[jj].GetHitDPhi());
      
      ez.push_back(hitvecR[jj].GetHitDZ());
      
    
      
      Srz += r[jj]*z[jj]/er[jj]/er[jj];
      Szz_r += z[jj]*z[jj]/er[jj]/er[jj];
      Sz_r += z[jj]/er[jj]/er[jj];
      Sr += r[jj]/er[jj]/er[jj];
      S0_r += 1.0/er[jj]/er[jj];


      if(hitvecR.at(jj).GetHitClass()!=9){
	phi.push_back(hitvecR[jj].GetHitPhi());
	ephi.push_back(0.8);
	Swphi += 1.0/ephi.at(ephi.size()-1)/ephi.at(ephi.size()-1);
	Sphiwphi+=phi.at(phi.size()-1)/ephi.at(ephi.size()-1)/ephi.at(ephi.size()-1);
	//	std::cout<<phi.at(phi.size()-1)<<"  ->"<<Sphiwphi<<std::endl;
	  //Swphi += 1.0/ephi[jj]/ephi[jj];
	  //Sphiwphi+= phi[jj]/ephi[jj]/ephi[jj];
      }
  }


double a_rz = (Srz*S0_r - Sz_r*Sr) / (Szz_r*S0_r - Sz_r*Sz_r);   // angular coefficient
double b_rz = (Sr*Szz_r - Sz_r*Srz) / (Szz_r*S0_r - Sz_r*Sz_r);  // intercept   R=(a_rz)Z + b_rz
double e_a_rz = sqrt( S0_r / (S0_r*Szz_r - Sz_r*Sz_r) );         
double e_b_rz = sqrt( Szz_r / (S0_r*Szz_r - Sz_r*Sz_r) );
double phim=Sphiwphi/Swphi;
double e_phim= 1.0/sqrt(Swphi);

//std::cout<<"->->---->"<<phim<<std::endl;
double covab= - (Sz_r) / (Szz_r*S0_r - Sz_r*Sz_r);

double chi2r = 0.0;
double chi2p = 0.0;
double chi2= 0.0;
double normchi2red=0.0;


 unsigned int cl1count=0;

 for(unsigned int jj =0; jj<sizeHitvR; jj++)
   {
     if(hitvecR[jj].GetHitClass()==1)
       {
	 chi2r += (a_rz*hitvecR[jj].GetHitZ()+b_rz - hitvecR[jj].GetHitR())*(a_rz*hitvecR[jj].GetHitZ()+b_rz - hitvecR[jj].GetHitR())/cl1err/cl1err;
	 cl1count++;
       }
     else
       chi2r += (a_rz*hitvecR[jj].GetHitZ()+b_rz - hitvecR[jj].GetHitR())*(a_rz*hitvecR[jj].GetHitZ()+b_rz - hitvecR[jj].GetHitR())/hitvecR[jj].GetHitDR()/hitvecR[jj].GetHitDR();
       
   }



  for(unsigned int jjj =0; jjj<sizeHitvR; jjj++){          
      
      chi2p += (phim - phi[jjj])*(phim - phi[jjj])/ephi[jjj]/ephi[jjj];
    }


  
  if(cl1count>=2) 
    {      
      chi2=  (chi2p+chi2r)/2.0;
      normchi2red=(chi2p/((double)(cl1count)-1)+ chi2r/((double)(cl1count)-2))/2.0 ;
    }
  else
    {
      chi2=1.0;
      normchi2red=1.0;
      std::cout<<"Error: T2SelectionCutUtils RPhiFit without at least 2 hits"<<std::endl;
    }


  if(phim>360.0)                // This could happen when num0>0 AND num360>0
    phim=phim-360.0;
    

    phim=phim*3.14159265/180.0;
    e_phim=e_phim*3.14159265/180.0;

    //std::cout<<"rphifit: "<<phim*180/3.14159265<<std::endl;
    

    TVectorD vect(4);
    for(int oo=0; oo<4; oo++)
      vect[oo]=0.;
   
  
    vect[0]=phim;    // phi intercept
    vect[1]=0.;      // phi ang. coeff.
    vect[2]=b_rz;    // r intercept
    vect[3]=a_rz;    // r ang. coeff.
   

  
    TMatrixD mat(4,4);
    for(unsigned int oo=0; oo<4; oo++)
      for(unsigned int ooo=0;ooo<4; ooo++)
	mat[oo][ooo]=0.;


      mat[0][0]=e_phim*e_phim;
      mat[0][1]=0.;
      mat[0][2]=0.;
      mat[0][3]=0.;
      mat[1][0]=0.;
      mat[1][1]=0.;
      mat[1][2]=0.;
      mat[1][3]=0.;
      mat[2][0]=0.;
      mat[2][1]=0.;
      mat[2][2]=e_b_rz*e_b_rz;
      mat[2][3]=covab;
      mat[3][0]=0.;
      mat[3][1]=0.;
      mat[3][2]=covab;
      mat[3][3]=e_a_rz*e_a_rz;
    



//   double theeta = fabs(-log(tan(thetheta/2.))) * (double)hemisphere;

   //std::cout<<"calc phim "<<phim*180.0/3.14159265<<" "<<hitvec2[0].GetHitPhi()<<" "<<hitvec[0].GetHitPhi()<<std::endl;
   
    T1T2Track fittedtrack(vect,mat,a_rz, b_rz, phim, e_a_rz, e_b_rz, e_phim, chi2,chi2r,chi2p, normchi2red, hemisphere,2);   

    //I don't need here a final output trk.
    /*
    for(unsigned int jj =0; jj<hitvec2.size(); jj++)
      {
	//fittedtrack.AddHit(hitvec[jj]);
	fittedtrack.AddHit(hitvec2[jj]);
      }
    */


    return fittedtrack; 


}



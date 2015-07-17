/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Mirko Berretti (mirko.berretti@gmail.com)
*
****************************************************************************/

#include "TotemAnalysis/TotemNtuplizer/interface/T2Ntuplizer.h"

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"

#include "CLHEP/Vector/LorentzVector.h"

#include "DataFormats/T2Digi/interface/T2StripDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2StripDigi.h"
#include "DataFormats/T2Digi/interface/T2PadDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2PadDigi.h"
#include "DataFormats/T2DetId/interface/T2DetId.h"
#include "DataFormats/T2Hit/interface/T2Hit.h"
#include "DataFormats/T2Hit/interface/T2HitCollection.h"
#include "DataFormats/T2Cluster/interface/T2PadClusterCollection.h"
#include "DataFormats/T2Cluster/interface/T2StripClusterCollection.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TTree.h"
#include "TMath.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenEvent.h"

T2Ntuplizer::T2Ntuplizer(const edm::ParameterSet& iConfig) : Ntuplizer(iConfig)
{
    /*
  Chicut = iConfig.getParameter<double>("Chicut");
  PhiChiProbCut = iConfig.getParameter<double>("PhiChiProbCut");
  RChiProbCut = iConfig.getParameter<double>("RChiProbCut");  
  energy=iConfig.getParameter<double>("FitdzEnergy");
  DZScale=iConfig.getParameter<double>("DZScale");
  tracketamin=iConfig.getParameter<double>("TrkEtamin");
  tracketamax=iConfig.getParameter<double>("TrkEtaMAX");
  singleparticle=iConfig.getParameter<bool>("singleparticle");
  HepMCProductLabel = iConfig.getParameter<std::string>("HepMCProductLabel");
  */

  //vtxpos=iConfig.getParameter<double>("vtxpos");
  t2StripClusterCollectionLabel = iConfig.getParameter<edm::InputTag>("T2StripClusterCollectionLabel");
  t2PadClusterCollectionLabel = iConfig.getParameter<edm::InputTag>("T2PadClusterCollectionLabel");
  HitLabel = iConfig.getParameter<edm::InputTag>("HitLabel");
  RoadLabel = iConfig.getParameter<edm::InputTag>("RoadLabel");
  TrackLabel= iConfig.getParameter<edm::InputTag>("TrackLabel");
  t2PadDigiCollectionLabel = iConfig.getParameter<edm::InputTag>("T2PadDigiCollectionLabel");
  t2StripDigiCollectionLabel = iConfig.getParameter<edm::InputTag>("T2StripDigiCollectionLabel");

  //RawDataName= iConfig.getParameter<std::string>("RawDataName");

  // T2CutsUtil.SetCuts(T2_TrkEtamin,T2_TrkEtaMAX,T2_trkMultimin,T2_trkMultiMAX,T2_DZMultiplicator,T2_PhiChiProbCut,T2_RChiProbCut,T2_QuarterUsed,IgnoredSmallAnglePrimarySlopeCut,_T2_usesXYtracking);
  std::vector<int> qused;qused.push_back(0);qused.push_back(1);qused.push_back(2);qused.push_back(3);
  
  T2CutsUtil.SetCuts(4.5,7.5,4,11,2.,0.01,0.01,qused,0.001,true);


}


//DZ curve obtained with 50 GeV pion
double sigmaZfunctE50(double eta){
  double theDZ=0.0;
  double nepero=10.0;

  if(fabs(eta)<5.2)
    theDZ=3000.0;

  if(fabs(eta)>6.6)
    theDZ=4000.0;

  if((fabs(eta)>5.2)&&(fabs(eta)<6.6))
    theDZ= 7.79425*pow(nepero,6.0)-  5.16635*pow(nepero,6.0)*eta  +  1.28288*pow(nepero,6.0)*eta*eta  -1.41458*pow(nepero,5.0)*eta*eta*eta+ 5.84678*pow(nepero,3.0)*eta*eta*eta*eta;

  if((fabs(eta)>5.5)&&(fabs(eta)<5.6))
    theDZ=2000.0;

  return theDZ;
}



//DZ curve obtained with 10 GeV pion
double sigmaZfunctE10(double eta)
{
  double theDZ=0.0; 
  double nepero=10.0;
 

 if(fabs(eta)<5.37)
     theDZ=3200.0;

  if(fabs(eta)>6.25)
     theDZ=4400.0;


   if((fabs(eta)>5.37)&&(fabs(eta)<6.25))
     {     
       theDZ=1.09980*pow(nepero,5.0) -3.84187*pow(nepero,4.0)*eta+  3.44411*pow(nepero,3.0)*eta*eta;
     }

   //if((fabs(eta)>5.42)&&(fabs(eta)<5.78))
   if((fabs(eta)>5.5)&&(fabs(eta)<5.75))  
     theDZ=3300.0;
 

  return theDZ;
}


//DZ curve obtained with 10-50  GeV pion
double sigmaZfunctE30(double eta){

  double theDZ=0.0;
  double nepero=10.0;

  
  if(fabs(eta)<5.37)
     theDZ=2200.0;
  
   if(fabs(eta)>6.4)
     theDZ=4000.0;

   if((fabs(eta)>5.37)&&(fabs(eta)<6.4)){
       theDZ=1.11966*pow(nepero,5.0)-3.91045*pow(nepero,4.0)*eta+3.47461*pow(nepero,3.0)*eta*eta;
     }
   
 
   if((fabs(eta)>5.45)&&(fabs(eta)<5.7))  
     theDZ=2400.0;

  return theDZ;
}

unsigned int T2Ntuplizer::RawtoSymb(uint32_t thedet)
{
  T2DetId converter;
  unsigned int pl=converter.plane(thedet);
  unsigned int pls=converter.planeSide(thedet);
  unsigned int ht=converter.halfTelescope(thedet);
  unsigned int arm=converter.arm(thedet);
  unsigned int symbolic=pl*2+pls+ht*10+20*arm;	  
 
  return symbolic;
}

double T2Ntuplizer::EtaFromAveragesVtxCorr(T1T2Track thetrk,double vtx_x,double vtx_y,double vtx_z)
{

double thetaRHitAverage=0.;double hemis=0.;
  double counter=0.;

  if(thetrk.GetHitT2(1).GetHitClass()!=9)
    hemis=thetrk.GetHitT2(1).GetHitZ()/fabs(thetrk.GetHitT2(1).GetHitZ());
  else
    std::cout<<"Error in T2SelectionCutUtils::EtaFromAverages"<<std::endl;

  for(unsigned int mm=0;mm<thetrk.GetHitEntries();mm++)
    {
      // hitvector.push_back(thetrk.GetHitT2(mm));
       if(thetrk.GetHitT2(mm).GetHitClass()!=9)
	 {

	   double realR= (thetrk.GetHitT2(mm).GetHitX()-vtx_x)*(thetrk.GetHitT2(mm).GetHitX()-vtx_x)+(thetrk.GetHitT2(mm).GetHitY()-vtx_y)*(thetrk.GetHitT2(mm).GetHitY()-vtx_y);
	   realR=sqrt(realR);
	   double realZ=(thetrk.GetHitT2(mm).GetHitZ()-vtx_z);
	   // std::cout<<"RecTrkHit-Eta:"<<(hemis*(-1.0)*log(fabs(realR / realZ)/2.0))<<std::endl;
	   thetaRHitAverage += fabs(realR / realZ);
	   counter=counter+1.0;
	 }
    }
  
  thetaRHitAverage=thetaRHitAverage/counter;
  double TrkEtaFromHitAverages_4Hits = (hemis*(-1.0)*log(thetaRHitAverage/2.0)); //low angle approximation
  
 return TrkEtaFromHitAverages_4Hits;

}

bool IsinT2(double eta, double tracketamin, double tracketamax){
  bool flag;
  if ((fabs(eta)>tracketamin)&&(fabs(eta)<tracketamax))
    flag=true;
  else
    flag=false;
  return flag;
}

//
// member functions
//
// ------------ method called to for each event  ------------
void T2Ntuplizer::FillEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  /* LOADING OF ALL THE RECORDS FROM THE EVENT */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */


    /* :::::::::::::TakeDigi::::::::::::*/

  edm::Handle<T2PadDigiCollection> t2paddigicoll;
  iEvent.getByLabel(t2PadDigiCollectionLabel, t2paddigicoll);
  t2paddigicoll.product();
//  const T2PadDigiCollection* PadDigiptr;

  edm::Handle<T2StripDigiCollection> t2stripdigicoll;
  iEvent.getByLabel(t2StripDigiCollectionLabel, t2stripdigicoll);
  t2stripdigicoll.product();
//  const T2StripDigiCollection* StripDigiptr;


  DigiContainerIterator<T2DetId, T2PadDigi> itp;
  DigiContainerIterator<T2DetId, T2StripDigi> its;

  /* :::::::::::::Take The Clusters::::::::::::*/

  //std::cout << " check --------- A " <<  std::endl;
  Handle<T2StripClusterCollection> t2strclcoll;
  iEvent.getByLabel(t2StripClusterCollectionLabel,t2strclcoll);

  Handle<T2PadClusterCollection> t2padclcoll;
  iEvent.getByLabel(t2PadClusterCollectionLabel,t2padclcoll);

/*::::::Take  T2  Hits::::::*/
  edm::Handle<T2HitCollection> t2hitcoll;
  iEvent.getByLabel(HitLabel, t2hitcoll);


  //:::::: Take T2 tracks ::::::
  Handle<T1T2TrackCollection> trackCollection;
  iEvent.getByLabel(TrackLabel, trackCollection);
  
 
 
 //****************************************************
 //***************************************************
 
 //Fill the T2 ntple

 t2obj.ev_no = iEvent.id().event();
 t2obj.run_no =iEvent.id().run();
 t2obj.timestamp = iEvent.time().unixTime();



 t2obj.Pad_row.clear();             //strip row   
 t2obj.Pad_col.clear();             //strip column   
 t2obj.Pad_det.clear();             //symbolic id of the detector containind the pad
 
 t2obj.Strip_row.clear();           //strip row
 t2obj.Strip_col.clear();           //strip column
 t2obj.Strip_det.clear();           //symbolic id of the detector containind the strip

 
 t2obj.TrkAx.clear();            // slope of the track projection in the XZ plane
 t2obj.TrkAy.clear();            // slope of the track projection in the YZ plane
 t2obj.TrkX0.clear();            // X at Z=0 for the XZ projected track  
 t2obj.TrkY0.clear();            // Y at Z=0 for the XZ projected track
 t2obj.TrkPhi.clear();           // Trk Phi (XY fit)
 t2obj.TrkEta_XY.clear();
 t2obj.TrkZmin_XY.clear();
 t2obj.TrkRmin_XY.clear();
 t2obj.TrkZImpact.clear();
 t2obj.TrkChi2XProb.clear();          //Chi2-X probability (goodnes of the XZ projection fit)
 t2obj.TrkChi2YProb.clear();          //Chi2-Y probability (goodnes of the YZ projection fit)
 t2obj.TrkClass1HitCounter.clear();   //Number of class1 Hit in the Trk
 t2obj.TrkHitCounter.clear();         //Number of class1 + cluster Pad hits in the Trk

 
 t2obj.TrkThetaR_RZFit.clear();// Trk Polar angle obtained tracking in the Rz plane
 t2obj.TrkEta_RZFit.clear();   // Trk Eta obtained tracking in the Rz plane 
 t2obj.TrkPhi_RZFit.clear();   // Trk Phi obtained with a constant fit.
 t2obj.TrkZ0_RZFit.clear();    // Crossing Point between Trk and Z Axis obtained tracking in the Rz plane 
 t2obj.TrkBX_RZFit.clear();    // X0 @ Z=0 obtained tracking in the Rz plane 
 t2obj.TrkBY_RZFit.clear();    // X0 @ Z=0 obtained tracking in the Rz plane 


 t2obj.HitID.clear(); 
 t2obj.HitPhi.clear();      // Phi position of all the Hits (deg)
 t2obj.HitR.clear();        // R position of all the Hits (mm)
 t2obj.HitType.clear();     // 0-> only pad; 1-> only strip 2->Class 1 Hit (superimposition Pad/Strip)
 t2obj.HitNumPad.clear();   // Cluster Pad Size 
 t2obj.HitNumStrip.clear(); // Cluster Strip Size
 

 t2obj.TrkEntryX.clear();   //Trk Entry point X
 t2obj.TrkEntryY.clear();   //Trk Entry point Y 
 t2obj.TrkEntryZ.clear();   //Trk Entry point Z
 t2obj.TrkExitX.clear();    //Trk Exit point X
 t2obj.TrkExitY.clear();    //Trk Exit point Y  
 t2obj.TrkExitZ.clear();    //Trk Exit point Z

 t2obj.NumPadCluH0=0;
 t2obj.NumPadCluH1=0;
 t2obj.NumPadCluH2=0;
 t2obj.NumPadCluH3=0;

 t2obj.TrkNumHitInH0.clear();  
 t2obj.TrkNumHitInH1.clear();  
 t2obj.TrkNumHitInH2.clear();  
 t2obj.TrkNumHitInH3.clear();  

 t2obj.TrkEta2.clear();   
 
 t2obj.ProbChi2R_rz.clear();
 t2obj.Chi2Rreduced_rz.clear();
 t2obj.TrkChiProb.clear();
 t2obj.TrkPrimaryGeneratorEta.clear();

 
 T2GeometryUtil::T2DetInfo planeinfo;
 T2GeometryUtil conv;
  
 /*
     for(itp= PadDigiptr->begin(); itp!=PadDigiptr->end(); ++itp)
       {
	 T2DetId mydet=(*itp).first; 
	 int detnumb = mydet.arm()*20 + mydet.halfTelescope()*10 + mydet.plane()*2 + mydet.planeSide();

	 for(std::vector<T2PadDigi>::const_iterator itpad =(*itp).second.first; itpad !=(*itp).second.second; itpad++)
	   {
	     t2obj.Pad_row.push_back((*itpad).getRow());
	     t2obj.Pad_col.push_back((*itpad).getCol());
	     t2obj.Pad_det.push_back(detnumb);
	   }

       }
 */

 T1T2TrackCollection::const_iterator TrkCit;
 	 
	 
     double Zimpact=0.;
     
     for(TrkCit=trackCollection->begin(); TrkCit!=trackCollection->end(); TrkCit++){
        std::vector<T2Hit> hitvector;	
       unsigned int class1Hitcounter=0;
       unsigned int numh0=0;
       unsigned int numh1=0;
       unsigned int numh2=0;
       unsigned int numh3=0;
	 for(unsigned int i=0; i<(*TrkCit).GetHitEntries();i++)
	   { 

	     uint32_t rawiddet=(*TrkCit).GetHitT2(i).GetHitDetRawId();
	       //unsigned int  symb=RawtoSymb(rawiddet);
	       unsigned int symb=10;

	       planeinfo=conv.GetT2Info(rawiddet);
	       symb=planeinfo.symb;
	       if(symb<10)
		 numh0++;
	       else{
		 if(symb<20)
		   numh1++;
		 else{
		   if(symb<30)
		     numh2++;
		   else
		     numh3++;
		 }
	       }
	     if((*TrkCit).GetHitT2(i).GetHitClass()==1)
	       class1Hitcounter++;
	     hitvector.push_back((*TrkCit).GetHitT2(i));
	   }
	 
	 t2obj.TrkNumHitInH0.push_back(numh0);
	 t2obj.TrkNumHitInH1.push_back(numh1);
	 t2obj.TrkNumHitInH2.push_back(numh2);
	 t2obj.TrkNumHitInH3.push_back(numh3);

	 double eta2=EtaFromAveragesVtxCorr((*TrkCit),0.,0.,0.);// 
	 t2obj.TrkEta2.push_back(eta2);
	  
	 t2obj.TrkHitCounter.push_back((*TrkCit).GetHitEntries());
	 t2obj.TrkClass1HitCounter.push_back(class1Hitcounter);
	 
	 t2obj.TrkEntryX.push_back(TrkCit->GetHitT2(0).GetHitX());
	 t2obj.TrkEntryY.push_back(TrkCit->GetHitT2(0).GetHitY());
	 t2obj.TrkEntryZ.push_back(TrkCit->GetHitT2(0).GetHitZ());
	 t2obj.TrkExitX.push_back(TrkCit->GetHitT2((*TrkCit).GetHitEntries()-1).GetHitX());
	 t2obj.TrkExitY.push_back(TrkCit->GetHitT2((*TrkCit).GetHitEntries()-1).GetHitY());
	 t2obj.TrkExitZ.push_back(TrkCit->GetHitT2((*TrkCit).GetHitEntries()-1).GetHitZ());

	
	 t2obj.TrkAx.push_back(TrkCit->GetTx());
	 t2obj.TrkAy.push_back(TrkCit->GetTy());
	 t2obj.TrkX0.push_back(TrkCit->X0());
	 t2obj.TrkY0.push_back(TrkCit->Y0());
	 t2obj.TrkPhi.push_back(TrkCit->Phi()*180/3.14159);
	 t2obj.TrkEta_XY.push_back(TrkCit->Eta());
	 t2obj.TrkZmin_XY.push_back(TrkCit->Z_at_Rmin());
	 t2obj.TrkRmin_XY.push_back(TrkCit->Rmin());



	 double C0=((*TrkCit).GetHitT2(0).GetHitX()*(*TrkCit).GetHitT2(0).GetHitX()+(*TrkCit).GetHitT2(0).GetHitY()*(*TrkCit).GetHitT2(0).GetHitY())/((*TrkCit).GetTx()*(*TrkCit).GetHitT2(0).GetHitX()+(*TrkCit).GetTy()*(*TrkCit).GetHitT2(0).GetHitY());
	 Zimpact=(*TrkCit).GetHitT2(0).GetHitZ()-C0;

	 
	 t2obj.TrkZImpact.push_back(Zimpact);
	 

	 
	 //R0xy=sqrt(((*TC_iter).bx_)*((*TC_iter).bx_)+((*TC_iter).by_)*((*TC_iter).by_));
	 //R0rz=sqrt(((*TC_iter).bx_)*((*TC_iter).bx_)+((*TC_iter).by_)*((*TC_iter).by_)); 
	 //Z0xy=(*TC_iter).Z_at_Rmin();
	 // Z0rz=(*TC_iter).Z_at_Rmin();
	 double ProbChi2X_xy=TMath::Prob((*TrkCit).ChiSquaredX(),((*TrkCit).GetHitEntries()-2));
	 double ProbChi2Y_xy=TMath::Prob((*TrkCit).ChiSquaredY(),((*TrkCit).GetHitEntries()-2));

	 double ProbChi2R_rz=TMath::Prob((*TrkCit).ChiSquaredR() , ((*TrkCit).GetHitEntries()-2));
	 double Chi2Rreduced_rz=(*TrkCit).ChiSquaredR()/((*TrkCit).GetHitEntries()-2);
	 
	 t2obj.ProbChi2R_rz.push_back(ProbChi2R_rz);
	 t2obj.Chi2Rreduced_rz.push_back(Chi2Rreduced_rz);
	 //double chiprob=0.;
	 unsigned int numhitacl1=(*TrkCit)._numCl1HitHit_t2;
	 double ProbChi2_XY=TMath::Prob((*TrkCit).ChiSquared(),(numhitacl1*2-4));
	 t2obj.TrkChiProb.push_back(ProbChi2_XY);
	 
	 //ProbChi2Phi_rz=TMath::Prob((*TC_iter).ChiSquaredPhi(),((*TC_iter).GetHitEntries()-1));
	 
	 t2obj.TrkChi2YProb.push_back(ProbChi2Y_xy);
	 t2obj.TrkChi2XProb.push_back(ProbChi2X_xy);
	    

	 T1T2Track trk2rz=T2CutsUtil.TrackFromHits(false,hitvector);//RZFIT

	 t2obj.TrkZ0_RZFit.push_back(trk2rz.Z_at_Rmin());
	 t2obj.TrkPhi_RZFit.push_back(trk2rz.Phi()*180/3.14159); 
	 t2obj.TrkEta_RZFit.push_back(trk2rz.Eta());
	 t2obj.TrkThetaR_RZFit.push_back(trk2rz.GetTy());
	 t2obj.TrkBX_RZFit.push_back(trk2rz.bx_);
	 t2obj.TrkBY_RZFit.push_back(trk2rz.by_);

	 int ActivePlH0[10]= {0,0,0,0,0,0,0,0,0,0};
	 int ActivePlH1[10]= {0,0,0,0,0,0,0,0,0,0};
	 int ActivePlH2[10]= {0,0,0,0,0,0,0,0,0,0}; 
	 int ActivePlH3[10]= {0,0,0,0,0,0,0,0,0,0};

	 unsigned int intplane=0;
	 unsigned int totHitAllH0=0;unsigned int totHitAllH1=0;unsigned int totHitAllH2=0;unsigned int totHitAllH3=0;
  
	 for(T2PadClusterCollection::const_iterator itpad= t2padclcoll->begin(); itpad != t2padclcoll->end(); itpad++){//Reft2padclcoll
	   vector<T2Cluster> padClv = itpad->second;
	   T2DetId *detID =new T2DetId(itpad->first);
	   
	   uint32_t cmsswdId= detID->calculateRawId(detID->arm(),detID->halfTelescope(),detID->plane(),detID->planeSide());

	   delete detID;

	   unsigned int symb=RawtoSymb(cmsswdId);
	   //unsigned int quarter=symbol/10;
	   //unsigned int plane=symbol%10;
	   
	   if((symb/10)==0){
	     totHitAllH0=totHitAllH0+padClv.size();
	     intplane=(symb%10);
	     ActivePlH0[intplane]=ActivePlH0[intplane]+padClv.size();
	   }
	   if((symb/10)==1){	
	     totHitAllH1=totHitAllH1+padClv.size();
	     intplane=(symb%10);
	     ActivePlH1[intplane]=ActivePlH1[intplane]+padClv.size();
	     
	   }
	   if((symb/10)==2){
	     totHitAllH2=totHitAllH2+padClv.size();
	     intplane=(symb%10);
	     ActivePlH2[intplane]=ActivePlH2[intplane]+padClv.size();
	   }
	   if((symb/10)==3){
	     totHitAllH3=totHitAllH3+padClv.size();     
	     intplane=(symb%10);
	     ActivePlH3[intplane]=ActivePlH3[intplane]+padClv.size();
	   }
	 }
	 
  
	 t2obj.NumPadCluH0=totHitAllH0; 
	 t2obj.NumPadCluH1=totHitAllH1; 
	 t2obj.NumPadCluH2=totHitAllH2; 
	 t2obj.NumPadCluH3=totHitAllH3;
	 
	 //TrkEtaALL->Fill(trk2rz.Eta());

       
     }
 

  
}
     




/*

//----------------------------------------------------------------\\
//----------------- Note about the histograms---------------------\\
//----------------------------------------------------------------\

The meaningfull histograms for one particle events are:

  clusterstripentries
  clusterpadentries
  diffphiGUNHIT
  DPhiGoodTrk
  DEtaGoodTrk
  RLocHRecoH
  diffRCluHit
  diffphiCluGun
  diffphiCluHit
  SingleParticleEfficiencyAllCuts
  SingleParticleEfficiency


The meaningfull histograms for one particle or multiparticle events are:
  Trketagood
  Trkphigood
  Trketa
  Trkphi
  MultiParticleEfficiencyCutNorm
  MultiParticleEfficiencyCut
  Chi2PhiProb
  Chi2RProb
  Chi2PhiProbLogy
  Chi2RProbLogy

*/

double T2Ntuplizer::PartCharge(int pdgcode){
  double carica=0;
 switch(pdgcode){
 case(11):
   carica=-1; //electron
   break;
case(-11):
   carica=1; 
   break;

 case(13):
	carica=-1; //muon
	break;
case(-13):
	carica=+1; 
	break;
 case(15):
	carica=-1; //tau
	break;
case(-15):
	carica=+1; //tau+
	break;

 case(211):
	carica=1; //picharg+
	break;
 case(-211):
	carica=-1; //picharg+
	break;

 case(213):
	carica=1; //rhocharg+
	break;
 case(-213):
	carica=-1; //rhocharg-
	break;

 case(321):
	carica=1; //kcharg+
	break;

 case(-321):
	carica=-1; //kcharg-
	break;

 case(2212):
	carica=+1; //proton
	break;
 case(-2212):
	carica=-1; //antiproton
	break;

 case(1114):
	carica=-1; //delta-
	break;
 case(2214):
	carica=-1; //delta+
	break;
 case(2224):
	carica=-1; //delta++
	break;
 case(411):
	carica=1; //D+
	break;
case(-411):
	carica=-1; //D-
	break;
case(521):
	carica=1; //B+
	break;
case(-521):
	carica=-1; //B-
	break;
 }
 return carica;
}

// ------------ method called once each job just before starting event loop  ------------
void T2Ntuplizer::CreateBranches(const edm::EventSetup&, TTree *tree)
{
  TH1::AddDirectory(kFALSE);
 
  R0distro=std::auto_ptr<TH1F>(new TH1F("R0distro","Track r0 only chi2 cut",32,-100,100));
  R0distro->SetDirectory(0);			       
  Z0distro=std::auto_ptr<TH1F>(new TH1F("Z0distro","Track Z at r min only chi2 cut",50,-8000,8000));
  Z0distro->SetDirectory(0);
  NumhitinTrack=std::auto_ptr<TH1F>(new TH1F("NumhitinTrack","# hit in Track",18,-0.5,17.5));
  NumhitinTrack->SetDirectory(0);
  NumhitinTrackGood=std::auto_ptr<TH1F>(new TH1F("NumhitinTrackGood","# hit in Track",18,-0.5,17.5));
  NumhitinTrackGood->SetDirectory(0);
  
  chPartinT2=std::auto_ptr<TH1F>(new TH1F("chPartinT2","# Charged Particle in T2",16,-0.5,15.5));
  chPartinT2->SetDirectory(0);
  HNumrecotrackcutmag0=std::auto_ptr<TH1F>(new TH1F("HNumrecotrackcutmag0std","# Track reconstructed, #eta >0, all cut",16,-0.5,15.5));
  HNumrecotrackcutmag0->SetDirectory(0);
  HTrackcounter=std::auto_ptr<TH1F>(new TH1F("HTrackcounter","# Track reconstructed",26,-0.5,25.5));
  HTrackcounter->SetDirectory(0);
  HTrackcounterCrit=std::auto_ptr<TH1F>(new TH1F("HTrackcounterCrit","# Track reconstructed 5.5<#eta Ch< 5.6",10,-0.5,9.5));
  HTrackcounterCrit->SetDirectory(0);
  HTrackcounterNonCrit=std::auto_ptr<TH1F>(new TH1F("HTrackcounterNonCrit","# Track reconstructed 5.7<#eta Ch< 6.4",10,-0.5,9.5));     
  HTrackcounterNonCrit->SetDirectory(0);
   
  tantheta=std::auto_ptr<TH1F>(new TH1F("tantheta","tantheta, arz",100,-0.05,0.1));     
  tantheta->SetDirectory(0);

  tanthetam0=std::auto_ptr<TH1F>(new TH1F("tanthetam0","tantheta, arz",100,-0.05,0.1));     
  tanthetam0->SetDirectory(0);

  RecoZHit = std::auto_ptr<TH1F>(new TH1F("RecoZHit","Reco Hit Z",2000, 13500, 14500));
  RecoZHit->SetDirectory(0);
  
  SimuZHit = std::auto_ptr<TH1F>(new TH1F("SimuZHit","Simulated Hit Z",2000, 13500, 14500));
  SimuZHit->SetDirectory(0);
  // Chi2PhiProbLogy= std::auto_ptr<TCanvas>(new TCanvas("Chi2PhiProblog","Azimuthal #chi^{2} probability",400,400));
  // Chi2PhiProbLogy->SetDirectory(0);
  //Chi2RProbLogy= std::auto_ptr<TCanvas>(new TCanvas("Chi2RProblog","Radial #chi^{2} probability",400,400));
  //Chi2RProbLogy->SetDirectory(0);

  clusterstripentries = std::auto_ptr<TH1F>(new TH1F("clusterstripentries","Number of strips in strip-clusters",10, -0.5, 9.5));
  clusterstripentries->SetDirectory(0); 
  clusterpadentries = std::auto_ptr<TH1F>(new TH1F("clusterpadentries","Number of pads in pad-clusters",10, -0.5, 9.5));
  clusterpadentries->SetDirectory(0); 
  diffphiGUNHIT= std::auto_ptr<TH1F>(new TH1F("diffphiGUNHIT","#Phi Generated particle - #Phi Reconstructed Hit",120, -30, 30));
  diffphiGUNHIT->SetDirectory(0); 
  RLocHRecoH= std::auto_ptr<TH1F>(new  TH1F("RLocHRecoH"," Reconstructed Hit R - GEANT4 Hit R",100,-1.5, 1.5));
  RLocHRecoH->SetDirectory(0); 
  diffRCluHit= std::auto_ptr<TH1F>(new  TH1F("diffRCluHit"," Cluster R - GEANT4 Hit R",100,-1.5, 1.5));
  diffRCluHit->SetDirectory(0); 
  SourceofSecondary= std::auto_ptr<TH1F>(new TH1F("SourceofSecondary","#eta which generate secondaries",1800,-9.0,9.0));
  SourceofSecondary->SetDirectory(0);
  SourceofReco= std::auto_ptr<TH1F>(new TH1F("SourceofReco","#eta which generate reconstruction",1800,-9.0,9.0));
  SourceofReco->SetDirectory(0);
  SymbIdHT0= std::auto_ptr<TH1F>(new TH1F("SymbIdHT0","Symb det id",21,-0.5,19.5));
  SymbIdHT0->SetDirectory(0);
  TrkEtaGoodH0= std::auto_ptr<TH1F>(new TH1F("TrkEtaGoodH0","Reconstructed Track #eta (all cut used)",1800,-9.0,9.0));
  TrkEtaGoodH0->SetDirectory(0); 
  TrkEtaGoodH1= std::auto_ptr<TH1F>(new TH1F("TrkEtaGoodH1","Reconstructed Track #eta (all cut used)",1800,-9.0,9.0));
  TrkEtaGoodH1->SetDirectory(0); 
  TrkEtaH0= std::auto_ptr<TH1F>(new TH1F("TrkEtaH0","Reconstructed Track #eta (no cut used)",1800,-9.0,9.0));
  TrkEtaH0->SetDirectory(0); 
  TrkEtaH1= std::auto_ptr<TH1F>(new TH1F("TrkEtaH1","Reconstructed Track #eta (no cut used)",1800,-9.0,9.0));
  TrkEtaH1->SetDirectory(0); 

  Trketa= std::auto_ptr<TH1F>(new TH1F("Trketa","Reconstructed Track #eta (no cut used)",1800,-9.0,9.0));
  Trketa->SetDirectory(0); 
  Trkphi= std::auto_ptr<TH1F>(new TH1F("Trkphi","Reconstructed Track #phi (no cut used)",361,0,360));
  Trkphi->SetDirectory(0); 
  DPhiGoodTrk= std::auto_ptr<TH1F>(new TH1F("DPhiGoodTrk","#phi Track - #phi generated particle (Only primary tracks are considered)",51,-25.5,25.5));
  DPhiGoodTrk->SetDirectory(0); 
  DEtaGoodTrk= std::auto_ptr<TH1F>(new TH1F("DEtaGoodTrk","#eta Track - #eta generated particle (Only primary tracks are considered)",200,-1.0,1.0));
  DEtaGoodTrk->SetDirectory(0); 
  diffphiCluGun = std::auto_ptr<TH1F>(new TH1F("diffphiCluGun","Cluster #Phi - #Phi generated particle",81, -40.5, 40.5));
  diffphiCluGun->SetDirectory(0); 
  diffphiCluHit= std::auto_ptr<TH1F>(new TH1F("diffphiCluHit","Cluster #Phi - GEANT4 Hit #Phi",50, -5, 5));
  diffphiCluHit->SetDirectory(0); 
  Trketagood= std::auto_ptr<TH1F>(new TH1F("Trketagood","Reconstructed Track #eta (all cut used)",1800,-9.0,9.0));
  Trketagood->SetDirectory(0); 
  Trkphigood= std::auto_ptr<TH1F>(new TH1F("Trkphigood","Reconstructed Track #phi (all cut used)",361,0,360));
  Trkphigood->SetDirectory(0); 

  TrKPhi0degstudyg= std::auto_ptr<TH1F>(new TH1F("TrKPhi0degstudyg","trk phi",720, -270.25, 89.75));
  TrKPhi0degstudyg->SetDirectory(0); 

  DEtaChiCutOneTrk= std::auto_ptr<TH1F>(new TH1F("DEtaChiCutOneTrk","#eta Track - #eta generated particle, chi cut",200,-1.0,1.0));
  DEtaChiCutOneTrk->SetDirectory(0); 
  DPhiChiCutOneTrk= std::auto_ptr<TH1F>(new TH1F("DPhiChiCutOneTrk","#phi Track - #phi generated particle, chi cut ",51,-25.5,25.5));
  DPhiChiCutOneTrk->SetDirectory(0); 

  Chi2RProb = std::auto_ptr<TH1F>(new TH1F("Chi2RProb","Radial #chi^{2} probability",105,0.0,1.05));
  Chi2RProb->SetDirectory(0); 
  Chi2PhiProb = std::auto_ptr<TH1F>(new TH1F("Chi2PhiProb","Azimuthal #chi^{2} probability",105,0.0,1.05));
  Chi2PhiProb->SetDirectory(0);

  chiRdistro = std::auto_ptr<TH1F>(new TH1F("chiRdistro","Radial #chi^{2} ",105,0.0,100));
  chiRdistro->SetDirectory(0);
  chiPhidistro = std::auto_ptr<TH1F>(new TH1F("chiPhidistro","Azimuthal #chi^{2}",105,0.0,100));
  chiPhidistro->SetDirectory(0);


  stablepdg = std::auto_ptr<TH1F>(new TH1F("stablepdg","Pdg Id stable particles",200001,-10000.5,10000.5));
  stablepdg->SetDirectory(0);
  energyP04575= std::auto_ptr<TH1F>(new TH1F("energyP04575","#pi^{0} energy",2000,0.,2000.));
  energyP04575->SetDirectory(0);
  etaP04575= std::auto_ptr<TH1F>(new TH1F("etaP04575","#pi^{0} #eta",60,4.5,7.5));
  etaP04575->SetDirectory(0);
  
  P0EtaEnergy= std::auto_ptr<TH2F>(new TH2F("P0EtaEnergy","#eta vs Energy #pi^{0}",30,4.5,7.5,40,1.,400.));
  P0EtaEnergy->SetDirectory(0);

  ChOutEtaEnergy= std::auto_ptr<TH2F>(new TH2F("ChOutEtaEnergy","#eta vs Energy CH PARTICLE OUTSIDE T2",30,4.5,8.5,40,1.,400.));
  ChOutEtaEnergy->SetDirectory(0);



  ChEnergyinT2= std::auto_ptr<TH1F>(new TH1F("ChEnergyinT2","Charged particle energy in T2",100,0,1000));
  ChEnergyinT2->SetDirectory(0);



  muEtaEnergy= std::auto_ptr<TH2F>(new TH2F("muEtaEnergy","#eta vs Energy muons",30,4.5,8.5,40,1.,400.));
  muEtaEnergy->SetDirectory(0);


  P0EtaEnergycorr= std::auto_ptr<TProfile>(new TProfile("P0EtaEnergycorr","#eta vs Energy #pi^{0}",60,4.5,7.5,""));
  P0EtaEnergycorr->SetDirectory(0);

  P0NumEtacorr= std::auto_ptr<TProfile>(new TProfile("P0NumEtacorr"," # #pi^{0} VS #pi^{0} <#eta>",21,-0.5,20.5,""));
  P0NumEtacorr->SetDirectory(0);

  

  P0NumEnergycorr= std::auto_ptr<TProfile>(new TProfile("P0NumEnergycorr"," # #pi^{0} VS #pi^{0} <energy>",21,-0.5,20.5,""));
  P0NumEnergycorr->SetDirectory(0);

  SingleParticleEfficiencyCutEta= std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiencyCutEta","# Tracks reconstructed vs generated particle #eta (cut #eta)",34,4.35,7.75,""));//Prima l'estr. era a 7.75
  SingleParticleEfficiencyCutEta->SetDirectory(0);
  SingleParticleEfficiencyCutDZ= std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiencyCutDZ","# Tracks reconstructed vs generated particle #eta (cut DZ)",34,4.35,7.75,""));
  SingleParticleEfficiencyCutDZ->SetDirectory(0);
  SingleParticleEfficiencyCutChi= std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiencyCutChi","# Tracks reconstructed vs generated particle #eta (cut #c)",34,4.35,7.75,""));
  SingleParticleEfficiencyCutChi->SetDirectory(0);

  SingleParticleEfficiencyAllCuts = std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiencyAllCuts","# Tracks reconstructed vs generated particle #eta (all cuts)",34,4.35,7.75,""));
  SingleParticleEfficiencyAllCuts->SetDirectory(0); 
  SingleParticleEfficiency = std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiency","# Tracks reconstructed vs generated particle #eta (all track)",34,4.35,8.75,""));
  SingleParticleEfficiency->SetDirectory(0); 
  MultiParticleEfficiencyCut= std::auto_ptr<TProfile>(new TProfile("MultiParticleEfficiencyCut","<# Tracks> vs # Charged Particles (all track cuts)",12,-0.5,11.5,""));
  MultiParticleEfficiencyCut->SetDirectory(0); 
  MultiParticleEfficiencyCutNorm= std::auto_ptr<TProfile>(new TProfile("MultiParticleEfficiencyCutNorm","<# Tracks / # Charged Particles > vs # Charged Particles (all track cuts)",12,-0.5,11.5,"")); 
  MultiParticleEfficiencyCutNorm->SetDirectory(0); 
  MeanT2HitvsEta = std::auto_ptr<TProfile>(new TProfile("MeanT2HitvsEta","<# Hit> reconstructed vs generated particle #eta ",34,4.35,7.75,""));
  MeanT2HitvsEta->SetDirectory(0);

  MeanT2PadDigivsEta= std::auto_ptr<TProfile>(new TProfile("MeanT2PadDigivsEta","<# detector> with at least 1 Digi-Pad vs generated particle #eta ",34,4.35,7.75,""));
  MeanT2PadDigivsEta->SetDirectory(0);
 MeanT2StripDigivsEta= std::auto_ptr<TProfile>(new TProfile("MeanT2StripDigivsEta","<# detector> with at least 1 Digi-Strip vs generated particle #eta ",34,4.35,7.75,""));
  MeanT2StripDigivsEta->SetDirectory(0);

  MeanT2HitCl1vsEta = std::auto_ptr<TProfile>(new TProfile("MeanT2HitCl1vsEta","<# Hit> reconstructed vs generated particle #eta ",34,4.35,7.75,""));
  MeanT2HitCl1vsEta->SetDirectory(0);

  MeanT2GeantHitvsEta= std::auto_ptr<TProfile>(new TProfile("MeanT2GeantHitvsEta","<# detector> with at least 1 Geant Hit vs generated particle #eta ",34,4.35,7.75,""));
  MeanT2GeantHitvsEta->SetDirectory(0);

  // SingleParticleEfficiencyAllCutsH0 = std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiencyAllCutsH0","# Tracks reconstructed vs generated particle #eta (all cuts) HalfTele0",34,4.35,7.75,""));
  //SingleParticleEfficiencyAllCutsH0->SetDirectory(0);

  // SingleParticleEfficiencyAllCutsH1 = std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiencyAllCutsH1","# Tracks reconstructed vs generated particle #eta (all cuts) HalfTele1",34,4.35,7.75,""));
  //SingleParticleEfficiencyAllCutsH1->SetDirectory(0);

  MeanT2StripvsEta= std::auto_ptr<TProfile>(new TProfile("MeanT2StripvsEta","<# Strip> reconstructed vs generated particle #eta ",34,4.35,7.75,"")); 
  MeanT2PadvsEta= std::auto_ptr<TProfile>(new TProfile("MeanT2PadvsEta","<#  Pad> reconstructed vs generated particle #eta ",34,4.35,7.75,""));
  
  MeanT2RoadEntriesvsEta= std::auto_ptr<TProfile>(new TProfile("MeanT2RoadEntriesvsEta","<# Road Entries> reconstructed vs generated particle #eta ",34,4.35,7.75,""));
  
  MeanEtaResvsPhi= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhi","<Track #eta>  vs #phi (no cuts)",24,0.0,360.0,""));
  MeanEtaResvsPhi->SetDirectory(0); 
  MeanEtaResvsPhiHem1= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhiHem1","<Track #eta>  vs #phi (no cuts) hemis. 1",24,0.0,360.0,""));//15 deg
  MeanEtaResvsPhiHem1->SetDirectory(0); 
  MeanEtaResvsPhiHem2= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhiHem2","<Track #eta>  vs #phi (no cuts) hemis. 2",24,0.0,360.0,""));
  MeanEtaResvsPhiHem2->SetDirectory(0); 
  EtaParticle=std::auto_ptr<TH1F>(new TH1F("EtaParticle","Generated Particle #eta",1800,-9.0,9.0));
  EtaParticle->SetDirectory(0); 
  EtaParticleAll=std::auto_ptr<TH1F>(new TH1F("EtaParticleAll","Generated Particle #eta",60,-15.0,15.0));
  EtaParticleAll->SetDirectory(0);
  MeanEtaResvsPhiPart= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhiPart","<Generated #eta>  vs #phi ",24,0.0,360.0,""));
  MeanEtaResvsPhiPart->SetDirectory(0); 
  MeanEtaResvsPhiHem1Part= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhiHem1Part","<Generated #eta>  vs #phi  hemis. 1",24,0.0,360.0,""));//15 deg
  MeanEtaResvsPhiHem1Part->SetDirectory(0); 
  MeanEtaResvsPhiHem2Part= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhiHem2Part","<Generated #eta>  vs #phi  hemis. 2",24,0.0,360.0,""));
  MeanEtaResvsPhiHem2Part->SetDirectory(0); 
  MeanEtaResvsPhiCut= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhiCut","<Track #eta>  vs #phi (all cuts)",24,0.0,360.0,""));
  MeanEtaResvsPhiCut->SetDirectory(0);   
  MeanEtaResvsPhiHem1Cut= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhiHem1Cut","<Track #eta>  vs #phi (all cuts) hemis. 1",24,0.0,360.0,""));
  MeanEtaResvsPhiHem1Cut->SetDirectory(0); 
  MeanEtaResvsPhiHem2Cut= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhiHem2Cut","<Track #eta>  vs #phi (all cuts) hemis. 2",24,0.0,360.0,""));
  MeanEtaResvsPhiHem2Cut->SetDirectory(0); 

  MinDz2Plane= std::auto_ptr<TH1F>(new TH1F("MinDz2Plane","Minimum Z Distance resp. to second plane",400,0.0,400.0));
  MinDz2Plane->SetDirectory(0); 

  Trketachicut= std::auto_ptr<TH1F>(new TH1F("Trketachicut","Reconstructed Track #eta (#chi^{2}-p cut)",1800,-9.0,9.0));
 Trketachicut->SetDirectory(0);
  Trketaetacut = std::auto_ptr<TH1F>(new TH1F("Trketaetacut","Reconstructed Track #eta (#eta cut)",1800,-9.0,9.0));
  Trketaetacut->SetDirectory(0);
  
  Trketaetachicut= std::auto_ptr<TH1F>(new TH1F("Trketaetachicut","Reconstructed Track #eta (#eta , #chi^{2}-p cut)",1800,-9.0,9.0));
  Trketaetachicut->SetDirectory(0);
 
  Trketadzcut= std::auto_ptr<TH1F>(new TH1F("Trketadzcut","Reconstructed Track #eta (#Delta Z cut)",1800,-9.0,9.0));
  Trketadzcut->SetDirectory(0);
  Trketadzetacut= std::auto_ptr<TH1F>(new TH1F("Trketadzetacut","Reconstructed Track #eta (#eta, #Delta Z  cut)",1800,-9.0,9.0));
  Trketadzetacut->SetDirectory(0); 
  Trketadzetachicut= std::auto_ptr<TH1F>(new TH1F("Trketadzetachicut","Reconstructed Track #eta (#eta, #Delta Z  cut)",1800,-9.0,9.0));
  Trketadzetachicut->SetDirectory(0);
  

  NumChPartVsNumP0 = std::auto_ptr<TProfile>(new TProfile("NumChPartVsNumP0","# P0 in 4.5-7.5 vs # Ch Part in T2",21,-0.5,20.5,""));
  NumChPartInVsNumChPartOut= std::auto_ptr<TProfile>(new TProfile("NumChPartInVsNumChPartOut","# Ch outside T2  vs # Ch Part in T2",21,-0.5,20.5,""));

  NumP0VsNumTracks= std::auto_ptr<TProfile>(new TProfile("NumP0VsNumTracks","# P0 in 4.5-7.5  vs # Tracks (all cuts)",15,-0.5,14.5,""));

  

  char sZname2[1024];
 char sZnamehist[1024];

  for(unsigned int m=0;m<15; m++)
   {
       sprintf(sZname2, "PimenoEtaEnergy%d", m); 
       sprintf(sZnamehist, "PimenoEtaEnergy%d", m);
       PimenoEtaEnergy[m] = std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,20,5.0,7.0,40,1.,400.));
       PimenoEtaEnergy[m]->SetDirectory(0);
       PimenoEtaEnergy[m]->SetOption("lego");
   }

  chPartinT2->SetXTitle("# ch. Particle in T2");
  ChEnergyinT2->SetXTitle("ch. Particle Energy in T2");
  clusterstripentries->SetXTitle("# strips in one strip-cluster");
  clusterstripentries->SetYTitle("# cluster-strip");
  clusterpadentries->SetXTitle("# pad in one pad-cluster");
  clusterpadentries->SetYTitle("# cluster-pad");
  clusterpadentries->GetYaxis()->SetTitleOffset(1.25);
  clusterstripentries->GetYaxis()->SetTitleOffset(1.25);
  diffphiGUNHIT->SetXTitle("#Delta #Phi (deg)");
  RLocHRecoH->SetXTitle("#Delta R (mm)");
  diffRCluHit->SetXTitle("#Delta R (mm)");
  Trketa->SetXTitle("Track #eta");

  Trketachicut->SetXTitle("Track #eta");
  Trketaetacut->SetXTitle("Track #eta");
  Trketaetachicut->SetXTitle("Track #eta");
  Trketadzcut->SetXTitle("Track #eta");
  Trketadzetacut->SetXTitle("Track #eta");
  Trketadzetachicut->SetXTitle("Track #eta");


  Trkphi->SetXTitle("Track #phi (deg)");
  TrKPhi0degstudyg->SetXTitle("Track #phi (deg)");
  DEtaGoodTrk->SetXTitle("#Delta #eta");
  diffphiCluGun->SetXTitle("#Delta #Phi (deg)");
  diffphiCluHit->SetXTitle("#Delta #Phi (deg)");
  DPhiGoodTrk->SetXTitle("#Delta #Phi (deg)");
  Trkphigood->SetXTitle("Track #phi (deg)");
  Trketagood->SetXTitle("Track #eta");
  SingleParticleEfficiency->SetXTitle("Particle #eta");
  SingleParticleEfficiencyAllCuts->SetXTitle("Particle #eta");
  SingleParticleEfficiencyCutEta->SetXTitle("Particle #eta");
  SingleParticleEfficiencyCutChi->SetXTitle("Particle #eta");
  SingleParticleEfficiencyCutDZ ->SetXTitle("Particle #eta"); 
  SingleParticleEfficiency->SetYTitle("Pseudo-Efficiency");
  SingleParticleEfficiencyAllCuts->SetYTitle("Efficiency");
  MultiParticleEfficiencyCutNorm->SetYTitle("Efficiency");
  MultiParticleEfficiencyCutNorm->SetXTitle("ch. particle in T2");  
  MultiParticleEfficiencyCut->SetYTitle("<Track>");
  MultiParticleEfficiencyCut->SetYTitle("ch. particle in T2");    
  Chi2RProb->SetXTitle("Probability");
  Chi2PhiProb->SetXTitle("Probability");
  HTrackcounter->SetXTitle("# Tracks"); 
  HNumrecotrackcutmag0->SetXTitle("# Tracks Reconstructed");
  HTrackcounterNonCrit->SetXTitle("# Tracks Reconstructed");
  HTrackcounterCrit->SetXTitle("# Tracks Reconstructed");  
  stablepdg->SetXTitle("# Particle ID");  
  NumhitinTrack->SetXTitle("# Hit");  
  NumhitinTrackGood->SetXTitle("# Hit");  

  P0EtaEnergy->SetXTitle("#eta"); 
  P0EtaEnergy->SetYTitle("Energy (GeV)"); 
  // P0EtaEnergy->Draw("");
  P0EtaEnergy->SetOption("lego");
  ChOutEtaEnergy->SetOption("lego");

  muEtaEnergy->SetOption("lego");
  //Chi2RProbLogy->cd();
  // Chi2RProb->Draw("");
  //Chi2RProbLogy->SetLogy();
  //Chi2PhiProbLogy->cd();
  // Chi2PhiProb->Draw("");
  //Chi2PhiProbLogy->SetLogy();

  // Create a ROOT Tree
  tree->Branch("branchT2EV.",&t2obj);
}

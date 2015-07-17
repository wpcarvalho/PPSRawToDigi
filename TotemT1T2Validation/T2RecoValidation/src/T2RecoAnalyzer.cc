#include "TotemT1T2Validation/T2RecoValidation/interface/T2RecoValidation.h"

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/T2Digi/interface/T2StripDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2StripDigi.h"
#include "DataFormats/T2Digi/interface/T2PadDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2PadDigi.h"
#include "DataFormats/T2DetId/interface/T2DetId.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"

#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMath.h"

#include "CLHEP/Vector/LorentzVector.h"


T2RecoAnalyzer::T2RecoAnalyzer(const edm::ParameterSet& iConfig){

  SimTrackContainerLabel = iConfig.getParameter<edm::InputTag>("SimTrackContainerLabel");
  Chicut = iConfig.getParameter<double>("Chicut");
  PhiChiProbCut = iConfig.getParameter<double>("PhiChiProbCut");
  RChiProbCut = iConfig.getParameter<double>("RChiProbCut");  
  energy=iConfig.getParameter<double>("FitdzEnergy");
  DZScale=iConfig.getParameter<double>("DZScale");
  tracketamin=iConfig.getParameter<double>("TrkEtamin");
  tracketamax=iConfig.getParameter<double>("TrkEtaMAX");
  singleparticle=iConfig.getParameter<bool>("singleparticle");

  outputFileName = iConfig.getUntrackedParameter<std::string>("OutputFile");
  HepMCProductLabel = iConfig.getParameter<std::string>("HepMCProductLabel");
  CluLabel = iConfig.getParameter<std::string>("CluLabel");
  HitLabel = iConfig.getParameter<std::string>("HitLabel");
  RoadLabel = iConfig.getParameter<std::string>("RoadLabel");
  TrackLabel= iConfig.getParameter<std::string>("TrackLabel");
}

unsigned int T2RecoAnalyzer::RawtoSymb(uint32_t thedet)
{
  T2DetId converter;
  unsigned int pl=converter.plane(thedet);
  unsigned int pls=converter.planeSide(thedet);
  unsigned int ht=converter.halfTelescope(thedet);
  unsigned int arm=converter.arm(thedet);
  unsigned int symbolic=pl*2+pls+ht*10+20*arm;	  
 
  return symbolic;
}



T2RecoAnalyzer::~T2RecoAnalyzer()
{
}


double PartCharge(int pdgcode){
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
void T2RecoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace HepMC;
if((numevent%10)==0)
  std::cout<<numevent<<std::endl;
	numevent++;
  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  /* LOADING OF ALL THE RECORDS FROM THE EVENT */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */


  /* :::::::::::::Take The Clusters::::::::::::*/

  Handle<T2StripClusterCollection> t2strclcoll;
  iEvent.getByLabel(CluLabel,"T2StripClusters",t2strclcoll);
  Handle<T2PadClusterCollection> t2padclcoll;
  iEvent.getByLabel(CluLabel,"T2PadClusters",t2padclcoll);

  /*::::::Take particle generated with Gun::::::*/
  Handle<HepMCProduct> EvtHandle ;
  iEvent.getByLabel(HepMCProductLabel, EvtHandle ) ;
  const GenEvent* evt = EvtHandle->GetEvent();

  /*::::::Take Geant local hit on T2::::::*/
  /*
  Handle<CrossingFrame<PSimHit> > cFrame;
  iEvent.getByLabel("mix", "g4SimHitsTotemHitsT2Gem", cFrame);
  // get hits from G4Sim
  const string nameOfHits("TotemHitsT2Gem");
  auto_ptr<MixCollection<PSimHit> >  T2SimHits( new MixCollection<PSimHit>( cFrame.product() ) );

  // map hits to T2 plane
  map<int, PSimHitContainer> hitMap;
  for(MixCollection<PSimHit>::MixItr hitItr = T2SimHits->begin(); hitItr != T2SimHits->end(); ++hitItr) {
    hitMap[hitItr->detUnitId()].push_back(*hitItr);
  }
  */
  /*::::::Take  T2  Hits::::::*/
  Handle<T2HitCollection> t2hitcoll;

  iEvent.getByLabel(HitLabel,"T2Hits",t2hitcoll);
  /*::::::Take  T2  Roads::::::*/

  Handle<T2RoadCollection> t2roadcoll;
  //event.getByLabel(RoadModuleLabel,RoadInstanceLabel,myRoadColl);
  //iEvent.getByLabel(RoadLabel,"T2RoadColl",t2roadcoll);
iEvent.getByLabel("T2RoadPadFinder","NewRoadFinderRELOAD",t2roadcoll);
  /*:::::: Take T2 tracks ::::::*/
  Handle<T1T2TrackCollection> trackCollection;
  iEvent.getByLabel(TrackLabel,"T2TrackColl",trackCollection);

  //--***************************************************        Reco TRACK        ************************************************************--/
  //std::cout<<"Track reco Only"<<std::endl; 
 

  T1T2TrackCollection::const_iterator TrkCit;
  double trketa=0.;
  double trkphi=0.;
  unsigned int trackcounter=0;
  unsigned int trackcountergood=0;
  double DZgood;
  double numrecotrackcutmag0=0.;
  double chiRProb=0.;
  double chiPhiProb=0.;
  bool chi2condition;
  unsigned int trackplus=0;unsigned int trackminus=0;

  for(TrkCit=trackCollection->begin(); TrkCit!=trackCollection->end(); TrkCit++){
    trketa= (*TrkCit).Eta();

    if(trketa>0)
      trackplus++;
    else
      trackminus++;

     

    trkphi=(*TrkCit).Phi()*180/3.14159265;
    Trketa->Fill(trketa);
    Trkphi->Fill(trkphi);
    chiRProb=TMath::Prob((*TrkCit).ChiSquaredR(),((*TrkCit).GetHitEntries()-2));
    chiPhiProb=TMath::Prob((*TrkCit).ChiSquaredPhi(),((*TrkCit).GetHitEntries()-1));  
    Chi2RProb->Fill(chiRProb);
    Chi2PhiProb->Fill(chiPhiProb);
    
    if((chiPhiProb<PhiChiProbCut)&&(chiRProb<RChiProbCut))  
      chi2condition=false;
    else
      chi2condition=true;

    DZgood=1000.;
    if (energy==30.)
      DZgood=sigmaZfunctE30(fabs(trketa));
    if (energy==50.)
      DZgood=sigmaZfunctE50(fabs(trketa));

    if((chi2condition)&&(IsinT2(trketa,tracketamin,tracketamax))&&(fabs((*TrkCit).Z_at_Rmin())<DZgood)){
      trackcountergood++;
      Trketagood->Fill(trketa);
      Trkphigood->Fill(trkphi);
      if(trketa>0)	       
	numrecotrackcutmag0= numrecotrackcutmag0+1;
      
      
    }
    trackcounter++;
  }

  double lastParticleeta=0.;
  double chp12mag0=0.;
  double chp12min0=0.;
  int num_somethinginplus=0;int num_somethinginminus=0;
  int num_somethinginPlusBelowT2=0; int num_somethinginPlusAboveT2=0;
  int num_somethinginplusCH=0;int num_somethinginminusCH=0;
  int somethingingapPlus=0; int somethingingapMinus=0; int allstableCharged=0;
  double ET2Plus=0.;double ET2Minus=0.;

  double forwInvMass=0.; double CentInvMass=0.;  double PomeronInvMass=0.;double AllRecoMass=0.;
  std::vector<int> protRP;
  double EFw=0.; double PxFw=0.; double PyFw=0.; double PzFw=0.;
  double ECent=0.; double PxCent=0.; double PyCent=0.; double PzCent=0.;
  double EProt=0.; double PxProt=0.; double PyProt=0.; double PzProt=0.;
  double EAll=0.; double PxAll=0.; double PyAll=0.; double PzAll=0.;
  double csi1=0.; double csi2=0.;double rap1=0.;
  double EAllNoProton=0.; double PxAllNoProton=0.; double PyAllNoProton=0.; double PzAllNoProton=0.;


  //Additional Loop xp= 1-E'/E. 
  int centralCMS_ch=0; 
  int numt1_ch=0; 
  int numfsc_ch=0;
  int numt2trk=trackplus+trackminus;
  int thebarc=0;
  for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
    lastParticleeta=(*p)->momentum().eta();
   
    
    
    
    if((*p)->pdg_id()==2212)
       if((*p)->status()==1)
	 {
	   //                                std::cout<<(*p)->momentum().e()<<" "<<((*p)->momentum().eta())<<std::endl;
	   if(fabs((*p)->momentum().eta())>8.){
	     
	     if((*p)->momentum().e()>2000.){
	       thebarc=(*p)->barcode();
	       protRP.push_back(thebarc);
	     }
	   }
	 }


    if(((*p)->status()==1)&&((*p)->momentum().eta()<6.5)&&((*p)->momentum().eta()>5.3))
      ET2Plus= ET2Plus+(*p)->momentum().e();
      
    
     if(((*p)->status()==1)&&((*p)->momentum().eta()>-6.5)&&((*p)->momentum().eta()<-5.3))
      ET2Minus= ET2Minus+(*p)->momentum().e();
      

    if ((PartCharge((*p)->pdg_id())!=0)&&((*p)->status()==1))
      allstableCharged++;

    if ((PartCharge((*p)->pdg_id())!=0)&&((*p)->status()==1)&&((*p)->momentum().eta()<5.3)&&((*p)->momentum().eta()>4.7))
      somethingingapPlus++;

    if ((PartCharge((*p)->pdg_id())!=0)&&((*p)->status()==1)&&((*p)->momentum().eta()>-5.3)&&((*p)->momentum().eta()<-4.7))
      somethingingapMinus++;


    if ((PartCharge((*p)->pdg_id())!=0)&&((*p)->status()==1)&&((*p)->momentum().eta()<6.5)&&((*p)->momentum().eta()>5.3))
      chp12mag0++;
    
    if ((PartCharge((*p)->pdg_id())!=0)&&((*p)->status()==1)&&((*p)->momentum().eta()<-5.3)&&((*p)->momentum().eta()>-6.5))
      chp12min0++;

    if(((*p)->status()==1)&&((*p)->momentum().eta()<6.5)&&((*p)->momentum().eta()>5.3)){
      num_somethinginplus++;
      if(PartCharge((*p)->pdg_id())!=0)
	 num_somethinginplusCH++;
    }

    if(((*p)->status()==1)&&((*p)->momentum().eta()>-6.5)&&((*p)->momentum().eta()<-5.3)){
      num_somethinginminus++;
      if(PartCharge((*p)->pdg_id())!=0)
	 num_somethinginminusCH++;

    }
 
    if(((*p)->status()==1)&&((*p)->momentum().eta()>6.5))
      num_somethinginPlusBelowT2++;

    if(((*p)->status()==1)&&((*p)->momentum().eta()<5.3)&&((*p)->momentum().eta()>0))
      num_somethinginPlusAboveT2++;

  }


  
  for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
    CumulativeEtaGeneratorALLParticles->Fill((*p)->momentum().eta());
    if(PartCharge((*p)->pdg_id())!=0)
      if((*p)->status()==1)
	CumulativeEtaGeneratorStableCHParticles->Fill((*p)->momentum().eta());
  
    if((*p)->status()==1)
      CumulativeEtaGeneratorStableParticles->Fill((*p)->momentum().eta());
    
    if((chp12mag0+chp12min0)>0)
      if((*p)->status()==1)
	CumulativeEtaGeneratorStableParticles_WhenACHARGEInT2->Fill((*p)->momentum().eta());
      
    
 }

  
  //std::cout<<protRP.size()<<std::endl;
 if(protRP.size()==2)
    {
      NumEvtSel2Proton->Fill(1.);

      for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
	
	if((*p)->status()==1){
	  EAll=EAll+(*p)->momentum().e();
	  PxAll=PxAll+ (*p)->momentum().px(); 
	  PyAll=PyAll+ (*p)->momentum().py(); 
	  PzAll=PzAll+ (*p)->momentum().pz(); 
	  if(((*p)->barcode())!=protRP.at(0))
	    if(((*p)->barcode())!=protRP.at(1)){
	        EAllNoProton=EAllNoProton+(*p)->momentum().e();
		PxAllNoProton=PxAllNoProton+ (*p)->momentum().px(); 
		PyAllNoProton=PyAllNoProton+ (*p)->momentum().py(); 
		PzAllNoProton=PzAllNoProton+ (*p)->momentum().pz();
	    }
	}
	
	
	if(((*p)->barcode())==protRP.at(0)) 
	  if((*p)->status()==1){
	    EProt=EProt+(*p)->momentum().e();
	    PxProt=PxProt+ (*p)->momentum().px(); 
	    PyProt=PyProt+ (*p)->momentum().py(); 
	    PzProt=PzProt+ (*p)->momentum().pz(); 
	    csi1=(4000.-(*p)->momentum().e())/4000.;
	    rap1=0.5*log((EProt+PzProt)/(EProt-PzProt));

	  }
	
	if(((*p)->barcode())==protRP.at(1))
	  if((*p)->status()==1){
	    EProt=EProt+(*p)->momentum().e();
	    PxProt=PxProt+ (*p)->momentum().px(); 
	    PyProt=PyProt+ (*p)->momentum().py(); 
	    PzProt=PzProt+ (*p)->momentum().pz();
	    csi2=(4000.-(*p)->momentum().e())/4000.;
	   // rap2=0.5*log((EProt+PzProt)/(EProt-PzProt));
	  }


	if((*p)->status()==1) 
	  if(PartCharge((*p)->pdg_id())!=0)
	    if(fabs((*p)->momentum().eta())>6)
	      if(fabs((*p)->momentum().eta())<8)
		numfsc_ch++;


	if((*p)->status()==1) 
	  if(PartCharge((*p)->pdg_id())!=0)
	    if(fabs((*p)->momentum().eta())>3.1)
	      if(fabs((*p)->momentum().eta())<4.7)
		numt1_ch++;
		
		

	      
	if(fabs((*p)->momentum().eta())<5.5)  //Default 5.5
	  if((*p)->status()==1){
	    ECent=ECent+(*p)->momentum().e();
	    PxCent=PxCent+ (*p)->momentum().px(); 
	    PyCent=PyCent+ (*p)->momentum().py(); 
	    PzCent=PzCent+ (*p)->momentum().pz(); 
	    if(PartCharge((*p)->pdg_id())!=0)
	      if((*p)->momentum().e()>5.)
		centralCMS_ch++;
	  }
	
	if(fabs((*p)->momentum().eta())>5.3)
	  if((*p)->status()==1)
	    if(((*p)->barcode())!=protRP.at(0))
	      if(((*p)->barcode())!=protRP.at(1)) 
		{
		  EFw=EFw+(*p)->momentum().e();
		  PxFw=PxFw+ (*p)->momentum().px(); 
		  PyFw=PyFw+ (*p)->momentum().py(); 
		  PzFw=PzFw+ (*p)->momentum().pz(); 	     
		}	
      }
    }
 
 double AllMassNoProton=0.;
 forwInvMass=0.;  
 CentInvMass=0.;   
 PomeronInvMass=0.;
 AllRecoMass=0.;
 //csi2Histo->Fill(csi2);
 // std::cout<<EAll<<std::endl;

 AllMassNoProton=sqrt(EAllNoProton*EAllNoProton-PxAllNoProton*PxAllNoProton-PyAllNoProton*PyAllNoProton-PzAllNoProton*PzAllNoProton);
 AllRecoMass=sqrt(EAll*EAll-PxAll*PxAll-PyAll*PyAll-PzAll*PzAll);
 forwInvMass=sqrt(EFw*EFw-PxFw*PxFw-PyFw*PyFw-PzFw*PzFw);
 CentInvMass=sqrt(ECent*ECent-PxCent*PxCent-PyCent*PyCent-PzCent*PzCent);
 PomeronInvMass=sqrt(csi1*csi2)*8000.0;//sqrt(EProt*EProt-PxProt*PxProt-PyProt*PyProt-PzProt*PzProt);

 DifferenceMassFromProtons_andAllOtherParticles->Fill(PomeronInvMass-AllMassNoProton);
 csi1Histo->Fill(csi1);




 
 if((csi1>0.05)&&(csi2>0.05)){
   CentralEnergyWhenHigh005Csi->Fill(ECent);
   if(ECent<100){
     T1chargeWhen_CentralEnergyWhenHigh005Csi_100->Fill(numt1_ch);
     FSCchargeWhen_CentralEnergyWhenHigh005Csi_100->Fill(numfsc_ch);  
     T2TrkWhen_CentralEnergyWhenHigh005Csi_100->Fill(numt2trk); 
     ForwardWhen_CentralEnergyWhenHigh005Csi_100->Fill(numt1_ch+numfsc_ch+numt2trk);
   }
   
   if(ECent<30){
     T1chargeWhen_CentralEnergyWhenHigh005Csi_30->Fill(numt1_ch);
     FSCchargeWhen_CentralEnergyWhenHigh005Csi_30->Fill(numfsc_ch);  
     T2TrkWhen_CentralEnergyWhenHigh005Csi_30->Fill(numt2trk);
     ForwardWhen_CentralEnergyWhenHigh005Csi_30->Fill(numt1_ch+numfsc_ch+numt2trk);
   }

   if((numt1_ch+numfsc_ch+numt2trk)==0){
      CentralEnergyWhenHigh005Csi_AndfscT1T2OFF->Fill(ECent);
      
   }

 }
  



 csi1_vs_csi2->Fill(csi1,csi2);

 ProtonRapidity->Fill(rap1);
 EtaMaxAroundLeftProton->Fill((rap1 -log(csi1)));


 if((numt1_ch+numfsc_ch+numt2trk)==0){
   CentralEnergyWhenfscT1T2OFF->Fill(ECent);
   PomeronInvMassWhenfscT1T2OFFvsCentralE->Fill(PomeronInvMass,ECent);
 }


 /*
 NumEvtSel2Proton->Fill(1.);
 AllRecoMassHist->Fill();
 forwInvMassHist->Fill();
 CentInvMassHist->Fill();
 PomeronInvMassHist->Fill();
 */

  AllRecoMassHist->Fill(AllRecoMass);
  forwInvMassHist->Fill(forwInvMass);
  CentInvMassHist->Fill(CentInvMass);
  PomeronInvMassHist->Fill(PomeronInvMass);
  double difference=(PomeronInvMass-CentInvMass);
  NonCentralMassFromPom->Fill(difference);






 bool cutcsi_003_008=false;


 if((csi1<0.08)&&(csi1>0.03))
   if((csi2<0.08)&&(csi2>0.03))
     cutcsi_003_008=true;



 
  if(cutcsi_003_008)
    {

      MultiplicityWhenT2ShouldHAVE->Fill((trackplus+trackminus));
      Pom_minus_Cent_003_008->Fill(PomeronInvMass-CentInvMass);
      
      if((numt1_ch+numfsc_ch+numt2trk)==0){
	
	if((PomeronInvMass-CentInvMass)>300.)
	  ForwEnergyWhenForwOFF_003_008_Cut_Discr300GeV->Fill(EFw);//forwInvMass
	
	if((PomeronInvMass-CentInvMass)>400.)
	  ForwEnergyWhenForwOFF_003_008_Cut_Discr400GeV->Fill(EFw);
	
	if((PomeronInvMass-CentInvMass)>100.)
	  ForwEnergyWhenForwOFF_003_008_Cut_Discr100GeV->Fill(EFw);
	
	if((PomeronInvMass-CentInvMass)>200.)
	  ForwEnergyWhenForwOFF_003_008_Cut_Discr200GeV->Fill(EFw);
      }      


       if((numt1_ch+numt2trk)==0){
	
	if((PomeronInvMass-CentInvMass)>300.)
	  ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC300GeV->Fill(EFw);//forwInvMass ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC100GeV
	
	if((PomeronInvMass-CentInvMass)>400.)
	  ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC400GeV->Fill(EFw);
	
	if((PomeronInvMass-CentInvMass)>100.)
	  ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC100GeV->Fill(EFw);
	
	if((PomeronInvMass-CentInvMass)>200.)
	  ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC200GeV->Fill(EFw);
      }      




    }
 




  bool analysis_consistent=false;

 
 if((csi1>0.05))
   if((csi2>0.05))
     analysis_consistent=true;
   


if(analysis_consistent){
  //MultiplicityWhenT2ShouldHAVE->Fill((trackplus+trackminus));

  if((numt1_ch+numfsc_ch+numt2trk)==0)
    CentralEnergyWhenfscT1T2OFF_An005Cut->Fill(ECent);

  //((trackplus+trackminus)==0)

  if((numt1_ch+numfsc_ch+numt2trk)==0){
   
    if((PomeronInvMass-CentInvMass)>100.)
      ForwEnergyWhenForwOFF_005Cut_Discr100GeV->Fill(EFw);

     if((PomeronInvMass-CentInvMass)>200.)
      ForwEnergyWhenForwOFF_005Cut_Discr200GeV->Fill(EFw);

     if((PomeronInvMass-CentInvMass)>300.)
      ForwEnergyWhenForwOFF_005Cut_Discr300GeV->Fill(EFw);//forwInvMass

    if((PomeronInvMass-CentInvMass)>400.)
      ForwEnergyWhenForwOFF_005Cut_Discr400GeV->Fill(EFw);

    bool cmsOFF=false;
    if((ECent>20)||(centralCMS_ch>0))
      cmsOFF=true;

    if(cmsOFF){
     if((PomeronInvMass-CentInvMass)>100.)
      ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr100GeV->Fill(EFw);

     if((PomeronInvMass-CentInvMass)>200.)
      ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr200GeV->Fill(EFw);

     if((PomeronInvMass-CentInvMass)>300.)
      ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr300GeV->Fill(EFw);//forwInvMass

    if((PomeronInvMass-CentInvMass)>400.)
      ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr400GeV->Fill(EFw);
    }

  }
 }



if(trackplus==0)
    NonTriggeredEnergyPlus->Fill(ET2Plus);

  if((trackplus+trackminus)==0)
    NonTriggeredEnergy->Fill(ET2Plus+ET2Minus);


  if((PomeronInvMass-CentInvMass)>100.){
    
    if(trackplus==0)
       NonTriggeredEnergyInPlusOnly_When100GeVDiscr->Fill(ET2Plus);
    
    if((trackplus+trackminus)==0)
      NonTriggeredEnergyInPlusAndMinus_When100GeVDiscr->Fill(ET2Minus+ET2Plus);
  }














  bool leakaget1t2=false;
  if(somethingingapPlus>0)
    if(somethingingapMinus>0)
      if((somethingingapMinus+somethingingapPlus)==allstableCharged)
	leakaget1t2=true;
  
  if(leakaget1t2)
    LeakT1T2->Fill(1);
  else
    LeakT1T2->Fill(0);


  if(num_somethinginPlusAboveT2==0)
    if(num_somethinginPlusBelowT2>0) 
      if(num_somethinginplus==0){
	if(trackplus>0){
	  ProbGetLowMass->Fill(1);
	  std::cout<<"******************Extra Low mass detected:"<<std::endl;
	  for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
	    if((*p)->status()==1)
	      if((*p)->momentum().eta()>0)
		std::cout<<"Eta-ID:"<<(*p)->momentum().eta()<<" "<<(*p)->pdg_id()<<std::endl;
	  }
	  
	}
	else
	  ProbGetLowMass->Fill(0);
      }

  if((num_somethinginminusCH>0)&&(num_somethinginplusCH==0))
    GeneratorEvtHemispCH->Fill(-1);
  if((num_somethinginminusCH==0)&&(num_somethinginplusCH>0))
    GeneratorEvtHemispCH->Fill(1);
  if((num_somethinginminusCH==0)&&(num_somethinginplusCH==0))
    GeneratorEvtHemispCH->Fill(0);
  if((num_somethinginminusCH>0)&&(num_somethinginplusCH>0))
    GeneratorEvtHemispCH->Fill(2);

  
  if((num_somethinginminus>0)&&(num_somethinginplus==0))
    GeneratorEvtHemisp->Fill(-1);
  if((num_somethinginminus==0)&&(num_somethinginplus>0))
    GeneratorEvtHemisp->Fill(1);
  if((num_somethinginminus==0)&&(num_somethinginplus==0))
    GeneratorEvtHemisp->Fill(0);
  if((num_somethinginminus>0)&&(num_somethinginplus>0))
    GeneratorEvtHemisp->Fill(2);
  
  for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
    if((*p)->status()==1){

      if((num_somethinginminus>0)&&(num_somethinginplus==0)){
	CumulativeEtaGeneratorStableParticlesLeftOnly->Fill((*p)->momentum().eta());
      }
      
      if((num_somethinginminus==0)&&(num_somethinginplus>0)){	
	CumulativeEtaGeneratorStableParticlesRightOnly->Fill((*p)->momentum().eta());
      }

      
    } 
  }
  


  if(num_somethinginplus==0){
    if(trackplus>0)
      TrackPlusWhenGenRapGapPlus->Fill(1);
    else
      TrackPlusWhenGenRapGapPlus->Fill(0);
  }
 
  T2RoadCollection::const_iterator RoadCit;
  unsigned int   roadplus=0;
  unsigned int   roadminus=0;
  unsigned int   roadH0=0;
  unsigned int   roadH1=0;
  unsigned int   roadH2=0;
  unsigned int   roadH3=0;
  T2GeometryUtil convert;
  for(RoadCit=t2roadcoll->begin();RoadCit!=t2roadcoll->end();RoadCit++)
    {
      T2Cluster roadclu=(*RoadCit).thisPadRoad.at(0);
      T2GeometryUtil::T2DetInfo planeinfo3;
      // double Z=roadclu.GetClusterZ();
      T2ROGeometry t2rogeo(roadclu.GetDetID());
      //std::cout<<"HEREb2: "<<std::endl;
      planeinfo3=convert.GetT2Info(roadclu.GetDetID());
      if(planeinfo3.symb<20){
	roadplus++;
	if(planeinfo3.symb<10)
	  roadH0++;
	else
	  roadH1++;
      }
      else{
	roadminus++;
	if(planeinfo3.symb<30)
	  roadH2++;
	else
	  roadH3++;

      }
    }


  NumTotRoadH0->Fill(roadH0);
  NumTotRoadH1->Fill(roadH1);
  NumTotRoadH2->Fill(roadH2);
  NumTotRoadH3->Fill(roadH3);

   if((roadplus==0)&&(roadminus>0))
    RoadHemisphereNorequirement->Fill(-1);
      
  if((roadplus>0)&&(roadminus==0))
    RoadHemisphereNorequirement->Fill(1);
  
  if((roadplus>0)&&(roadminus>0))
    RoadHemisphereNorequirement->Fill(2);
  
  if((roadplus==0)&&(roadminus==0))
    RoadHemisphereNorequirement->Fill(0);
      
  //  std::cout<<"Roads+-: "<<roadplus<<" "<<roadminus<<std::endl;






  if((trackplus==0)&&(trackminus>0))
    TrackHemisphereNorequirement->Fill(-1);
      
  if((trackplus>0)&&(trackminus==0))
    TrackHemisphereNorequirement->Fill(1);
  
  if((trackplus>0)&&(trackminus>0))
    TrackHemisphereNorequirement->Fill(2);
  
  if((trackplus==0)&&(trackminus==0))
    TrackHemisphereNorequirement->Fill(0);
    
  double energyEfficiency=0.;
  if(trackplus>0)
    energyEfficiency=1.;
  else
    energyEfficiency=0.;
  

  EfficiencyTrigT2PlusVSEplus->Fill(ET2Plus,energyEfficiency);
  
  
  if((chp12min0+chp12mag0)>0)
    {
      if((trackplus==0)&&(trackminus>0))
	TrackHemisphereWhenChargedParticle->Fill(-1);
      
      if((trackplus>0)&&(trackminus==0))
	TrackHemisphereWhenChargedParticle->Fill(1);

       if((trackplus>0)&&(trackminus>0))
	TrackHemisphereWhenChargedParticle->Fill(2);

       if((trackplus==0)&&(trackminus==0))
	TrackHemisphereWhenChargedParticle->Fill(0);
      
    }

  EventEffiXChpartInT2YNumTrk->Fill((chp12mag0+chp12min0),trackcounter);

  NumRecoTrkPerEvent->Fill(trackcounter);

  if((chp12mag0+chp12min0)>0)
    NumRecoTrkWhenAtLeastAstableInT2->Fill(trackcounter);

  SingleParticleEfficiency->Fill(lastParticleeta,trackcounter);
  SingleParticleEfficiencyCut->Fill(lastParticleeta,trackcountergood);


  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  GEANT TRACK @@@@@@@@@@@@@@@@@@@@@@@@@@
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  GEANT TRACK @@@@@@@@@@@@@@@@@@@@@@@@@@
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  GEANT TRACK @@@@@@@@@@@@@@@@@@@@@@@@@@

  double SingleTrackEvent_EtaRec=0;
  double PI = 3.1415927;
  double SingleTrackEvent_PhiRec=0;
  

  //get G4 tracks
  
//std::cout<<"Track Reco + Geant"<<std::endl;
  Handle<SimTrackContainer> G4TrkContainer;
  iEvent.getByLabel(SimTrackContainerLabel, G4TrkContainer);
  if (!G4TrkContainer.isValid()) {
    LogError("TrackerHitAnalyzer::analyze") << "Unable to find SimTrack in event!";
    return;
  }

  if((trackcountergood==1)&&(singleparticle)) 
  {
    vector< pair<double,double> > PrimarySimTracks;
    for (SimTrackContainer::const_iterator itTrk = G4TrkContainer->begin(); itTrk != G4TrkContainer->end(); ++itTrk) {
      double eta =0, phi =0, p =0;
      const CLHEP::HepLorentzVector G4Trk(itTrk->momentum().x(),itTrk->momentum().y(),itTrk->momentum().z(), itTrk->momentum().e() ) ;
      p = sqrt(G4Trk[0]*G4Trk[0]+G4Trk[1]*G4Trk[1]+G4Trk[2]*G4Trk[2]);
      if ( p == 0){
	LogError("TrackerHitAnalyzer::analyze") << "TrackerTest::INFO: Primary has p = 0 ";
      } else {
	double costheta  = G4Trk[2]/p;
	double theta = acos(TMath::Min(TMath::Max(costheta, -1.),1.));
	eta = -log(tan(theta/2));
	if ( G4Trk[0] != 0 || G4Trk[1] != 0)
	  phi = atan2(G4Trk[1],G4Trk[0]);

	if(phi<0)
	  phi = 2*PI + phi;
	phi=phi*180/PI;

	pair<double,double> *myprimarytrack = new pair<double,double>(eta,phi);
	PrimarySimTracks.push_back(*myprimarytrack);
      }
    }

  // comparison reconstructed track  - GEANT track

    for(T1T2TrackCollection::const_iterator TC_it=trackCollection->begin(); TC_it!=trackCollection->end(); TC_it++)
      {
	DZgood=1000.;
	if (energy==30.){
	  DZgood=sigmaZfunctE30(fabs((*TC_it).Eta()));
	}
	if (energy==50.) {
	  DZgood=sigmaZfunctE50(fabs((*TC_it).Eta()));
	}

	SingleTrackEvent_EtaRec = (*TC_it).Eta();
	SingleTrackEvent_PhiRec = (*TC_it).Phi()*180.0/3.14159265;

	for(vector< pair<double,double> >::iterator pair_it=PrimarySimTracks.begin(); pair_it!=PrimarySimTracks.end(); pair_it++){
	  double DE_temp = SingleTrackEvent_EtaRec - (*pair_it).first;
	  double DF_temp = SingleTrackEvent_PhiRec - (*pair_it).second;

	  ////std::cout << "Track DE, DF = " << DE_temp<< " , " << DF_temp << std::endl;
	  ////std::cout << "Phi reco: "<< SingleTrackEvent_PhiRec << "Phi G4: "  << (*pair_it).second << std::endl;
	  chiRProb=TMath::Prob((*TC_it).ChiSquaredR(),((*TC_it).GetHitEntries()-2));
	  chiPhiProb=TMath::Prob((*TC_it).ChiSquaredPhi(),((*TC_it).GetHitEntries()-1));          
	  if((chiPhiProb<PhiChiProbCut)&&(chiRProb<RChiProbCut))  
	    chi2condition=false;
	  else
	    chi2condition=true;

	  if((chi2condition)&&(IsinT2((*TC_it).Eta(),tracketamin,tracketamax))){
	    double DZgood=1000.;

	    if (energy==30.)
	      DZgood=sigmaZfunctE30(fabs((*TC_it).Eta()));
	    if (energy==50.)
	      DZgood=sigmaZfunctE50(fabs((*TC_it).Eta()));
	    
	    if(fabs((*TC_it).Z_at_Rmin())<DZgood){
	      DEtaGoodTrk->Fill(DE_temp);
	      DPhiGoodTrk->Fill(DF_temp);
	    }

	  }
	}
      }
  }

  //::::::::::Local variable Declaration:::::::::





  T2DetId myT2Det;






  //int zmm[2][10]; //2 row for halftelescope. 10 coloumn for planes
       //14035.605; //first Gem first drift gas zone (mm)

  //25.0;

        //4.5; From first Drift zone to RO board = 9mm


  // **************************************************** STRIP STUDIES ************************************************ //
  // **************************************************** STRIP STUDIES ************************************************ //
  // **************************************************** STRIP STUDIES ************************************************ //
  // **************************************************** STRIP STUDIES ************************************************ //
  // **************************************************** STRIP STUDIES ************************************************ //
//std::cout<<"strip Only"<<std::endl;

  for(T2StripClusterCollection::const_iterator itstrip = t2strclcoll->begin(); itstrip != t2strclcoll->end(); itstrip++){
    vector<T2Cluster> stripClv = itstrip->second;
    for(unsigned int k=0;k<stripClv.size();k++){
      clusterstripentries->Fill(stripClv[k].GetNoOfEntries());
    }
  }


//std::cout<<"strip+geant"<<std::endl;
/*
 if((trackcountergood==1)&&(singleparticle))  
 for(T2StripClusterCollection::const_iterator itstrip = t2strclcoll->begin(); itstrip != t2strclcoll->end(); itstrip++){
      myT2Det = itstrip->first;
      myhalftele= myT2Det.halfTelescope();
      myplaneside= myT2Det.planeSide();
      myplane= myT2Det.plane();
      vector<T2Cluster> stripClv = itstrip->second;
      for(unsigned int k=0;k<stripClv.size();k++){
    

      for(map<int, PSimHitContainer>::const_iterator hitMapItr = hitMap.begin(); hitMapItr != hitMap.end(); ++hitMapItr){

          const PSimHitContainer & planeSimHits = hitMapItr->second;
          T2DetId *theT2DetId =new T2DetId(hitMapItr->first);

          if((theT2DetId->plane()==myplane)&&(theT2DetId->planeSide()==myplaneside)&&(planeSimHits.size()==1)&&(stripClv.size()==1)){
            for (unsigned int l=0; l<planeSimHits.size(); l++){

              PSimHit myplanehits= planeSimHits[l];
              lmyx= myplanehits.localPosition().x();
              lmyy= myplanehits.localPosition().y();
              lmyr= sqrt(lmyx*lmyx + lmyy*lmyy);
	      diffRCluHit->Fill(stripClv[k].GetClusterR()-lmyr);
	      
	    }
	  }
      }

    }
  }
*/


  // **************************************************** PAD STUDIES ************************************************ //
  // **************************************************** PAD STUDIES ************************************************ //
  // **************************************************** PAD STUDIES ************************************************ //
  // **************************************************** PAD STUDIES ************************************************ //
  // **************************************************** PAD STUDIES ************************************************ //
//std::cout<<"pad+geant"<<std::endl;



 int ActivePlH0[10]= {0,0,0,0,0,0,0,0,0,0};
 int ActivePlH1[10]= {0,0,0,0,0,0,0,0,0,0};
 int ActivePlH2[10]= {0,0,0,0,0,0,0,0,0,0}; 
 int ActivePlH3[10]= {0,0,0,0,0,0,0,0,0,0};

  unsigned int intplane=0;
  unsigned int totHitAllH0=0;unsigned int totHitAllH1=0;unsigned int totHitAllH2=0;unsigned int totHitAllH3=0;
  
  for(T2PadClusterCollection::const_iterator itpad= t2padclcoll->begin(); itpad != t2padclcoll->end(); itpad++){//Reft2padclcoll
    vector<T2Cluster> padClv = itpad->second;
    
    for(unsigned int k=0;k<padClv.size();k++)
      clusterpadentries->Fill(padClv[k].GetNoOfEntries());
    
    T2DetId *detID =new T2DetId(itpad->first);
    
    uint32_t cmsswdId= detID->calculateRawId(detID->arm(),detID->halfTelescope(),detID->plane(),detID->planeSide());
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

if((trackplus>0)&&(trackminus>0))
  for(unsigned int i=0;i<10;i++){
    PadCluMultiplVsPlane->Fill((i+30),ActivePlH3[i]);
    PadCluMultiplVsPlane->Fill((i+20),ActivePlH2[i]);
    PadCluMultiplVsPlane->Fill((i+10),ActivePlH1[i]);
    PadCluMultiplVsPlane->Fill((i),ActivePlH0[i]);    
  }


 int ActivePlH0str[10]= {0,0,0,0,0,0,0,0,0,0};
 int ActivePlH1str[10]= {0,0,0,0,0,0,0,0,0,0};
 int ActivePlH2str[10]= {0,0,0,0,0,0,0,0,0,0}; 
 int ActivePlH3str[10]= {0,0,0,0,0,0,0,0,0,0};

 totHitAllH0=0;totHitAllH1=0; totHitAllH2=0;totHitAllH3=0;
 intplane=0;

 for(T2StripClusterCollection::const_iterator itstrip = t2strclcoll->begin(); itstrip != t2strclcoll->end(); itstrip++){
    vector<T2Cluster> stripClv = itstrip->second;
    
    
    T2DetId *detID =new T2DetId(itstrip->first);
    
    uint32_t cmsswdId= detID->calculateRawId(detID->arm(),detID->halfTelescope(),detID->plane(),detID->planeSide());
    unsigned int symb=RawtoSymb(cmsswdId);
   
    if((symb/10)==0){
	totHitAllH0=totHitAllH0+stripClv.size();
	intplane=(symb%10);
	ActivePlH0str[intplane]=ActivePlH0str[intplane]+stripClv.size();
      }
      if((symb/10)==1){	
	totHitAllH1=totHitAllH1+stripClv.size();
	intplane=(symb%10);
	ActivePlH1str[intplane]=ActivePlH1str[intplane]+stripClv.size();

      }
      if((symb/10)==2){
	totHitAllH2=totHitAllH2+stripClv.size();
	intplane=(symb%10);
	ActivePlH2str[intplane]=ActivePlH2str[intplane]+stripClv.size();
      }
      if((symb/10)==3){
	totHitAllH3=totHitAllH3+stripClv.size();     
	intplane=(symb%10);
	ActivePlH3str[intplane]=ActivePlH3str[intplane]+stripClv.size();
      }
 }

for(unsigned int i=0;i<10;i++){
    StripCluMultiplVsPlane->Fill((i),ActivePlH0str[i]);
    StripCluMultiplVsPlane->Fill((i+10),ActivePlH1str[i]);
    StripCluMultiplVsPlane->Fill((i+20),ActivePlH2str[i]);
    StripCluMultiplVsPlane->Fill((i+30),ActivePlH3str[i]);    
  }

/*
  if((trackcountergood==1)&&(singleparticle))                       //Comupte the differences only if there is no secondaries
    for(T2PadClusterCollection::const_iterator itpad = t2padclcoll->begin(); itpad != t2padclcoll->end(); itpad++){

      myT2Det = itpad->first;
      myhalftele= myT2Det.halfTelescope();
      myplaneside= myT2Det.planeSide();
      myplane= myT2Det.plane();

      vector<T2Cluster> padClv = itpad->second;

      //int symbplaneName=
      //padClv.size()

      for(unsigned int k=0;k<padClv.size();k++){

	// clusterpadentries->Fill(padClv[k].GetNoOfEntries());

        // now we run for each plane the digi-simulation
        for(map<int, PSimHitContainer>::const_iterator hitMapItr = hitMap.begin(); hitMapItr != hitMap.end(); ++hitMapItr){

          const PSimHitContainer & planeSimHits = hitMapItr->second;
          T2DetId *theT2DetId =new T2DetId(hitMapItr->first);

          if((theT2DetId->plane()==myplane)&&(theT2DetId->planeSide()==myplaneside)&&(planeSimHits.size()==1)&&(padClv.size()==1)){
            for (unsigned int l=0; l<planeSimHits.size(); l++){

              PSimHit myplanehits= planeSimHits[l];
              lmyx= myplanehits.localPosition().x();
              lmyy= myplanehits.localPosition().y();
              lmyr= sqrt(lmyx*lmyx + lmyy*lmyy);

              //strange GEANT x-y sign convention
              if((myplaneside==1)&&(myhalftele==0)){
                phipart=atan2(lmyy,lmyx);
                phipart=phipart*180.0/3.14159265;
                if (phipart<0)
                  phipart= 360.0 - fabs(phipart);
              }
              if((myplaneside==0)&&(myhalftele==0)){
                phipart=atan2(-lmyy,lmyx);
                phipart=phipart*180.0/3.14159265;
                if (phipart<0)
                  phipart= 360.0 - fabs(phipart);
              }
              if((myplaneside==1)&&(myhalftele==1)){
                phipart=atan2(-lmyy,-lmyx);
                phipart=phipart*180.0/3.14159265;
                if (phipart<0)
                  phipart= 360.0 - fabs(phipart);
              }
              if((myplaneside==0)&&(myhalftele==1)){
                phipart=atan2(lmyy,-lmyx);
                phipart=phipart*180.0/3.14159265;
                if (phipart<0)
                  phipart= 360.0 - fabs(phipart);
                ////std::cout<<"PS0 HT1  Reco: "<<padClv[k].GetClusterPhi()<<"; Geant: "<<phipart<<";    X-Y: "<<lmyx<<"-"<<lmyy<<std::endl;
              }
              diffphiCluHit->Fill(padClv[k].GetClusterPhi()-phipart);
            }
          }
        }                 // RPhiHitGeant

        for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
          phipart=(*p)->momentum().phi()*180.0/3.14159265;

          if (phipart<0)
            phipart= 360.0 - fabs(phipart);  // From particleGun Phi to Clusterizator Phi
          diffphiCluGun->Fill(padClv[k].GetClusterPhi()-phipart);
        }
      }
    }
*/

  // **************************************************** HIT STUDIES *********************************************** //
  // **************************************************** HIT STUDIES ************************************************ //
  // **************************************************** HIT STUDIES ************************************************ //
  // **************************************************** HIT STUDIES ************************************************ //
  // **************************************************** HIT STUDIES ************************************************ //

  // COMPARISON RECOHIT LOCALHIT
  // now we run for each plane the digi-simulation
//std::cout<<"hit+geant"<<std::endl;
/*
 if((trackcountergood==1)&&(singleparticle))
    for(map<int, PSimHitContainer>::const_iterator hitMapItr = hitMap.begin();  hitMapItr != hitMap.end(); ++hitMapItr){

      const PSimHitContainer & planeSimHits = hitMapItr->second;
      T2DetId *theT2DetId =new T2DetId(hitMapItr->first);

      myhalftele= theT2DetId->halfTelescope();
      myplaneside= theT2DetId->planeSide();
      myplane= theT2DetId->plane();
      zdetshift=myplane*planedist+myplaneside*btbdist;
      if (myhalftele==0)
        zdetshift=zdetshift+ovdist;
      zglobhit=zdetshift+z1+zinsidedet;

      for(T2HitCollection::const_iterator ithit = t2hitcoll->begin(); ithit != t2hitcoll->end(); ithit++){
        if (fabs(zglobhit-(ithit->GetHitZ()))<0.5)
          for (unsigned int l=0; l<planeSimHits.size(); l++){
            PSimHit myplanehits= planeSimHits[l];
            lmyx= myplanehits.localPosition().x();
            lmyy= myplanehits.localPosition().y();  //prima entryPoint()
            lmyr= sqrt(lmyx*lmyx + lmyy*lmyy);

            RLocHRecoH->Fill((ithit->GetHitR())-lmyr);
          }
      }
    }
*/
  //Studies on Reco Hit


  double hitphi=0.0;




  double phipar;
//std::cout<<"hit onlyt"<<std::endl;
  if(trackcountergood==1)
    for(T2HitCollection::const_iterator ithit = t2hitcoll->begin(); ithit != t2hitcoll->end(); ithit++) {
      ithit->GetHitR();
      hitphi=ithit->GetHitPhi();
      ithit->GetHitDR();
     ithit->GetHitDPhi();
      ithit->GetHitNumStrip();
     ithit->GetHitNumPad();

      for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
        phipar=(*p)->momentum().phi()*180.0/3.14159265;

        if (phipar<0)
          phipar= 360.0 + phipar;  // From particleGun Phi to Clusterizator Phi

        if ((fabs((*p)->momentum().eta())<6.5)&&(fabs((*p)->momentum().eta())>5.3)){
          if (fabs(phipar-hitphi)<25.0){
            //cout<<"fabs(phipar-hitphi)  ="<<fabs(phipar-hitphi)<<"   particle gen phi: "<<phipar<<"   Hit phi"<< hitphi<<endl;
            diffphiGUNHIT->Fill(phipar-hitphi);
          }
        }
      }
    }



if(chp12mag0>0)
  {
 MultiParticleEfficiencyCutNorm->Fill(chp12mag0,((double)numrecotrackcutmag0)/((double)chp12mag0));
 //cout<<" eta>0 Ch -- Tr:   "<<chp12mag0 <<" --  "<<numrecotrackcutmag0<<endl;
 MultiParticleEfficiencyCut->Fill(((double)chp12mag0),((double)numrecotrackcutmag0));
}


        //DECOMMENT

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
  SingleParticleEfficiencyCut
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



// ------------ method called once each job just before starting event loop  ------------const edm::EventSetup&
void T2RecoAnalyzer::beginJob()
{
  TH1::AddDirectory(kFALSE);

       
  // Chi2PhiProbLogy= std::auto_ptr<TCanvas>(new TCanvas("Chi2PhiProblog","Azimuthal #chi^{2} probability",400,400));
  // Chi2PhiProbLogy->SetDirectory(0);
  // Chi2RProbLogy= std::auto_ptr<TCanvas>(new TCanvas("Chi2RProblog","Radial #chi^{2} probability",400,400));
  //Chi2RProbLogy->SetDirectory(0);

  clusterstripentries = std::auto_ptr<TH1F>(new TH1F("clusterstripentries","Number of strips in strip-clusters",10, -0.5, 9.5));
  clusterstripentries->SetDirectory(0); 

  clusterpadentries = std::auto_ptr<TH1F>(new TH1F("clusterpadentries","Number of pads in pad-clusters",10, -0.5, 9.5));
  clusterpadentries->SetDirectory(0);
 
  diffphiGUNHIT= std::auto_ptr<TH1F>(new TH1F("diffphiGUNHIT","#Phi Generated particle - #Phi Reconstructed Hit",120, -30, 30));
  diffphiGUNHIT->SetDirectory(0); 
  NonTriggeredEnergyPlus= std::auto_ptr<TH1F>(new  TH1F("NonTriggeredEnergyPlus","NonTriggeredEnergyPlus",1000,0, 5000));

  NonTriggeredEnergyInPlusOnly_When100GeVDiscr= std::auto_ptr<TH1F>(new  TH1F("NonTriggeredEnergyInPlusOnly_When100GeVDiscr"," NonTriggeredEnergyInPlusOnly_When100GeVDiscr",1000,0, 5000));
  NonTriggeredEnergyInPlusAndMinus_When100GeVDiscr= std::auto_ptr<TH1F>(new  TH1F("NonTriggeredEnergyInPlusAndMinus_When100GeVDiscr","NonTriggeredEnergyInPlusAndMinus_When100GeVDiscr ",1000,0, 5000));

  RLocHRecoH= std::auto_ptr<TH1F>(new  TH1F("RLocHRecoH"," Reconstructed Hit R - GEANT4 Hit R",100,-1.5, 1.5));
  RLocHRecoH->SetDirectory(0); 

  diffRCluHit= std::auto_ptr<TH1F>(new  TH1F("diffRCluHit"," Cluster R - GEANT4 Hit R",100,-1.5, 1.5));
  diffRCluHit->SetDirectory(0); 

  Trketa= std::auto_ptr<TH1F>(new TH1F("Trketa","Reconstructed Track #eta (no cut used)",600,-8.0,8.0));
  Trketa->SetDirectory(0); 

  Trkphi= std::auto_ptr<TH1F>(new TH1F("Trkphi","Reconstructed Track #phi (no cut used)",361,0,360));
  Trkphi->SetDirectory(0); 

  DPhiGoodTrk= std::auto_ptr<TH1F>(new TH1F("DPhiGoodTrk","#phi Track - #phi GEANT4 Track (Only primary tracks are considered)",51,-25.5,25.5));
  DPhiGoodTrk->SetDirectory(0); 

  DEtaGoodTrk= std::auto_ptr<TH1F>(new TH1F("DEtaGoodTrk","#eta Track - #phi GEANT4 Track (Only primary tracks are considered)",200,-1.0,1.0));
  DEtaGoodTrk->SetDirectory(0); 

  diffphiCluGun = std::auto_ptr<TH1F>(new TH1F("diffphiCluGun","Cluster #Phi - #Phi Generated particle",81, -40.5, 40.5));
  diffphiCluGun->SetDirectory(0);
 
  diffphiCluHit= std::auto_ptr<TH1F>(new TH1F("diffphiCluHit","Cluster #Phi - GEANT4 Hit #Phi",50, -5, 5));
  diffphiCluHit->SetDirectory(0); 

  Trketagood= std::auto_ptr<TH1F>(new TH1F("Trketagood","Reconstructed Track #eta (all cut used)",300,3.0,8.0));
  Trketagood->SetDirectory(0); 

  Trkphigood= std::auto_ptr<TH1F>(new TH1F("Trkphigood","Reconstructed Track #phi (all cut used)",361,0,360));
  Trkphigood->SetDirectory(0); 

 
  Chi2RProb = std::auto_ptr<TH1F>(new TH1F("Chi2RProb","Radial #chi^{2} probability",15000,0.0,1.05));
  Chi2RProb->SetDirectory(0); 

  Chi2PhiProb = std::auto_ptr<TH1F>(new TH1F("Chi2PhiProb","Azimuthal #chi^{2} probability",15000,0.0,1.05));
  Chi2PhiProb->SetDirectory(0);

  SingleParticleEfficiencyCut = std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiencyCut","# Tracks reconstructed vs generated particle #eta (all cuts)",34,4.35,7.75,""));
  SingleParticleEfficiencyCut->SetDirectory(0);
 
  EfficiencyTrigT2PlusVSEplus= std::auto_ptr<TProfile>(new TProfile("EfficiencyTrigT2PlusVSEplus","T2+ Trigger Efficiency vs E ",1000,-0.5,999.5,""));

  NonCentralMassFromPom= std::auto_ptr<TH1D>(new TH1D("NonCentralMassFromPom","NonCentralMassFromPom",1000,-0.5,999.5));



ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr100GeV= std::auto_ptr<TH1D>(new TH1D("ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr100GeV","ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr100GeV",4000,-0.5,3999.5));
ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr200GeV= std::auto_ptr<TH1D>(new TH1D("ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr200GeV","ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr200GeV",4000,-0.5,3999.5));
 ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr300GeV= std::auto_ptr<TH1D>(new TH1D("ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr300GeV","ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr300GeV",4000,-0.5,3999.5));
  ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr400GeV= std::auto_ptr<TH1D>(new TH1D("ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr400GeV","ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr400GeV",4000,-0.5,3999.5));




CentralEnergyWhenfscT1T2OFF_An005Cut= std::auto_ptr<TH1D>(new TH1D("CentralEnergyWhenfscT1T2OFF_An005Cut","CentralEnergyWhenfscT1T2OFF_An005Cut",4000,-0.5,3999.5));
ForwEnergyWhenForwOFF_005Cut_Discr100GeV= std::auto_ptr<TH1D>(new TH1D("ForwEnergyWhenForwOFF_005Cut_Discr100GeV","ForwEnergyWhenForwOFF_005Cut_Discr100GeV",4000,-0.5,3999.5));
ForwEnergyWhenForwOFF_005Cut_Discr200GeV= std::auto_ptr<TH1D>(new TH1D("ForwEnergyWhenForwOFF_005Cut_Discr200GeV","ForwEnergyWhenForwOFF_005Cut_Discr200GeV",4000,-0.5,3999.5));
ForwEnergyWhenForwOFF_005Cut_Discr300GeV= std::auto_ptr<TH1D>(new TH1D("ForwEnergyWhenForwOFF_005Cut_Discr300GeV","ForwEnergyWhenForwOFF_005Cut_Discr300GeV",4000,-0.5,3999.5)); 
ForwEnergyWhenForwOFF_005Cut_Discr400GeV= std::auto_ptr<TH1D>(new TH1D("ForwEnergyWhenForwOFF_005Cut_Discr400GeV","ForwEnergyWhenForwOFF_005Cut_Discr400GeV",4000,-0.5,3999.5));


ForwEnergyWhenForwOFF_003_008_Cut_Discr100GeV= std::auto_ptr<TH1D>(new TH1D("ForwEnergyWhenForwOFF_003_008_Cut_Discr100GeV","ForwEnergyWhenForwOFF_003_008_Cut_Discr100GeV",4000,-0.5,3999.5));
ForwEnergyWhenForwOFF_003_008_Cut_Discr200GeV= std::auto_ptr<TH1D>(new TH1D("ForwEnergyWhenForwOFF_003_008_Cut_Discr200GeV","ForwEnergyWhenForwOFF_003_008_Cut_Discr200GeV",4000,-0.5,3999.5));
ForwEnergyWhenForwOFF_003_008_Cut_Discr300GeV= std::auto_ptr<TH1D>(new TH1D("ForwEnergyWhenForwOFF_003_008_Cut_Discr300GeV","ForwEnergyWhenForwOFF_003_008_Cut_Discr300GeV",4000,-0.5,3999.5)); 
ForwEnergyWhenForwOFF_003_008_Cut_Discr400GeV= std::auto_ptr<TH1D>(new TH1D("ForwEnergyWhenForwOFF_003_008_Cut_Discr400GeV","ForwEnergyWhenForwOFF_003_008_Cut_Discr400GeV",4000,-0.5,3999.5));

ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC100GeV= std::auto_ptr<TH1D>(new TH1D("ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC100GeV","ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC100GeV",4000,-0.5,3999.5));
ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC200GeV= std::auto_ptr<TH1D>(new TH1D("ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC200GeV","ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC200GeV",4000,-0.5,3999.5));
ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC300GeV= std::auto_ptr<TH1D>(new TH1D("ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC300GeV","ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC300GeV",4000,-0.5,3999.5)); 
ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC400GeV= std::auto_ptr<TH1D>(new TH1D("ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC400GeV","ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC400GeV",4000,-0.5,3999.5));


 Pom_minus_Cent_003_008= std::auto_ptr<TH1D>(new TH1D("Pom_minus_Cent_003_008","Pom_minus_Cent_003_008",4000,-0.5,3999.5));



  ForwRegionMassWhenT2OFF_Discr300GeV= std::auto_ptr<TH1D>(new TH1D("ForwRegionMassWhenT2OFF_Discr300GeV","ForwRegionMassWhenT2OFF_Discr300GeV",4000,-0.5,3999.5));

  ForwRegionMassWhenT2OFF_Discr200GeV= std::auto_ptr<TH1D>(new TH1D("ForwRegionMassWhenT2OFF_Discr200GeV","ForwRegionMassWhenT2OFF_Discr200GeV",4000,-0.5,3999.5));

  ForwRegionMassWhenT2OFF_Discr100GeV= std::auto_ptr<TH1D>(new TH1D("ForwRegionMassWhenT2OFF_Discr100GeV","ForwRegionMassWhenT2OFF_Discr100GeV",4000,-0.5,3999.5));
  
  ForwRegionMassWhenT2OFF_Discr400GeV= std::auto_ptr<TH1D>(new TH1D("ForwRegionMassWhenT2OFF_Discr400GeV","ForwRegionMassWhenT2OFF_Discr400GeV",4000,-0.5,3999.5));



 ForwRegionMassWhenT2OFF_Discr300GeVNolargeEta= std::auto_ptr<TH1D>(new TH1D("ForwRegionMassWhenT2OFF_Discr300GeVNolargeEta","ForwRegionMassWhenT2OFF_Discr300GeVNolargeEta",4000,-0.5,3999.5));

  ForwRegionMassWhenT2OFF_Discr200GeVNolargeEta= std::auto_ptr<TH1D>(new TH1D("ForwRegionMassWhenT2OFF_Discr200GeVNolargeEta","ForwRegionMassWhenT2OFF_Discr200GeVNolargeEta",4000,-0.5,3999.5));

  ForwRegionMassWhenT2OFF_Discr100GeVNolargeEta= std::auto_ptr<TH1D>(new TH1D("ForwRegionMassWhenT2OFF_Discr100GeVNolargeEta","ForwRegionMassWhenT2OFF_Discr100GeVNolargeEta",4000,-0.5,3999.5));
  
  ForwRegionMassWhenT2OFF_Discr400GeVNolargeEta= std::auto_ptr<TH1D>(new TH1D("ForwRegionMassWhenT2OFF_Discr400GeVNolargeEta","ForwRegionMassWhenT2OFF_Discr400GeVNolargeEta",4000,-0.5,3999.5));


  CentralEnergyWhenHigh005Csi= std::auto_ptr<TH1D>(new TH1D("CentralEnergyWhenHigh005Csi","CentralEnergyWhenHigh005Csi",4000,-0.5,3999.5));
CentralEnergyWhenHigh005Csi_AndfscT1T2OFF= std::auto_ptr<TH1D>(new TH1D("CentralEnergyWhenHigh005Csi_AndfscT1T2OFF","CentralEnergyWhenHigh005Csi_AndfscT1T2OFF",4000,-0.5,3999.5));
  
  CentralEnergyWhenfscT1T2OFF= std::auto_ptr<TH1D>(new TH1D("CentralEnergyWhenfscT1T2OFF","CentralEnergyWhenfscT1T2OFF",4000,-0.5,3999.5));
  PomeronInvMassWhenfscT1T2OFFvsCentralE= std::auto_ptr<TH2D>(new TH2D("PomeronInvMassWhenfscT1T2OFFvsCentralE","PomeronInvMassWhenfscT1T2OFFvsCentralE",4000,-0.5,3999.5,4000,-0.5,3999.5));

  
  T1chargeWhen_CentralEnergyWhenHigh005Csi_100= std::auto_ptr<TH1D>(new TH1D("T1chargeWhen_CentralEnergyWhenHigh005Csi_100","T1chargeWhen_CentralEnergyWhenHigh005Csi_100",50,-0.5,49.5));
  FSCchargeWhen_CentralEnergyWhenHigh005Csi_100= std::auto_ptr<TH1D>(new TH1D("FSCchargeWhen_CentralEnergyWhenHigh005Csi_100","FSCchargeWhen_CentralEnergyWhenHigh005Csi_100",50,-0.5,49.5));  
  T2TrkWhen_CentralEnergyWhenHigh005Csi_100= std::auto_ptr<TH1D>(new TH1D("T2TrkWhen_CentralEnergyWhenHigh005Csi_100","T2TrkWhen_CentralEnergyWhenHigh005Csi_100",50,-0.5,49.5));

  T1chargeWhen_CentralEnergyWhenHigh005Csi_30= std::auto_ptr<TH1D>(new TH1D("T1chargeWhen_CentralEnergyWhenHigh005Csi_30","T1chargeWhen_CentralEnergyWhenHigh005Csi_30",50,-0.5,49.5));
  FSCchargeWhen_CentralEnergyWhenHigh005Csi_30= std::auto_ptr<TH1D>(new TH1D("FSCchargeWhen_CentralEnergyWhenHigh005Csi_30","FSCchargeWhen_CentralEnergyWhenHigh005Csi_30",50,-0.5,49.5));  
  T2TrkWhen_CentralEnergyWhenHigh005Csi_30= std::auto_ptr<TH1D>(new TH1D("T2TrkWhen_CentralEnergyWhenHigh005Csi_30","T2TrkWhen_CentralEnergyWhenHigh005Csi_30",50,-0.5,49.5));  

  ForwardWhen_CentralEnergyWhenHigh005Csi_30= std::auto_ptr<TH1D>(new TH1D("ForwardWhen_CentralEnergyWhenHigh005Csi_30","ForwardWhen_CentralEnergyWhenHigh005Csi_30",72,-0.5,71.5));  
  ForwardWhen_CentralEnergyWhenHigh005Csi_100= std::auto_ptr<TH1D>(new TH1D("ForwardWhen_CentralEnergyWhenHigh005Csi_100","ForwardWhen_CentralEnergyWhenHigh005Csi_100",72,-0.5,71.5));  

  EtaMaxAroundLeftProton= std::auto_ptr<TH1D>(new TH1D("EtaMaxAroundLeftProton","EtaMaxAroundLeftProton",72,-18,18));
  NumEvtSel2Proton= std::auto_ptr<TH1D>(new TH1D("NumEvtSel2Proton","NumEvtSel2Proton",2,-0.5,1.5));
  AllRecoMassHist= std::auto_ptr<TH1D>(new TH1D("AllRecoMassHist","AllRecoMassHist",9000,-0.5,8999.5));
  forwInvMassHist= std::auto_ptr<TH1D>(new TH1D("forwInvMassHist","forwInvMassHist",1000,-0.5,999.5));
  CentInvMassHist= std::auto_ptr<TH1D>(new TH1D("CentInvMassHist","CentInvMassHist",6000,-0.5,5999.5));
  PomeronInvMassHist= std::auto_ptr<TH1D>(new TH1D("PomeronInvMassHist","PomeronInvMassHist",6000,-0.5,5999.5));
  ProtonRapidity= std::auto_ptr<TH1D>(new  TH1D("ProtonRapidity"," ProtonRapidity",100,-20.5, 20.5));
  csi1Histo= std::auto_ptr<TH1D>(new TH1D("csi1Histo","csi1Histo",500,0.,1.0));
  
  csi1_vs_csi2= std::auto_ptr<TH2D>(new TH2D("csi1_vs_csi2","csi1_vs_csi2",99,0.,1.,99,0.,1.));


  MultiplicityWhenT2ShouldHAVE= std::auto_ptr<TH1D>(new TH1D("MultiplicityWhenT2ShouldHAVE","MultiplicityWhenT2ShouldHAVE",200,-0.5,199.5));

CumulativeEtaGeneratorALLParticles= std::auto_ptr<TH1D>(new TH1D("CumulativeEtaGeneratorALLParticles","CumulativeEtaGeneratorALLParticles",1000,-9.,9.));

  CumulativeEtaGeneratorStableCHParticles= std::auto_ptr<TH1D>(new TH1D("CumulativeEtaGeneratorStableCHParticles","CumulativeEtaGeneratorStableCHParticles",1000,-9.,9.));
  CumulativeEtaGeneratorStableParticles= std::auto_ptr<TH1D>(new TH1D("CumulativeEtaGeneratorStableParticles","CumulativeEtaGeneratorStableParticles",1000,-9.,9.));
  CumulativeEtaGeneratorStableParticles_WhenACHARGEInT2= std::auto_ptr<TH1D>(new TH1D("CumulativeEtaGeneratorStableParticles_WhenACHARGEInT2","CumulativeEtaGeneratorStableParticles_WhenACHARGEInT2",1000,-9.,9.));


  CumulativeEtaGeneratorStableParticlesLeftOnly= std::auto_ptr<TH1D>(new TH1D("CumulativeEtaGeneratorStableParticlesLeftOnly","CumulativeEtaGeneratorStableParticlesLeftOnly",1000,-9.,9.)); 
  CumulativeEtaGeneratorStableParticlesRightOnly= std::auto_ptr<TH1D>(new TH1D("CumulativeEtaGeneratorStableParticlesRightOnly","CumulativeEtaGeneratorStableParticlesRightOnly",1000,-9.,9.));

  SingleParticleEfficiency = std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiency","# Tracks reconstructed vs generated particle #eta (all track)",34,4.35,7.75,""));
  SingleParticleEfficiency->SetDirectory(0); 

  MultiParticleEfficiencyCut= std::auto_ptr<TProfile>(new TProfile("MultiParticleEfficiencyCut","<# Tracks> vs # Charged Particles (all track cuts)",12,-0.5,11.5,""));
  MultiParticleEfficiencyCut->SetDirectory(0); 

  MultiParticleEfficiencyCutNorm= std::auto_ptr<TProfile>(new TProfile("MultiParticleEfficiencyCutNorm","<# Tracks / # Charged Particles > vs # Charged Particles (all track cuts)",12,-0.5,11.5,"")); 
  MultiParticleEfficiencyCutNorm->SetDirectory(0); 
   
  EventEffiXChpartInT2YNumTrk= std::auto_ptr<TH2D>(new TH2D("EventEffiXChpartInT2YNumTrk","EventEffiXChpartInT2YNumTrk",40,-0.5,39.5,40,-0.5,39.5));  

  TrackHemisphereNorequirement= std::auto_ptr<TH1F>(new TH1F("TrackHemisphereNorequirement","TrackHemisphereNorequirement",4,-1.5,2.5));
  TrackHemisphereNorequirement->SetDirectory(0); 
  
  RoadHemisphereNorequirement= std::auto_ptr<TH1F>(new TH1F("RoadHemisphereNorequirement","RoadHemisphereNorequirement",4,-1.5,2.5));
  RoadHemisphereNorequirement->SetDirectory(0); 
  
  TrackHemisphereWhenChargedParticle= std::auto_ptr<TH1F>(new TH1F("TrackHemisphereWhenChargedParticle","TrackHemisphereWhenChargedParticle",4,-1.5,2.5));
  TrackHemisphereWhenChargedParticle->SetDirectory(0); 

  PadCluMultiplVsPlane= std::auto_ptr<TProfile>(new TProfile("PadCluMultiplVsPlane","Pad Cluster multiplicity vs Plane",40,-0.5,39.5,""));
  PadCluMultiplVsPlane->SetDirectory(0);

  StripCluMultiplVsPlane= std::auto_ptr<TProfile>(new TProfile("StripCluMultiplVsPlane","Strip Cluster multiplicity vs Plane",40,-0.5,39.5,""));
  StripCluMultiplVsPlane->SetDirectory(0);

  NumRecoTrkWhenAtLeastAstableInT2= std::auto_ptr<TH1D>(new TH1D("NumRecoTrkWhenAtLeastAstableInT2","NumRecoTrkWhenAtLeastAstableInT2",51,-0.5,50.5));
  NumRecoTrkPerEvent= std::auto_ptr<TH1D>(new TH1D("NumRecoTrkPerEvent","NumRecoTrkPerEvent",51,-0.5,50.5));
numevent=0;

  NumTotRoadH0= std::auto_ptr<TH1D>(new TH1D("NumTotRoadH0","NumTotRoadH0",351,-0.5,350.5));
  NumTotRoadH1= std::auto_ptr<TH1D>(new TH1D("NumTotRoadH1","NumTotRoadH1",351,-0.5,350.5));
  NumTotRoadH2= std::auto_ptr<TH1D>(new TH1D("NumTotRoadH2","NumTotRoadH2",351,-0.5,350.5));
  NumTotRoadH3= std::auto_ptr<TH1D>(new TH1D("NumTotRoadH3","NumTotRoadH3",351,-0.5,350.5));
  TrackPlusWhenGenRapGapPlus= std::auto_ptr<TH1D>(new TH1D("TrackPlusWhenGenRapGapPlus","0: rap gap Preserved. 1: rap Gap destroyed",2,-0.5,1.5));
  GeneratorEvtHemisp= std::auto_ptr<TH1D>(new TH1D("GeneratorEvtHemisp","GeneratorEvtHemisp",4,-1.5,2.5));
  ProbGetLowMass= std::auto_ptr<TH1D>(new TH1D("ProbGetLowMass","ProbGetLowMass",2,-0.5,1.5));

  DifferenceMassFromProtons_andAllOtherParticles= std::auto_ptr<TH1D>(new TH1D("DifferenceMassFromProtons_andAllOtherParticles","DifferenceMassFromProtons_andAllOtherParticles",2000,-100.,100.));

  GeneratorEvtHemispCH= std::auto_ptr<TH1D>(new TH1D("GeneratorEvtHemispCH","GeneratorEvtHemispCH",4,-1.5,2.5));
  NonTriggeredEnergy= std::auto_ptr<TH1D>(new  TH1D("NonTriggeredEnergy"," NonTriggeredEnergy",1000,-0.5, 999.5));
  LeakT1T2= std::auto_ptr<TH1D>(new TH1D("LeakT1T2","LeakT1T2",2,-0.5,1.5));
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
  Trkphi->SetXTitle("Track #phi (deg)");
  DEtaGoodTrk->SetXTitle("#Delta #eta");
  diffphiCluGun->SetXTitle("#Delta #Phi (deg)");
  diffphiCluHit->SetXTitle("#Delta #Phi (deg)");
  DPhiGoodTrk->SetXTitle("#Delta #Phi (deg)");
  Trkphigood->SetXTitle("Track #phi (deg)");
  Trketagood->SetXTitle("Track #eta");
  SingleParticleEfficiency->SetXTitle("Particle #eta");
  SingleParticleEfficiencyCut->SetXTitle("Particle #eta");
  SingleParticleEfficiency->SetYTitle("Pseudo-Efficiency");
  SingleParticleEfficiencyCut->SetYTitle("Efficiency");
  MultiParticleEfficiencyCutNorm->SetYTitle("Efficiency");
  MultiParticleEfficiencyCutNorm->SetXTitle("ch. particle in T2");  
  MultiParticleEfficiencyCut->SetYTitle("<Track>");
  MultiParticleEfficiencyCut->SetYTitle("ch. particle in T2");    
  Chi2RProb->SetXTitle("Probability");
  Chi2PhiProb->SetXTitle("Probability");


  //Chi2RProbLogy->cd();
  //Chi2RProb->Draw("");
  //Chi2RProbLogy->SetLogy();
  //Chi2PhiProbLogy->cd();
  //Chi2PhiProb->Draw("");
  //Chi2PhiProbLogy->SetLogy();
  
}

// ------------ method called once each job just after ending the event loop  ------------

void T2RecoAnalyzer::endJob()
{
  
TFile *f = TFile::Open(outputFileName.c_str(), "recreate");
 if( !f || !f->IsWritable() ){
   std::cout << "Output file not opened correctly !!" << std::endl;
 }
 
   CumulativeEtaGeneratorALLParticles->Write();
CumulativeEtaGeneratorStableParticles_WhenACHARGEInT2->Write("");
CumulativeEtaGeneratorStableParticles->Write("");
CumulativeEtaGeneratorStableCHParticles->Write("");
 ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr100GeV->Write("");
ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr200GeV->Write("");
 ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr300GeV->Write("");
  ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr400GeV->Write("");

  CumulativeEtaGeneratorStableParticlesLeftOnly->Write(); 
  CumulativeEtaGeneratorStableParticlesRightOnly->Write();
  NumEvtSel2Proton->Write("");
  AllRecoMassHist->Write("");
  forwInvMassHist->Write("");
  CentInvMassHist->Write("");
  PomeronInvMassHist->Write("");
  csi1Histo->Write("");EtaMaxAroundLeftProton->Write("");csi1_vs_csi2->Write("");
  MultiplicityWhenT2ShouldHAVE->Write("");
  DifferenceMassFromProtons_andAllOtherParticles->Write("");
  T1chargeWhen_CentralEnergyWhenHigh005Csi_100->Write();
  FSCchargeWhen_CentralEnergyWhenHigh005Csi_100->Write();  
  T2TrkWhen_CentralEnergyWhenHigh005Csi_100->Write();

  T1chargeWhen_CentralEnergyWhenHigh005Csi_30->Write();
  FSCchargeWhen_CentralEnergyWhenHigh005Csi_30->Write();  
  T2TrkWhen_CentralEnergyWhenHigh005Csi_30->Write();  
  ForwardWhen_CentralEnergyWhenHigh005Csi_30->Write(); 
  ForwardWhen_CentralEnergyWhenHigh005Csi_100->Write(); 
    
  ForwRegionMassWhenT2OFF_Discr100GeVNolargeEta->Write("");
  ForwRegionMassWhenT2OFF_Discr200GeVNolargeEta->Write("");
  ForwRegionMassWhenT2OFF_Discr300GeVNolargeEta->Write("");
  ForwRegionMassWhenT2OFF_Discr400GeVNolargeEta->Write("");
  
  ForwRegionMassWhenT2OFF_Discr100GeV->Write("");
  ForwRegionMassWhenT2OFF_Discr200GeV->Write("");
  ForwRegionMassWhenT2OFF_Discr300GeV->Write("");
  ForwRegionMassWhenT2OFF_Discr400GeV->Write("");
  ProtonRapidity->Write("");

  CentralEnergyWhenfscT1T2OFF->Write("");CentralEnergyWhenHigh005Csi_AndfscT1T2OFF->Write("");
  PomeronInvMassWhenfscT1T2OFFvsCentralE->Write("");
NonCentralMassFromPom->Write("");
  LeakT1T2->Write("");
  clusterstripentries->Write("");
  clusterpadentries->Write("");
  diffphiGUNHIT->Write("");
  DPhiGoodTrk->Write("");
  DEtaGoodTrk->Write("");
  RLocHRecoH->Write("");
  CentralEnergyWhenHigh005Csi->Write("");
  diffRCluHit->Write("");
  diffphiCluGun->Write("");
  diffphiCluHit->Write("");
  SingleParticleEfficiencyCut->Write("");
  SingleParticleEfficiency->Write("");GeneratorEvtHemisp->Write("");
  Trketagood->Write("");
  Trkphigood->Write("");
  Trketa->Write("");
  Trkphi->Write("");
  TrackPlusWhenGenRapGapPlus->Write("");
  NumRecoTrkWhenAtLeastAstableInT2->Write("");
  NumRecoTrkPerEvent->Write("");
  EventEffiXChpartInT2YNumTrk->Write("");
  PadCluMultiplVsPlane->Write("");StripCluMultiplVsPlane->Write("");
  TrackHemisphereWhenChargedParticle->Write("");
  RoadHemisphereNorequirement->Write("");
  TrackHemisphereNorequirement->Write("");
  MultiParticleEfficiencyCutNorm->Write("");
  MultiParticleEfficiencyCut->Write("");
  Chi2PhiProb->Write(""); 
  CentralEnergyWhenHigh005Csi->Write(""); 
  Chi2RProb->Write(""); 
  //Chi2PhiProbLogy->Write(""); 
  //  Chi2RProbLogy->Write("");

  NonTriggeredEnergyInPlusOnly_When100GeVDiscr->Write(""); 
NonTriggeredEnergyInPlusAndMinus_When100GeVDiscr->Write(""); 
  GeneratorEvtHemisp->Write();GeneratorEvtHemispCH->Write();
 ProbGetLowMass->Write();
  NonTriggeredEnergy->Write("");  EfficiencyTrigT2PlusVSEplus->Write("");
  NumTotRoadH0->Write();
  NumTotRoadH1->Write();
  NumTotRoadH2->Write();
  NumTotRoadH3->Write();
  


  Pom_minus_Cent_003_008->Write();
  CentralEnergyWhenfscT1T2OFF_An005Cut->Write();
  ForwEnergyWhenForwOFF_005Cut_Discr100GeV->Write();
  ForwEnergyWhenForwOFF_005Cut_Discr200GeV->Write();
  ForwEnergyWhenForwOFF_005Cut_Discr300GeV->Write();
  ForwEnergyWhenForwOFF_005Cut_Discr400GeV->Write();
  ForwEnergyWhenForwOFF_003_008_Cut_Discr100GeV->Write();
  ForwEnergyWhenForwOFF_003_008_Cut_Discr200GeV->Write();
  ForwEnergyWhenForwOFF_003_008_Cut_Discr300GeV->Write();
  ForwEnergyWhenForwOFF_003_008_Cut_Discr400GeV->Write();
  ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC100GeV->Write();
  ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC200GeV->Write();
  ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC300GeV->Write();
  ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC400GeV->Write();
















  f->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(T2RecoAnalyzer);

 

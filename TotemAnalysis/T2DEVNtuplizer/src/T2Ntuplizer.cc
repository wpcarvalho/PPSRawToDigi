/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Mirko Berretti (mirko.berretti@gmail.com)
*
****************************************************************************/

#include "TotemAnalysis/T2DEVNtuplizer/interface/T2Ntuplizer.h"

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
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
//#include "FastSimulation/Particle/interface/ParticleTable.h"

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
  simTrackContainerLabel = iConfig.getParameter<edm::InputTag>("SimTrackContainerLabel");
  simVertexContainerLabel = iConfig.getParameter<edm::InputTag>("SimVertexContainerLabel");
  t2PadDigiCollectionLabel = iConfig.getParameter<edm::InputTag>("T2PadDigiCollectionLabel");
  t2StripDigiCollectionLabel = iConfig.getParameter<edm::InputTag>("T2StripDigiCollectionLabel");
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
  CluLabel = iConfig.getParameter<std::string>("CluLabel");
  HitLabel = iConfig.getParameter<std::string>("HitLabel");
  RoadLabel = iConfig.getParameter<std::string>("RoadLabel");
  TrackLabel= iConfig.getParameter<std::string>("TrackLabel");
  //RawDataName= iConfig.getParameter<std::string>("RawDataName");

  // T2CutsUtil.SetCuts(T2_TrkEtamin,T2_TrkEtaMAX,T2_trkMultimin,T2_trkMultiMAX,T2_DZMultiplicator,T2_PhiChiProbCut,T2_RChiProbCut,T2_QuarterUsed,IgnoredSmallAnglePrimarySlopeCut,_T2_usesXYtracking);
  //std::vector<int> qused;qused.push_back(0);qused.push_back(1);qused.push_back(2);qused.push_back(3);
  //
    //T2CutsUtil.SetCuts(4.5,7.5,4,11,2.,0.01,0.01,qused,0.001,true);
  /*
  Chicut = iConfig.getParameter<double>("Chicut");
  PhiChiProbCut = iConfig.getParameter<double>("PhiChiProbCut");
  RChiProbCut = iConfig.getParameter<double>("RChiProbCut");  
  energy=iConfig.getParameter<double>("FitdzEnergy");
  DZScale=iConfig.getParameter<double>("DZScale");
  tracketamin=iConfig.getParameter<double>("TrkEtamin");
  tracketamax=iConfig.getParameter<double>("TrkEtaMAX");
  singleparticle=iConfig.getParameter<bool>("singleparticle");
  IncludeT1=iConfig.getParameter<bool>("IncludeT1");
  */
  
  std::vector<int> qused;qused.push_back(0);qused.push_back(1);qused.push_back(2);qused.push_back(3);
  
  T2CutsUtil.SetCuts(4.5,7.5,4,11,2.,0.01,0.01,qused,0.001,true);
  LoadPrimaryGeantInformation=iConfig.getParameter<bool>("LoadPrimaryGeantInformation");
  EasyGeneratorPrimaryDef=iConfig.getParameter<bool>("EasyGeneratorPrimaryDef");
  HepMCProductLabel = iConfig.getParameter<std::string>("HepMCProductLabel");
  



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
   carica=1; //picharg
   break;
 
 case(-211):
   carica=-1; //picharg
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
   
   //Sigmas  and antisigmas
   
 case(3222):
   carica=1; 
   break;
   
 case(-3222):
   carica=-1;
   break;
   
 case(3112):
   carica=1; 
   break;
   
 case(-3112):
   carica=-1; 
   break;
  
   //CSI:
   
 case(3312):
   carica=-1; 
   break;
   
   
 case(-3312):
   carica=-1; 
   break;
  
   
   //Omegas:
 case(3334):
   carica=-1; 
   break;


 case(-3334):
   carica=-1; 
   break;
   
   
 }
 
 return carica;
}




bool IsinT2(double eta, double tracketamin, double tracketamax){
  bool flag;
  if ((fabs(eta)>tracketamin)&&(fabs(eta)<tracketamax))
    flag=true;
  else
    flag=false;
  return flag;
}




double T2Ntuplizer::GetEta2FromGeantHits(std::vector<std::pair<double,double> > SimHitRZ_ofGeantPrimaryTrk, double vtxZ)
{
  
  double thetaHitAverage=0.;     
  double hemis=0.; 
  for(unsigned int u=0;u<SimHitRZ_ofGeantPrimaryTrk.size();u++)
    {
      std::pair<double, double> SimHitRZ_PrimTrk=SimHitRZ_ofGeantPrimaryTrk.at(u);
      thetaHitAverage += fabs(SimHitRZ_PrimTrk.first/(SimHitRZ_PrimTrk.second-vtxZ));
      hemis=fabs(SimHitRZ_PrimTrk.second)/SimHitRZ_PrimTrk.second;
    }

  thetaHitAverage=thetaHitAverage/SimHitRZ_ofGeantPrimaryTrk.size();  
  double TrkEtaFromHitAverages_4Hits = (hemis*(-1.0)*log(thetaHitAverage/2.0)); //low angle approximation
  
  return TrkEtaFromHitAverages_4Hits;
    
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











T1T2Track T2Ntuplizer::RPhiFit(std::vector<T2Hit> &hitvec2)
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
    // if((hitvecR[jj].GetHitPhi()>0.)&&(hitvecR[jj].GetHitPhi()<10.))
    if((hitvecR[jj].GetHitPhi()>0.)&&(hitvecR[jj].GetHitPhi()<80.))
      num0++;

    //    if((hitvecR[jj].GetHitPhi()>350.)&&(hitvecR[jj].GetHitPhi()<360.))
    if((hitvecR[jj].GetHitPhi()>280.)&&(hitvecR[jj].GetHitPhi()<360.))
      num360++;
  }


 if((num0>0) && (num360>0)) 
   {
    
     for(unsigned int jj =0; jj<sizeHitvR; jj++)
       {	
	 if((hitvecR[jj].GetHitPhi()>0.)&&(hitvecR[jj].GetHitPhi()<80.))
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
  

  //std::cout<<"calc phim "<<phim*180.0/3.14159265<<" "<<hitvec2[0].GetHitPhi()<<" "<<hitvec[0].GetHitPhi()<<std::endl;
  
  T1T2Track fittedtrack(vect,mat,a_rz, b_rz, phim, e_a_rz, e_b_rz, e_phim, chi2,chi2r,chi2p, normchi2red, hemisphere,2);   
  return fittedtrack; 
}


double T2Ntuplizer::CalculateGeneratorEtaFromTrkId(std::vector<int> &ThePrimaryTrkIdBarcode,const HepMC::GenEvent* evt)
{

  double geneta=0.; 

  for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
    if((*p)->barcode()==ThePrimaryTrkIdBarcode.at(0))
      {
	geneta=((*p)->momentum().eta());
      }
  }

  return geneta;
  
}

int T2Ntuplizer::CaluclateCommonPDGSECMoth(std::vector<int> &TheSecondaryTrkIdFromRecHits,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, const HepMC::GenEvent* evt,bool printing, int& motherbarcode){
  int motherpdg=0;
  SimTrack atrk;
  int motherbarcodeRet=0;
  int thePidGeneratorOfThisTrkRet=0;
  
  int firstone=TheSecondaryTrkIdFromRecHits.at(0);
  int lastone=TheSecondaryTrkIdFromRecHits.at(TheSecondaryTrkIdFromRecHits.size()-1);
  int matchfirst=0;int matchlast=0;
  int trkid=0;
  for(unsigned int u=0;u<TheSecondaryTrkIdFromRecHits.size();u++){
    if(TheSecondaryTrkIdFromRecHits.at(u)==firstone)
      matchfirst++;
    if(TheSecondaryTrkIdFromRecHits.at(u)==lastone)
      matchlast++;
  }
  if(matchfirst>matchlast)
    trkid=firstone;
  else
   trkid=lastone;

  if((matchlast==0)&&(matchfirst==0)){
    motherpdg=0;
    return motherpdg;
  }

  //getTrack: Returns track corresponding to DirectMotherTrk_Id. Returns true if found, false otherwise.
  if(getTrack(trkid,theSimTracks,atrk) == true){
    motherpdg=GetTrkOlderMotherPid(atrk,theSimVertices,theSimTracks, evt,motherbarcodeRet,thePidGeneratorOfThisTrkRet,printing);
    motherbarcode=motherbarcodeRet;//atrk.genpartIndex();
  }
  else{
    motherpdg=0;
    motherbarcode=0;
    return motherpdg;
  }
  return motherpdg;
}
// int CaluclateCommonPDGSEC(std::vector<int> &TheSecondaryTrkIdPDGID);
//int GetTrkOlderMotherPid(SimTrack atrk,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, const HepMC::GenEvent* evt,int &barcodeMother,int &DirectTrkPID);


int T2Ntuplizer::CaluclateCommonPDGSEC(std::vector<int> &TheSecondaryTrkIdPDGID){
  int toret=-1;
  
  int firstone=TheSecondaryTrkIdPDGID.at(0);
  int lastone=TheSecondaryTrkIdPDGID.at(TheSecondaryTrkIdPDGID.size()-1);
  int matchfirst=0;int matchlast=0;
  for(unsigned int u=0;u<TheSecondaryTrkIdPDGID.size();u++){
    if(TheSecondaryTrkIdPDGID.at(u)==firstone)
      matchfirst++;
    if(TheSecondaryTrkIdPDGID.at(u)==lastone)
      matchlast++;
  }
  if(matchfirst>matchlast)
    toret=firstone;
  else
    toret=lastone;

  if((matchlast==0)&&(matchfirst==0))
    toret=-1;


  return toret;
}


T1T2Track T2Ntuplizer::MyLinearfitCorrDEV(std::vector<T2Hit> hitvec2,TMatrixD &par_covariance_matrix,double &chi2_,bool Usestrip, int RoadID){

  std::vector<T2Hit> hitvec; hitvec.clear();
   
 
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
 
  //if(verbosity){
  // std::cout<<" MyLinearfitCorrDEV start with "<<hitvec2.size()<<std::endl;
    
  //}

  /////////////////////////////////////////////////////////////////////////////////
  ///           SAVE IN HITVEC ONLY WANTED HITS
  ///////////////////////////////////////////////////////////////////////////////

  for(unsigned int m=0;m<hitvec2.size();m++)
    {
      inserted=false;
      if(hitvec2.at(m).GetHitNumPad()>0)
	{
	  if((hitvec2.at(m).GetHitPhi()<85)&&(hitvec2.at(m).GetHitPhi()>0))//was <90
	    numphimin20++;
	
	  if((hitvec2.at(m).GetHitPhi()>280))//was>340
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

	//if(verbosity)
	// std::cout<<"Z-Phi: "<<hitvec2.at(m).GetHitZ()<<" "<<hitvec2.at(m).GetHitPhi()<<" Num Pad-Strip:"<<hitvec2.at(m).GetHitNumPad()<<" "<<hitvec2.at(m).GetHitNumStrip()<<std::endl;


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

  //if(verbosity)
  // std::cout<<"MyLinearfitCorr Start computation  .. "<<std::endl;
  
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

  //if(verbosity)
  // std::cout<<"MyLinearfitCorr Parameter calculation  .. "<<std::endl;

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

 
 

 double chi2X=0.; double chi2Y=0.; 

   


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


  
  T1T2Track trk2rz=RPhiFit(hitvec);//RZFIT  
  
  double chi2Phi=trk2rz.ChiSquaredPhi();
  double chi2R=trk2rz.ChiSquaredR();
  
  double phiRZ=trk2rz._phiRZ;               //trk2rz.Phi(); 
  double e_phiRZ=trk2rz._e_phiRZ;
  
  double TanthetaRZ=trk2rz._TanthetaRZ;     //trk2rz.GetTy(); 
  double e_TanthetaRZ=trk2rz._e_TanthetaRZ; //trk2rz.GetTySigma();
  
  double bRZ=trk2rz._bRZ;          //trk2rz.GetTx(); 
  double e_bRZ=trk2rz._e_bRZ;      //trk2rz.GetTxSigma();
  
  if((fabs(vect[0])>1.8)||(fabs(vect[2])>1.8))
    {
      std::cout<<"WARNING in T2TrackProducer3 Ax:"<<(vect[0])<<"  Ay:"<<(vect[2])<<" too big!! Trk Hit (all) print :"<<std::endl;
      for(unsigned int jj =0; jj<hitvec.size(); jj++)
	{
	  std::cout<<"x-y-z-phi:  "<<hitvec[jj].GetHitX()<<"  "<<hitvec[jj].GetHitY()<<"  "<<hitvec[jj].GetHitZ()<<"  "<<hitvec[jj].GetHitPhi()<<"     DX:"<<hitvec[jj].GetHitDX()<<"     DY:"<<hitvec[jj].GetHitDY()<<" NumPad:"<<hitvec[jj].GetHitNumPad()<<std::endl;
	}
    }

  //if(verbosity){
  //  std::cout<<"Before save: XY:"<<std::endl;
    // std::cout<<"ax.bx.ay.by:"<<vect[0]<<" "<<vect[1]<<" "<<vect[2]<<" "<<vect[3]<<std::endl;
  
  //  std::cout<<"Before save: RZ:"<<std::endl;
  //  std::cout<<"Phi "<<phiRZ*180.0/3.14159<<" TanTheta"<<TanthetaRZ<<" Brz:"<<bRZ<<" Eta:"<<trk2rz.Eta()<<std::endl;

  //  std::cout<<" TotHit:"<<numHit_t2<<" cl1Hit:"<<numCl1HitHit_t2<<" StripOnly:"<<numStripOnlyHit_t2<<" PadOnly:"<<numPadOnlyHit_t2<<std::endl;
  // }


  T1T2Track fittedtrack(FittedParamConv,ParCov_matrConv,chi2,chi2X,chi2Y,hemisphere,2 , TanthetaRZ,  bRZ,   phiRZ,   e_TanthetaRZ,  e_bRZ,  e_phiRZ, chi2R, chi2Phi, RoadID, numHit_t2, numCl1HitHit_t2, numStripOnlyHit_t2, numPadOnlyHit_t2);
  

  //if(verbosity){      
  // std::cout<<"After save: :"<<std::endl;
  //std::cout<<"PhiRZ:"<<fittedtrack._phiRZ*180.0/3.14159<<"  EtaRZ:"<< fittedtrack._etaRZ<<"  EtaXY:"<<fittedtrack.Eta()<<" ChirR/N:"<<(chi2R/(numCl1HitHit_t2+numPadOnlyHit_t2))<<" Chi2/N; "<<chi2/(numCl1HitHit_t2+numPadOnlyHit_t2)<<" NumHit:"<<fittedtrack._numHit_t2<<std::endl;
    
  // }


  sizeHitv=hitvec.size();
  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      fittedtrack.AddHit(hitvec[jj]);
    }
  
 // std::cout<<"Trk with "<<sizeHitv<<" hit saved"<<std::endl;
  return fittedtrack;
  
}





T1T2Track T2Ntuplizer::CalculateRecoTrkFromGeant(std::vector<Local3DPoint> PrimaryGeantHitPos, std::vector<double> &theZ,  std::vector<uint32_t> &thecmsswId){
  
  // std::cout<<"Inside CalculateRecoTrkFromGeantA"<<PrimaryGeantHitPos.size()<<std::endl;
  
    TMatrixD covmat(4,4);
    covmat.Zero();
    double chi2corr=0.;		
    
    T2Hit worsthit;
    unsigned int roadId=0;//theroad.RoadID;
    std::vector<T2Hit> hitvector;
    for(unsigned int i=0;i<PrimaryGeantHitPos.size();i++)
      {
	double X=PrimaryGeantHitPos.at(i).x();
	double Y=PrimaryGeantHitPos.at(i).y();
	double Z=theZ.at(i);
	double EX=2.;
	double EY=2.;
	double EZ=1.;
	uint32_t mycmsswid=thecmsswId.at(i);
	T2Hit arechit(X,  Y,  Z,  EX,  EY, EZ, mycmsswid);
	arechit.SetHitClass(1);
	arechit.SetHitNumStrip(1);
	arechit.SetHitNumPad(1);
	//set other properties of pad-strip composition
	hitvector.push_back(arechit);
      }
    // std::cout<<"Inside CalculateRecoTrkFromGeantB"<<hitvector.size()<<std::endl;
     
    T1T2Track firstTrkRaw=MyLinearfitCorrDEV(hitvector,covmat,chi2corr,false,roadId);  
    return firstTrkRaw;
    //std::cout<<"CDEBUG"<<std::endl;
}






void T2Ntuplizer::GeantTrkAnalyzer(const std::auto_ptr<edm::PSimHitContainer>& theSimHits,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, const HepMC::GenEvent* evt, std::vector<double> &PrimGeantTrkEta2InT2, std::vector<int> &PrimGeantTrkHemiInT2,std::vector<double> &PrimGeantTrkPTInT2,std::vector<int> &PrimGeantTrkBarcodeInT2,double vtxZ){
  
  unsigned int Trkhitcounter=0;
  unsigned int TrkCounter=0;
  
  for(edm::SimTrackContainer::const_iterator ittrack = theSimTracks->begin();ittrack != theSimTracks->end(); ++ittrack){
    
    int trackId = ittrack->trackId();
    //std::map<unsigned int, std::list<Local3DPoint> >::iterator checkIsT2Trk;
    Trkhitcounter=0;
    TrkCounter++;
   
  // std::cout<<"New Gtrack ysis:"<<std::endl;
  // std::cout<<"a)Geant Trk PID Printing: "<<std::endl;

  
    //This Trk Id has some PSimHit Found and saved in trackHitList
    //If It is not true, like in the SimHit generated in Pi0 Gun, you cannot do more

    //double trkenergy=0.1; 
    double trkpt=0.01;
    trkpt=ittrack->momentum().pt();
    //trkenergy=
    ittrack->momentum().e();
    //  std::cout<<"Gtrk E:"<<trkenergy<<" PID:"<<ittrack->type()<<std::endl;
    std::vector<std::pair<double,double> > SimHitRZ_ofGeantTrk_H0;
    std::vector<std::pair<double,double> > SimHitRZ_ofGeantTrk_H1;
   std::vector<std::pair<double,double> > SimHitRZ_ofGeantTrk_H2;
    std::vector<std::pair<double,double> > SimHitRZ_ofGeantTrk_H3;
 
    std::vector<std::pair<double,double> > SimHitRZ_ofGeantTrk_ForEff;
    std::vector<unsigned int> PadsRowForSameTrk;
    std::vector<unsigned int> PadsColForSameTrk;
    std::vector<unsigned int> StripRowForSameTrk;
    std::vector<unsigned int> StripColForSameTrk;
    T2GeometryUtil::T2DetInfo detInfo2;
    unsigned int NumT2HitinGeantTrkPrimary=0;
    //  std::cout<<"B"<<std::endl;
    // std::cout<<" trackId: "<<trackId<<std::endl;

    
    SimTrack PythiaIPChParticleTrack=(*ittrack);
    int motherbarcode=99; 
   
  
    bool trkprimcond=PrimaryTrackCondition(PythiaIPChParticleTrack,theSimVertices,theSimTracks,evt,motherbarcode);
    // std::cout<<k:"<<trkprimcond<<std::endl;

    if(trkprimcond){
      // std::cout<<"New Prim Track ysis"<<std::endl;

      int GeneratorBarcode=PythiaIPChParticleTrack.genpartIndex();
     

    for(edm::PSimHitContainer::const_iterator ithit = theSimHits->begin();ithit != theSimHits->end(); ++ithit){
	        
      int ParenttrackId = ithit->trackId();  
      // if(TrkCounter==1)
      //std::cout<<"ALL HitTrkID vs first TrkID:"<<ParenttrackId<<"-"<<  trackId<<std::endl;
      //      std::cout<<"New Prim Track ysis-A"<<std::endl;
      if(ParenttrackId == trackId)
	{
	  // T2GeometryUtil::T2DetInfo detInfo2;
 
	  Trkhitcounter++;
	  NumT2HitinGeantTrkPrimary++;
 
	  T2GeometryUtil conv;  
	  Local3DPoint localPosition = ithit->localPosition();   
	  unsigned int detUnitId = ithit->detUnitId();
	  T2Geometry t2geom;t2geom.setPlane(detUnitId);
	  detInfo2=conv.GetT2Info(detUnitId);
	  SimVertex simuvertex; 
	   

	  //std::cout<<"New Prim Track analysis-B1"<<std::endl;
	  double theX=localPosition.x();
	  double theY=localPosition.y();
	  double theR=sqrt(theX*theX+theY*theY);
	  //std::cout<<"Hit xyr:"<<theX<<" "<<theY<<" "<<theR<<std::endl;
	  
	  double theZ=detInfo2.Zdet;
	  std::pair<double, double> SimHitRZ_PrimTrk_ForEff(theR,theZ);
	  SimHitRZ_ofGeantTrk_ForEff.push_back(SimHitRZ_PrimTrk_ForEff);
	  //if (getVertex(ittrack->vertIndex(), theSimVertices,simuvertex) == true) //This vertex exists
	  //if (simuvertex.parentIndex() < 0){
	  // Local3DPoint localPositione = ithit->localPosition();		  		  	       
	  std::pair<double, double> SimHitRZ_PrimTrk(theR,theZ);
	  if((theZ>0)&&((detInfo2.symb/10)==0))
	    SimHitRZ_ofGeantTrk_H0.push_back(SimHitRZ_PrimTrk);

	  if((theZ>0)&&((detInfo2.symb/10)==1))
	    SimHitRZ_ofGeantTrk_H1.push_back(SimHitRZ_PrimTrk);
	  
	  if((theZ<0)&&((detInfo2.symb/10)==2))
	    SimHitRZ_ofGeantTrk_H2.push_back(SimHitRZ_PrimTrk);

	  if((theZ<0)&&((detInfo2.symb/10)==3))
	    SimHitRZ_ofGeantTrk_H3.push_back(SimHitRZ_PrimTrk);
	  	
	  // std::cout<<"New Prim Track analysis-B2"<<std::endl;
	  
	     	      	      	      
	}//if(ParenttrackId == trackId)

      //std::cout<<"New Prim Track analysis-C"<<std::endl;  
   	  
    }//end simHit loop Now you are again in the simTrk loop


    
    
    bool assoctogen=false;
    for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
      if(GeneratorBarcode==(*p)->barcode()){
	assoctogen=true;
	//	GeneratorEtaAssocToGeantPrim_DEBUG->Fill((*p)->momentum().eta());
      }
    }

    if(assoctogen){

      if((SimHitRZ_ofGeantTrk_H0.size()>=4)&&(NumT2HitinGeantTrkPrimary>=4)){
	//include this trk in the dn/deta of the prim geant trk.
	double eta2geant=GetEta2FromGeantHits(SimHitRZ_ofGeantTrk_H0,vtxZ);
	PrimGeantTrkEta2InT2.push_back(eta2geant); PrimGeantTrkPTInT2.push_back(trkpt); 
	PrimGeantTrkHemiInT2.push_back(0);  PrimGeantTrkBarcodeInT2.push_back(GeneratorBarcode);
      }
      
      if((SimHitRZ_ofGeantTrk_H1.size()>=4)&&(NumT2HitinGeantTrkPrimary>=4)){
	//include this trk in the dn/deta of the prim geant trk.
	double eta2geant=GetEta2FromGeantHits(SimHitRZ_ofGeantTrk_H1,vtxZ);
	PrimGeantTrkEta2InT2.push_back(eta2geant); PrimGeantTrkPTInT2.push_back(trkpt); 
	PrimGeantTrkHemiInT2.push_back(1);  PrimGeantTrkBarcodeInT2.push_back(GeneratorBarcode);
      }
    
    
      if((SimHitRZ_ofGeantTrk_H2.size()>=4)&&(NumT2HitinGeantTrkPrimary>=4)){
	//include this trk in the dn/deta of the prim geant trk.
	double eta2geant=GetEta2FromGeantHits(SimHitRZ_ofGeantTrk_H2,vtxZ);
	PrimGeantTrkEta2InT2.push_back(eta2geant); PrimGeantTrkPTInT2.push_back(trkpt); 
	PrimGeantTrkHemiInT2.push_back(2);  PrimGeantTrkBarcodeInT2.push_back(GeneratorBarcode);
      }

      if((SimHitRZ_ofGeantTrk_H3.size()>=4)&&(NumT2HitinGeantTrkPrimary>=4)){
	//include this trk in the dn/deta of the prim geant trk.
	double eta2geant=GetEta2FromGeantHits(SimHitRZ_ofGeantTrk_H3,vtxZ);
	PrimGeantTrkEta2InT2.push_back(eta2geant); PrimGeantTrkPTInT2.push_back(trkpt); 
	PrimGeantTrkHemiInT2.push_back(3);  PrimGeantTrkBarcodeInT2.push_back(GeneratorBarcode);
      }
    }
    
    
    }
  }
  
}



bool T2Ntuplizer::CalculateGeantEtaFromTrkId(std::vector<int> &ThePrimaryTrkIdFromRecHits,const std::auto_ptr<edm::PSimHitContainer>& theSimHits,
						       const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,T1T2Track &TheGeantTrkReconstr)
{

  unsigned int geantHitCount=0;
  T2GeometryUtil conv;      T2Geometry t2geom;
  std::vector<unsigned int> GeantHiCmsswId;//Utilized for discrepancy studies.
  std::vector<double> GeantHitRefZ;//Utilized for discrepancy studies.
  std::vector<Local3DPoint> PrimaryGeantHitPos;

  for(edm::PSimHitContainer::const_iterator simHitIt = theSimHits->begin(); simHitIt != theSimHits->end(); ++simHitIt){
  
     int trackId = simHitIt->trackId();
   
    if(trackId==ThePrimaryTrkIdFromRecHits.at(0))   
      {
	geantHitCount++;
	unsigned int detUnitId = simHitIt->detUnitId();    
	T2GeometryUtil::T2DetInfo detInfo = conv.GetT2Info(detUnitId);
	t2geom.setPlane(detUnitId); //getNearestPadRow and col does not use the det info. Are refered to one particular plane.
    
	Local3DPoint hitPos = simHitIt->localPosition();    
	Local3DPoint corPos = CalculateCorrectXYZPos(Local3DPoint(hitPos.x(),hitPos.y(),hitPos.z()),detUnitId);
    
	//int padRow = t2geom.getNearestPadRow(&corPos);
	//int padCol = t2geom.getNearestPadCol(&corPos);
	//int padRowReal = t2geom.getNearestPadRow(&hitPos);
	//int padColReal = t2geom.getNearestPadCol(&hitPos);
	//int padSymbDET =detInfo.symb;
	//PadRowHitsGeantTrk.push_back(padRow); 
	//PadColHitsGeantTrk.push_back(padCol); 
	//PadSymbdet_vect.push_back(padSymbDET);
	


	  //It is assumed that the default track collection is the one you do the test
	  //The program only look at the first two and last two planes, this is assumed to be the reference.
	//int TherefSymbplane=RawtoSymb(detUnitId);
	//std::vector <int>::iterator iter;
	//iter = std::find( refSymbplane.begin(), refSymbplane.end(), TherefSymbplane);
	//double radiusGhit=sqrt(hitPos.x()*hitPos.x()+hitPos.y()*hitPos.y());
	//if(iter != refSymbplane.end())	    
	{
	  //Store the hitPosIn a vector.
	  //GeantHitRefX;
	  //GeantHitRefY;
	  //GeantHitRefZ;
	  double theZ=detInfo.Zdet;
	  //std::cout<<"Inside CalculateGeantEtaFromTrkId.  Geant-Hit eta2:"<<(-1.0)*log((radiusGhit/theZ)/2.0)<<std::endl;
	  GeantHitRefZ.push_back(theZ);
	  GeantHiCmsswId.push_back(detUnitId);
	  PrimaryGeantHitPos.push_back(corPos);
	  // std::cout<<TherefSymbplane<<" Included    |";
	}
	// else
	// std::cout<<TherefSymbplane<<"Not Inc    |";

	

      }

    
    
  }

  bool GeantTrkFound=false;
   if(PrimaryGeantHitPos.size()>=4){
     //   std::cout<<"Calling From GeantTrkReconstructed the recon:"<<PrimaryGeantHitPos.size()<<std::endl;
     TheGeantTrkReconstr =CalculateRecoTrkFromGeant(PrimaryGeantHitPos,GeantHitRefZ,GeantHiCmsswId);
    //std::cout<<"Geant Trk recon with "<<TheGeantTrkReconstr.GetHitEntries()<<std::endl;
     GeantTrkFound=true;
   }
   else
     GeantTrkFound=false;
  
   //std::cout<<"End CalculateGeantEtaFromTrkId"<<std::endl;
   return GeantTrkFound;
   
}



int T2Ntuplizer::GetTrkOlderMotherPid(SimTrack aTrack,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
					const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,const HepMC::GenEvent* evt, int &motherbarcode,int &thePidGeneratorOfThisTrk, bool printing)
 {
  
  /////////////////////////////////
  int OlderMotherParticle_Id=0; //This will be returned.
  /////////////////////////////////

  unsigned int count=0;
  int OlderMotherId=0; 
  SimTrack DirectMotherTrk=aTrack;
  int DirectMotherTrk_Id=aTrack.trackId();

  if(fabs(aTrack.type())!=11)
    printing=false;

  //genpartIndex : index of the corresponding Generator particle in the Event container (-1 if no Genpart) 
  if(printing){
    std::cout<<"HepMC event Print for GeantTrk Pdg:"<<aTrack.type()<<"(HepMC) Barcode:"<<aTrack.genpartIndex()<<std::endl;
    //    evt->print();
  }
  
  for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
    if(printing)
      std::cout<<"Barcode: "<<((*p)->barcode())<<" PdgID:"<<(*p)->pdg_id()<<" status_code:"<<(*p)->status()<<std::endl;
    if((*p)->barcode()==aTrack.genpartIndex())
      {
	// thePidGenerator=(*p)->pdg_id();
	OlderMotherId=((*p)->barcode());
	motherbarcode=((*p)->barcode());
	thePidGeneratorOfThisTrk=(*p)->pdg_id();	
	OlderMotherParticle_Id=(*p)->pdg_id();		
      }
  }
  if(printing)
    std::cout<<"GetTrkOlderMotherPid begin. Trk Id: "<<DirectMotherTrk_Id<<" Pythia pdg (0 if not yet found in HepMC): "<<OlderMotherParticle_Id<<" Trk_PID:"<<(DirectMotherTrk.type())<<std::endl;

  if(printing)
    std::cout<<"Try to find mother for Particle_ID: "<<(aTrack.type())<<" with Trk_ID:"<<DirectMotherTrk_Id<<", Barcode="<<DirectMotherTrk.genpartIndex()<<std::endl;

 while(DirectMotherTrk_Id!=(-1))//loop until you dont' find the stable mother or the source
   {
     count++;     

     SimVertex simuvertex;
 
     //getTrack: Returns track corresponding to DirectMotherTrk_Id. Returns true if found, false otherwise.
  
     if(getTrack(DirectMotherTrk_Id,theSimTracks,DirectMotherTrk) == true)
       {
	 if(printing)
	   std::cout<<"Retrieved geant Trk with Particle_PID: "<<(DirectMotherTrk.type())<<" and Trk_ID:"<<DirectMotherTrk_Id<<", Barcode="<<DirectMotherTrk.genpartIndex()<<std::endl;
	 motherbarcode=DirectMotherTrk.genpartIndex();
	 OlderMotherId=DirectMotherTrk_Id;
	 OlderMotherParticle_Id=DirectMotherTrk.type();
	 //std::cout<<"Mother Trk id= "<<DirectMotherTrk_Id<<std::endl;
	 
	 //This look only if parentVtx is primary (Index<0) 
	 if(isTrackPrimary(DirectMotherTrk_Id, theSimTracks, theSimVertices)==false)
	   {	     
	     if(printing)
	       std::cout<<"Not a primary track. Try to find originating vertex..."<<std::endl;
	     
	     //Take the vtx generating the trk.
	     //Note :vertIndex: index of the vertex in the Event container (-1 if no vertex)

	     if(getVertex(DirectMotherTrk.vertIndex(), theSimVertices,simuvertex) == true) 
	       {	    
		 if(printing)
		   std::cout<<"Originating Geant vertex found. Try to found the ancestors of this track:"<<std::endl;
		 // Note that DirectMotherTrk.genpartIndex() is the barcode, 
		 // not the index within GenParticleCollection, so I have to search the particle
		 

		 
		  //Find the mother track originating this vtx 
		  //This instruction allows to go back in the three to found the origin
		   //THIS IS THE ONLY IMPORTANT INSTRUCTION
		 //parentIndex:G4 TrackId of the parent in the Event  SimTrack  container (-1 if no parent) BE CAREFUL this is not a vector index 
		  DirectMotherTrk_Id=simuvertex.parentIndex(); 
		  if(printing)		  
		    std::cout<<"Mother Geant-Trk Id originating this vtx:"<<DirectMotherTrk_Id<<std::endl;
		
	       }
	     else
	        if(printing)
		  std::cout<<"Originating vertex NOT found:"<<std::endl;
	   }
	 else //Track reached backword (or prompt) is primary
	   { 
	    //--  std::cout<<"isTrackPrimary==true"<<std::endl;
	      //In this case DirectMotherTrk_Id should be <0 .
	     if(printing)
	       std::cout<<"isTrackPrimary==true so Mother Trk is < 0, Any other Geant link is present. I want to see if some other parent can be found From generator:"<<std::endl;
	     // Note that DirectMotherTrk.genpartIndex() is the barcode, not the index within GenParticleCollection, so I have to search the particle
	     for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
	       if((*p)->barcode()==DirectMotherTrk.genpartIndex())
		 {
		   (*p)->pdg_id();
		   OlderMotherParticle_Id=(*p)->pdg_id();
		   motherbarcode=DirectMotherTrk.genpartIndex();
		   
		     //   std::cout<<"Pythia Ancestors Particle-ID : ";
		   HepMC::GenVertex::particle_iterator ancestor;
		   ancestor = (*p)->production_vertex()->particles_begin(HepMC::parents);
		   for(;ancestor != (*p)->production_vertex()->particles_end(HepMC::parents); ++ancestor) {
		     
		     motherbarcode=(*ancestor)->barcode();
		 
		       
		       if(((std::find(knowPart.begin(),knowPart.end(),fabs((*ancestor)->pdg_id()))!=knowPart.end()))||((*ancestor)->pdg_id()==111))//Warning: Pi0 111 accepted from 25/10/2011
		       {
			 OlderMotherId=((*ancestor)->barcode());
			 OlderMotherParticle_Id=(*ancestor)->pdg_id();
			 motherbarcode=(*ancestor)->barcode();
			 if(printing)
			   std::cout<<"Ancestor is a Known particle: Barcode:"<<motherbarcode<<" PID:"<<OlderMotherParticle_Id<<std::endl;
		       }
		       else
			 if(printing)
			   std::cout<<"Ancestor Pdg-id "<<(*ancestor)->pdg_id()<<" rejected"<<std::endl;
		     //otherwise keep what you already have, since pythia is not giving any true particle
		   }
		   		   		   		   
		 }
	       //if(printing)
	       //std::cout<<"Begin-2"<<motherbarcode<<std::endl;
	     }
	     
	     //std::cout<<"Exit condition active since Mother Trk was primary"<<std::endl;
	     DirectMotherTrk_Id=-1;//Exit condition active
	   }
       }//If you are able to Get the track from the input or suggested as a parent-vtxID.
     else
      DirectMotherTrk_Id=-1;//Exit condition active  but you have not recovered the trk  
     
   }//while end

 if(printing)
   std::cout<<"GetTrkOlderMotherPid Return Trk_Id:"<<OlderMotherId<<std::endl;

 return OlderMotherParticle_Id;

}



double T2Ntuplizer::MostProbableTrkEnergy(std::vector<double> &PrimaryTrkenergy_ForTheRecoHits){

  double theE2=0.;
  for (unsigned int h=0;h<PrimaryTrkenergy_ForTheRecoHits.size();h++)
    {
      theE2=theE2+PrimaryTrkenergy_ForTheRecoHits.at(h);
    }
  if(theE2>0.)
    theE2=theE2/PrimaryTrkenergy_ForTheRecoHits.size();
  else
    theE2=-1.;


  std::vector<std::vector<double> > SameEnergyHits;

  for(unsigned int y=0; y<PrimaryTrkenergy_ForTheRecoHits.size();y++){

    bool asameenergyHitfound=false;

    if(y==0){
      std::vector<double> EnergyV;
      EnergyV.push_back(PrimaryTrkenergy_ForTheRecoHits.at(y));
      SameEnergyHits.push_back(EnergyV);
    }
    else{
      
      for(unsigned int i=0;i<SameEnergyHits.size();i++)
	for(unsigned int j=0;j<SameEnergyHits.at(i).size();j++)
	  if(fabs(PrimaryTrkenergy_ForTheRecoHits.at(y)-SameEnergyHits.at(i).at(j))<0.001){ //was 1 until 12/12/2011
	    
	    SameEnergyHits.at(i).push_back(PrimaryTrkenergy_ForTheRecoHits.at(y));
	    asameenergyHitfound=true;
	    goto nexthit;
	  }
    }
    
   
  nexthit: // goto label
   
    if(asameenergyHitfound==false)
      {
	std::vector<double> EnergyV;
	EnergyV.push_back(PrimaryTrkenergy_ForTheRecoHits.at(y)); 
	SameEnergyHits.push_back(EnergyV);
      }
    
    
  }

  
for(unsigned int i=0;i<SameEnergyHits.size();i++)
  if(SameEnergyHits.at(i).size()>(0.7*PrimaryTrkenergy_ForTheRecoHits.size())) { }
   // theE=SameEnergyHits.at(i).at(0);



  return theE2;
}


bool T2Ntuplizer::PrimaryTrackCondition(SimTrack atrk,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, const HepMC::GenEvent* evt,int &barcodeMother){


  bool primarycondition=false;
  
  int trackId=atrk.trackId();
  //SimTrack is stored and the generating vertex is at IP.
  
  bool isTrackPrimaryIP=false;
  if(fabs(atrk.charge())>0){
    isTrackPrimaryIP = isTrackPrimary(trackId, theSimTracks, theSimVertices);
  }
    
    //SimTrack PythiaIPChParticleTrack=(*ittrack);      
   int thePidGenerator=0; int thePidGeneratorOfThisTrk=0;
  
   thePidGenerator=GetTrkOlderMotherPid(atrk,theSimVertices,theSimTracks,evt,barcodeMother,thePidGeneratorOfThisTrk,false);

  if(thePidGeneratorOfThisTrk==0) //It means that Geant trk was not associated to a pythia particle
    thePidGeneratorOfThisTrk=atrk.type();


  bool alloweddecayparticle=false;
  if(std::find(AllowedDecayToCount.begin(),AllowedDecayToCount.end(),fabs(thePidGenerator))!=AllowedDecayToCount.end())
    alloweddecayparticle=true;

  //AllowedDecayToCount should only contains  particles produced by decays after ctau>1cm
  //I should remove decay products of neutral hadrons which decay in charged particles before 1 cm.
 alloweddecayparticle=false;
  //I'm Excluding 211 also because I assume that a nonprimary older mother 211 have done Nuclear interaction


  
 bool assoctogen=false;
 for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
   if(atrk.genpartIndex()==(*p)->barcode())
     if((*p)->status()==1)
       {
	 assoctogen=true;
	 //	GeneratorEtaAssocToGeantPrim_DEBUG->Fill((*p)->momentum().eta());
       }
 }


  if(fabs(atrk.charge())>0)
    {
     
      if((isTrackPrimaryIP)||((isTrackPrimaryIP==false)&&(alloweddecayparticle==true)))
	if(assoctogen)
	  primarycondition=true;
    }

  // std::cout<<"GeantTrkID:"<<trackId<<" PID_Mother:"<<thePidGenerator<<" Track PID: "<<thePidGeneratorOfThisTrk<<" Direct Primary_Condition: "<<isTrackPrimaryIP<<" Primary_Condition: "<<primarycondition<<std::endl;

 
  return primarycondition;

  

}


void T2Ntuplizer::TracksInfoAssociationDEBUG(const edm::Event& event,
							   const std::auto_ptr<edm::PSimHitContainer>& theSimHits,
							   const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
							   const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, 
							   const HepMC::GenEvent* evt,	      
							   T1T2Track recTrack,
							   std::vector<double> &PrimaryTrkenergy_ForTheRecoHits,
							   std::vector<int> &ThePrimaryTrkIdBarcode,std::vector<int> &ThePrimaryTrkIdFromRecHits,std::vector<int> &TheSecondaryTrkIdPDGID,std::vector<int> &TheSecondaryTrkIdFromRecHits
							   )
{
 
  //Convert simhits into information about plane, strip and pad.
  T2Geometry t2geom;
  T2GeometryUtil conv;     
  /*
  if (verbosity > 1){
    std::cout << "Analysing simhits.." << std::endl;
  }
  */
  
  int hitCount = recTrack.GetHitEntries();
  
//  unsigned int nonmatchingHits=0;
//  unsigned int numprimary=0;
  

  for (int h = 0;h< hitCount;++h){
    T2Hit hit = recTrack.GetHitT2(h);
    unsigned int detUnitId = hit.GetHitDetRawId();
    t2geom.setPlane(detUnitId);
    std::cout<<"----------------Process RecHit: X-Y-DX-DY: "<<hit.GetHitX()<<" "<<hit.GetHitY()<<" "<<hit.GetHitDX()<<" "<<hit.GetHitDY()<<std::endl;
 
    //Load simulated hits and associate for each corrisponding
    // T2-pad the  information about primary/secondary which generates it
    unsigned int numprimaryTEMP=0;unsigned int numsecondaryTEMP=0;

    unsigned int nummatchSimHit=0;unsigned int nummatchSimHitSaved=0;
    unsigned int nummatchSimHitSavedPrim=0;unsigned int nummatchSimHitSavedSec=0;
    for(edm::PSimHitContainer::const_iterator simHitIt = theSimHits->begin(); simHitIt != theSimHits->end(); ++simHitIt){
     unsigned int detUnitIdSIM = simHitIt->detUnitId();     

      if(detUnitIdSIM==detUnitId){
	

	/*
	T2GeometryUtil::T2DetInfo detInfo = conv.GetT2Info(detUnitId);
	t2geom.setPlane(detUnitId); //getNearestPadRo and col does not use the det info. Are refered to one particular plane.
    
	Local3DPoint hitPos = simHitIt->localPosition();    
	Local3DPoint corPos = CalculateCorrectXYZPos(Local3DPoint(hitPos.x(),hitPos.y(),hitPos.z()),detUnitId);

	//padRow = t2geom.getNearestPadRo(&corPos);
	//padCol = t2geom.getNearestPadCo(&corPos); 
	double radiusGhit=sqrt(hitPos.x()*hitPos.x()+hitPos.y()*hitPos.y()); 
	double phiGHit=atan(fabs(hitPos.y())/fabs(hitPos.x()));//atan2(hitPos.y(),hitPos.x());

	
	phiGHit=phiGHit*180./3.14159;


	if((corPos.y()<0)&&(corPos.x()>0))
	  phiGHit=360.0-phiGHit;
	
	if((corPos.y()>0)&&(corPos.x()<0))
	  phiGHit=180.0-phiGHit;

	if((corPos.y()<0)&&(corPos.x()<0))
	  phiGHit=phiGHit+180.;

	padRow = t2geom.getNearestPadRow2_(&corPos,detUnitId);
	padCol =t2geom.getNearestPadCol_(radiusGhit,phiGHit,detUnitId)
	*/


	
	Local3DPoint hitPos = simHitIt->localPosition();    
	Local3DPoint corPos = CalculateCorrectXYZPos(Local3DPoint(hitPos.x(),hitPos.y(),hitPos.z()),detUnitId);
    
	//int padRow = t2geom.getNearestPadRow(&corPos);
	//int padCol = t2geom.getNearestPadCol(&corPos);
	
	if(hit.GetHitNumPad()>0)
	  if(fabs(corPos.x()-hit.GetHitX())<3.0*hit.GetHitDX())
	    if(fabs(corPos.y()-hit.GetHitY())<3.0*hit.GetHitDY())
	      {
		nummatchSimHit++;
		unsigned int trackIdFromHit = simHitIt->trackId();
		bool isHitPrimary=false; SimTrack Atrack;
		int motherbarcode=99; 
		bool trkstored=false;
		if(getTrack(trackIdFromHit,theSimTracks,Atrack) == true)
		  {
		    nummatchSimHitSaved++;
		    isHitPrimary=PrimaryTrackCondition(Atrack,theSimVertices,theSimTracks,evt,motherbarcode);
		    PrimaryTrkenergy_ForTheRecoHits.push_back(Atrack.momentum().e());
		    trkstored=true;
		    std::cout<<"Hit PdgID:"<<simHitIt->particleType()<<" HAVE the track-ID:"<<simHitIt->trackId()<<" stored. ProcessType:"<<simHitIt->processType()<<" XY: "<<corPos.x()<<" "<<corPos.y()<<std::endl;	
		  }else
		  std::cout<<"Hit PdgID:"<<simHitIt->particleType()<<" doesn't have the track-ID:"<<simHitIt->trackId()<<" stored. ProcessType:"<<simHitIt->processType()<<" XY: "<<corPos.x()<<" "<<corPos.y()<<std::endl;	
		
		if(isHitPrimary){
		  nummatchSimHitSavedPrim++;
		  numprimaryTEMP++;		  
		  ThePrimaryTrkIdBarcode.push_back(Atrack.genpartIndex());//trackIdFromHit
		  ThePrimaryTrkIdFromRecHits.push_back(trackIdFromHit);
		}
		else{
		  
		  numsecondaryTEMP++; 
		  if(trkstored){
		    nummatchSimHitSavedSec++;
		    TheSecondaryTrkIdPDGID.push_back(Atrack.type()); 
		    TheSecondaryTrkIdFromRecHits.push_back(trackIdFromHit);
		  }
		  
		}
		    
	      }
      }

    }//end SimHitLoop
    std::cout<<"---->ReportMatch:"<<nummatchSimHit<<" "<<nummatchSimHitSaved<<" "<<nummatchSimHitSavedSec<<" "<<nummatchSimHitSavedPrim<<std::endl;
  }


     //      return trkinformation;

}





 

std::vector<int> T2Ntuplizer::RecotracksInfoFromSimu_BIS(const edm::Event& event,
							   const std::auto_ptr<edm::PSimHitContainer>& theSimHits,
							   const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
							   const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, 
							   const HepMC::GenEvent* evt,	      
							   T1T2Track recTrack,
							   std::vector<double> &PrimaryTrkenergy_ForTheRecoHits,
							   std::vector<double> &PrimaryTrkpt_ForTheRecoHits,
							   std::vector<int> &ThePrimaryTrkIdBarcode,std::vector<int> &ThePrimaryTrkIdFromRecHits,std::vector<int> &TheSecondaryTrkIdPDGID,std::vector<int> &TheSecondaryTrkIdFromRecHits
							   )
{
 
  //Convert simhits into information about plane, strip and pad.
  T2Geometry t2geom;
  T2GeometryUtil conv;     
  /*
  if (verbosity > 1){
    std::cout << "Analysing simhits.." << std::endl;
  }
  */
  
  int hitCount = recTrack.GetHitEntries();
  
  for (int h = 0;h< hitCount;++h){
    T2Hit hit = recTrack.GetHitT2(h);
    unsigned int detUnitId = hit.GetHitDetRawId();
    t2geom.setPlane(detUnitId);
    
 
    //Load simulated hits and associate for each corrisponding
    // T2-pad the  information about primary/secondary which generates it
    unsigned int numprimaryTEMP=0;unsigned int numsecondaryTEMP=0;
    for(edm::PSimHitContainer::const_iterator simHitIt = theSimHits->begin(); simHitIt != theSimHits->end(); ++simHitIt){
     unsigned int detUnitIdSIM = simHitIt->detUnitId();     

      if(detUnitIdSIM==detUnitId){
	

	/*
	T2GeometryUtil::T2DetInfo detInfo = conv.GetT2Info(detUnitId);
	t2geom.setPlane(detUnitId); //getNearestPadRo and col does not use the det info. Are refered to one particular plane.
    
	Local3DPoint hitPos = simHitIt->localPosition();    
	Local3DPoint corPos = CalculateCorrectXYZPos(Local3DPoint(hitPos.x(),hitPos.y(),hitPos.z()),detUnitId);

	//padRow = t2geom.getNearestPadRo(&corPos);
	//padCol = t2geom.getNearestPadCo(&corPos); 
	double radiusGhit=sqrt(hitPos.x()*hitPos.x()+hitPos.y()*hitPos.y()); 
	double phiGHit=atan(fabs(hitPos.y())/fabs(hitPos.x()));//atan2(hitPos.y(),hitPos.x());

	
	phiGHit=phiGHit*180./3.14159;


	if((corPos.y()<0)&&(corPos.x()>0))
	  phiGHit=360.0-phiGHit;
	
	if((corPos.y()>0)&&(corPos.x()<0))
	  phiGHit=180.0-phiGHit;

	if((corPos.y()<0)&&(corPos.x()<0))
	  phiGHit=phiGHit+180.;

	padRow = t2geom.getNearestPadRow2_(&corPos,detUnitId);
	padCol =t2geom.getNearestPadCol_(radiusGhit,phiGHit,detUnitId)
	*/


	
	Local3DPoint hitPos = simHitIt->localPosition();    
	Local3DPoint corPos = CalculateCorrectXYZPos(Local3DPoint(hitPos.x(),hitPos.y(),hitPos.z()),detUnitId);
    
	//int padRow = t2geom.getNearestPadRow(&corPos);
	//int padCol = t2geom.getNearestPadCol(&corPos);

	if(hit.GetHitNumPad()>0)
	  if(fabs(corPos.x()-hit.GetHitX())<3.0*hit.GetHitDX())
	    if(fabs(corPos.y()-hit.GetHitY())<3.0*hit.GetHitDY())
	      {

		unsigned int trackIdFromHit = simHitIt->trackId();
		bool isHitPrimary=false; SimTrack Atrack;
		int motherbarcode=99; 
		bool trkstored=false;
		if(getTrack(trackIdFromHit,theSimTracks,Atrack) == true)
		  {
		    isHitPrimary=PrimaryTrackCondition(Atrack,theSimVertices,theSimTracks,evt,motherbarcode);
		    PrimaryTrkenergy_ForTheRecoHits.push_back(Atrack.momentum().e());
		    PrimaryTrkpt_ForTheRecoHits.push_back(Atrack.momentum().pt());
		    trkstored=true;
		  }
	
		
		if(isHitPrimary){
		  numprimaryTEMP++;		  
		  ThePrimaryTrkIdBarcode.push_back(Atrack.genpartIndex());//trackIdFromHit
		  ThePrimaryTrkIdFromRecHits.push_back(trackIdFromHit);
		}
		else{
		  numsecondaryTEMP++; 
		  if(trkstored){
		    TheSecondaryTrkIdPDGID.push_back(Atrack.type());//(Atrack.genpartIndex());//(Atrack.type()); 
		    TheSecondaryTrkIdFromRecHits.push_back(trackIdFromHit);
		  }
		  
		}
		    
	      }
      }

    }//psimhit loop end

    
      //Look If you have found unique prim/sec association:
    if((numprimaryTEMP==0)&&(numsecondaryTEMP!=0)){}
//      numsecondary++;
    
    if((numprimaryTEMP!=0)&&(numsecondaryTEMP==0)){}
//      numprimary++;

    //AND WHAT IF BOTH ARE >0?
   }

     std::vector<int> trkinformation;
//     trkinformation.push_back(numdifferentsimutrks);
//     trkinformation.push_back(numprimary);
//     trkinformation.push_back(numsecondary);
//     trkinformation.push_back(nonmatchingHits);



      return trkinformation;

}
 

/*
  Reads all simulated tracks into theSimTracks.
  Returns number of tracks read.
 */ 
unsigned int
T2Ntuplizer::loadSimulatedTracks(const edm::Event& iEvent,
					  std::auto_ptr<edm::SimTrackContainer>& theSimTracks)
{
 

  edm::Handle<edm::SimTrackContainer> simTrackCollection;
  iEvent.getByLabel(simTrackContainerLabel, simTrackCollection);  
  theSimTracks.get()->clear();

  for(edm::SimTrackContainer::const_iterator ittrack = simTrackCollection->begin();ittrack != simTrackCollection->end(); ++ittrack){ 
    SimTrack thetrack=(*ittrack);
    theSimTracks->push_back(thetrack);
  }


  return theSimTracks->size();
}

/*
  Reads all simulated vertices into theSimVertices.
  Returns number of vertices read.
 */
unsigned int
T2Ntuplizer::loadSimulatedVertices(const edm::Event& iEvent,
					    std::auto_ptr<edm::SimVertexContainer>& theSimVertices)
{
  //if (verbosity > 1){
  //std::cout << "Loading all simvertices from root file..." << std::endl;
  // }

  edm::Handle<edm::SimVertexContainer> simVertexCollection;
  iEvent.getByLabel(simVertexContainerLabel, simVertexCollection);  
  theSimVertices.get()->clear();

  for(edm::SimVertexContainer::const_iterator itvertex = simVertexCollection->begin();itvertex != simVertexCollection->end(); ++itvertex){ 
    SimVertex thevertex=(*itvertex);   
    theSimVertices->push_back(thevertex);
  }

  //if (verbosity > 1){
  //std::cout << "Done reading simvertices." << std::endl << std::endl;    
  //}

  return theSimVertices->size(); //bbb

}

/*
  Reads all simulated hits into theSimHits.
  Returns number of simhits read.
 */
unsigned int T2Ntuplizer::loadSimulatedHits(const edm::Event& iEvent, 
					       const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
					       const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
					       std::auto_ptr<edm::PSimHitContainer>& theSimHits)
{
  
  edm::Handle<edm::PSimHitContainer> psimHitCollection;
  iEvent.getByLabel("g4SimHits","TotemHitsT2Gem",psimHitCollection);
  theSimHits.get()->clear();
  for(edm::PSimHitContainer::const_iterator ithit = psimHitCollection->begin();ithit != psimHitCollection->end(); ++ithit){ 
    theSimHits->push_back(*ithit); 
  }
  return theSimHits.get()->size();
}




/*
  Returns track corresponding to trackId 
  Returns true if found, false otherwise.
 */
bool
T2Ntuplizer::getTrack(const /*unsigned*/ int trackId,  const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, SimTrack &track)
{
  if(trackId < 1){
    return false;
  }

  //Guesss
  if (trackId <= (int) theSimTracks.get()->size()){
    track = (*theSimTracks.get())[trackId-1];
    if ((int) track.trackId() ==  trackId){
      return true;
    }
  }

  //Iterate all
  for(edm::SimTrackContainer::const_iterator ittrack = theSimTracks->begin();ittrack != theSimTracks->end(); ++ittrack){     
    if ((int) ittrack->trackId() ==  trackId){
      track  = (*ittrack);
      return true;
    }
  }
  
  return false;  
}

/*
  Returns vertex corresponding to vertexId 
  Returns true if found, false otherwise.
*/

bool
T2Ntuplizer::getVertex(const /*unsigned*/ int vertexId,  const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, SimVertex &vertex)
{
  if ((vertexId < (int) theSimVertices.get()->size()) && vertexId >=0){
    vertex = (*theSimVertices.get())[vertexId];
    return true;
  }
  return false;  
}


bool 
T2Ntuplizer::isTrackPrimary(const /*unsigned*/ int trackId,
					 const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
					 const std::auto_ptr<edm::SimVertexContainer>& theSimVertices)
{
  SimTrack track;
  if (getTrack(trackId,theSimTracks,track) == true){
    SimVertex vertex;    
    if (getVertex(track.vertIndex(), theSimVertices,vertex) == true){
      if (vertex.parentIndex() < 0){
	return true;
      }
    }
  }
  return false;
}



/*
  Calculates correct XYZ pos from simhit pos
 */
Local3DPoint
T2Ntuplizer::CalculateCorrectXYZPos(const Local3DPoint &pos, const unsigned int detUnitId)
{
  T2GeometryUtil conv;     
  T2GeometryUtil::T2DetInfo detInfo = conv.GetT2Info(detUnitId);  

  double corr_x_pos = pos.x();
  double corr_y_pos = pos.y();
  double corr_z_pos = pos.z();

  if (detInfo.ht == 1 && detInfo.plside == 1){
    corr_x_pos = -corr_x_pos;      
     if (detInfo.arm == 0){
      corr_y_pos = -corr_y_pos;       
       }
  } else if (detInfo.ht == 1 && detInfo.plside == 0){
    corr_x_pos = -corr_x_pos;
    if (detInfo.arm == 1){
      corr_y_pos = -corr_y_pos;       
        }
  } else if (detInfo.ht == 0 && detInfo.plside == 1){  //REFERENCE
     if (detInfo.arm == 1){
     corr_y_pos = -corr_y_pos;       
      }
  } else if (detInfo.ht == 0 && detInfo.plside == 0){
    if (detInfo.arm == 0){
      corr_y_pos = -corr_y_pos;       
     }
  }    
  
  return Local3DPoint(corr_x_pos,corr_y_pos,corr_z_pos);
}


//
// member functions
//
// ------------ method called to for each event  ------------
void T2Ntuplizer::FillEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace HepMC;

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  /* LOADING OF ALL THE RECORDS FROM THE EVENT */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  /*
  std::cout<<"A"<<std::endl;
  edm::ESHandle<HepPDT::ParticleDataTable> fPDGTable ;
  iSetup.getData( fPDGTable );
  std::cout<<"B"<<std::endl;
  */
  //121   edm::ESHandle < HepPDT::ParticleDataTable > pdt;
  //122   //  edm::ESHandle < DefaultConfig::ParticleDataTable > pdt;
  //123   es.getData(pdt);

    /* :::::::::::::TakeDigi::::::::::::*/

  edm::Handle<T2PadDigiCollection> t2paddigicoll;
  iEvent.getByLabel(t2PadDigiCollectionLabel, t2paddigicoll);
  const T2PadDigiCollection* PadDigiptr;
  PadDigiptr= t2paddigicoll.product();
  edm::Handle<T2StripDigiCollection> t2stripdigicoll;
  iEvent.getByLabel(t2StripDigiCollectionLabel, t2stripdigicoll);
  t2stripdigicoll.product(); 
DigiContainerIterator<T2DetId, T2PadDigi> itp;
  DigiContainerIterator<T2DetId, T2StripDigi> its;

  /* :::::::::::::Take The Clusters::::::::::::*/

  //std::cout << " check --------- A " <<  std::endl;
  Handle<T2StripClusterCollection> t2strclcoll;
  iEvent.getByLabel(CluLabel,"T2StripClusters",t2strclcoll);
  Handle<T2PadClusterCollection> t2padclcoll;
  iEvent.getByLabel(CluLabel,"T2PadClusters",t2padclcoll);

/*::::::Take  T2  Hits::::::*/
  edm::Handle<T2HitCollection> t2hitcoll;
  iEvent.getByLabel("T2Hits","T2Hits",t2hitcoll);
  
   
  std::auto_ptr<edm::SimTrackContainer> theSimTracks (new edm::SimTrackContainer()); 
  std::auto_ptr<edm::SimVertexContainer> theSimVertices (new edm::SimVertexContainer());
  std::auto_ptr<edm::PSimHitContainer> theSimHits (new edm::PSimHitContainer());
  Handle<HepMCProduct> EvtHandle ;
  
  //iEvent.getByLabel("source", EvtHandle ) ;

  if(LoadPrimaryGeantInformation){   
    //LOAD SIMULATED Tracks
    loadSimulatedTracks(iEvent,theSimTracks);
    //LOAD SIMULATED VERTICES
    loadSimulatedVertices(iEvent,theSimVertices);
    //LOAD SIMULATED HITS
    loadSimulatedHits(iEvent,theSimVertices,theSimTracks,theSimHits);
    
    iEvent.getByLabel(HepMCProductLabel, EvtHandle ) ;
  }
   std::vector<double> PrimaryTrkenergy_ForTheRecoHits_Bis;
   std::vector<double> PrimaryTrkpt_ForTheRecoHits_Bis;
   std::vector<int> trkinfo_bis;
   std::vector<int> ThePrimaryTrkIdBarcode;
   std::vector<int> ThePrimaryTrkIdFromRecHits;
   std::vector<int> TheSecondaryTrkIdPDGID;
   std::vector<int> TheSecondaryTrkIdFromRecHits;
  //::::::Take particle generated with Gun::::::
  /*
  Handle<HepMCProduct> EvtHandle ;
  iEvent.getByLabel(HepMCProductLabel, EvtHandle ) ;
  //iEvent.getByLabel("source", EvtHandle ) ;
  const GenEvent* evt = EvtHandle->GetEvent();

 // ::::::Take Geant local hit on T2::::::
  Handle<CrossingFrame<PSimHit> > cFrame;
  iEvent.getByLabel("mix", "g4SimHitsTotemHitsT2Gem", cFrame);
  // iEvent.getByLabel("mix", "TotemHitsT2Gem", cFrame);
  // get hits from G4Sim
  const string nameOfHits("TotemHitsT2Gem");
  auto_ptr<MixCollection<PSimHit> >  T2SimHits( new MixCollection<PSimHit>( cFrame.product() ) );

  // map hits to T2 plane
  map<int, PSimHitContainer> hitMap;
  for(MixCollection<PSimHit>::MixItr hitItr = T2SimHits->begin(); hitItr != T2SimHits->end(); ++hitItr) {
    hitMap[hitItr->detUnitId()].push_back(*hitItr);
  }
  */
  /*
  //::::::Take  T2  Hits::::::
  Handle<T2HitCollection> t2hitcoll;

  iEvent.getByLabel(HitLabel,"T2Hits",t2hitcoll);
  //::::::Take  T2  Roads::::::

  Handle<T2RoadCollection> t2roadcoll;
  iEvent.getByLabel(RoadLabel,"T2RoadColl",t2roadcoll);
  */
  //:::::: Take T2 tracks ::::::
  Handle<T1T2TrackCollection> trackCollection;
  iEvent.getByLabel(TrackLabel,"T2TrackColl",trackCollection);
  
 
  //std::cout<<"New Event"<<std::endl;

  std::map<int,std::vector<int> > Pi0Map;//Containing left:barcode of pion, right geanttrk-id
  std::map<int,std::vector<double> > Pi0MapZImp;
  std::map<int,std::vector<int> > Pi0MapTrkIndex;
  
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
 t2obj.TrkNature.clear();      // 1=primary 2=secondary
 t2obj.Trk_NumHitPrimary.clear();      // 1=primary 2=secondary


 t2obj.PrimaryGeantTrkEta.clear(); 
 t2obj.PrimaryGeantTrkQuarter.clear(); 
 t2obj.PrimGeantTrkPT.clear(); 
 t2obj.TrkPrimaryGeneratorEta.clear(); 
 t2obj.TrkPrimaryGeantEta.clear(); 
 t2obj.TrkPrimaryEnergy.clear();  // E primary;  E=-1 associated to secondary
 t2obj.TrkPrimaryPt.clear();
 t2obj.TrkPrimaryBarcode.clear();  

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

 t2obj.TrkEta_XY.clear();  
 
 t2obj.TrkZmin_XY.clear();   
 //t2obj.TrkZ0_XY.clear();  
 //t2obj.TrkR0_XY.clear();   
 t2obj.TrkRmin_XY.clear();   
 t2obj.TrkEta2.clear();   
 
 t2obj.TrkChiProb.clear();   
 t2obj.Trknumpadonly.clear();   
 t2obj.Trknumstriponly.clear();   
 t2obj.Trknumhitall.clear();   
 t2obj.Trknumhitacl1.clear();   


 t2obj.TrkNumHitInH0.clear();  
 t2obj.TrkNumHitInH1.clear();  
 t2obj.TrkNumHitInH2.clear();  
 t2obj.TrkNumHitInH3.clear();  
 t2obj.GeneratorStableCHPartileEta.clear();
 t2obj.GeneratorStableCHPartilePhi.clear();

 t2obj.NumPadCluH0=0;
 t2obj.NumPadCluH1=0;
 t2obj.NumPadCluH2=0;
 t2obj.NumPadCluH3=0;

 t2obj.TrkSecondaryRecoIndex.clear();
 t2obj.TrkSecondaryPDG.clear();
 t2obj.TrkSecondaryMotherPDG.clear();
 // t2obj.TrkSecondaryGeantRecTrk.clear();
 t2obj.GammaVtxSize.clear();
 t2obj.GammaVtxSizeZImp.clear();
 t2obj.GammaVtxSizeTrkIndex.clear();

 t2obj.GeneratorAllPartileEta.clear();
 t2obj.GeneratorAllPartilePhi.clear();

 t2obj.GeneratorAllPartilePDG.clear();
 t2obj.GeneratorAllPartileSTATUS.clear(); 
 t2obj.GeneratorAllPartileE.clear();
 t2obj.GeneratorAllPartilePt.clear();
 t2obj.GeneratorAllPartileCHARGE.clear();
 t2obj.GeneratorAllPartileAncestorPID.clear(); 
 t2obj.GeneratorAllPartileAncestorSize.clear();
 t2obj.GeneratorAllPartileCharge1.clear();
 t2obj.GeneratorAllPartileCharge2.clear();
 std::map<int,double> barcodestableChp;


 t2obj.PrimaryGeantTrkEta.clear(); 
 t2obj.PrimaryGeantTrkQuarter.clear(); 
 t2obj.PrimaryGeantTrkBarcode.clear(); 

 std::vector<double> PrimGeantTrkEta2InT2; 
 std::vector<int> PrimGeantTrkHemiInT2;
 std::vector<double> PrimGeantTrkPTInT2; 
 std::vector<int> PrimGeantTrkBarcodeInT2;
 iEvent.getByLabel("generator", EvtHandle ) ;
 const HepMC::GenEvent* evt = EvtHandle->GetEvent();
 
 GeantTrkAnalyzer(theSimHits, theSimVertices, theSimTracks,  evt,PrimGeantTrkEta2InT2, PrimGeantTrkHemiInT2,PrimGeantTrkPTInT2,PrimGeantTrkBarcodeInT2,0.0);

 for(unsigned int yy=0; yy<PrimGeantTrkHemiInT2.size();yy++)
   {
     t2obj.PrimaryGeantTrkEta.push_back(PrimGeantTrkEta2InT2.at(yy));
     t2obj.PrimaryGeantTrkQuarter.push_back(PrimGeantTrkHemiInT2.at(yy));
     t2obj.PrimGeantTrkPT.push_back(PrimGeantTrkPTInT2.at(yy));
     t2obj.PrimaryGeantTrkBarcode.push_back(PrimGeantTrkBarcodeInT2.at(yy));
   }
 

 
 if(LoadPrimaryGeantInformation){ 
  
   for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){	       
     int charge=PartCharge((*p)->pdg_id());  //OBSOLETE WAY (*p)->charge();//
     charge=abs(charge);

     t2obj.GeneratorAllPartileCharge1.push_back(charge);
     /*
      const HepPDT::ParticleData* PData = fPDGTable->particle(HepPDT::ParticleID((*p)->pdg_id()));
      if(PData==0) {
         std::cout << "Unknown id for charge= " << (*p)->pdg_id() << '\n';
      }
      else
	charge = PData->charge();
     */

     // t2obj.GeneratorAllPartileCharge2.push_back(charge);
     t2obj.GeneratorAllPartileEta.push_back((*p)->momentum().eta());
     t2obj.GeneratorAllPartilePhi.push_back((*p)->momentum().phi());
     t2obj.GeneratorAllPartileE.push_back((*p)->momentum().e());
     t2obj.GeneratorAllPartilePt.push_back((*p)->momentum().perp());
     t2obj.GeneratorAllPartilePDG.push_back((*p)->pdg_id());
     t2obj.GeneratorAllPartileSTATUS.push_back((*p)->status());
     t2obj.GeneratorAllPartileCHARGE.push_back(charge);
     t2obj.GeneratorAllPartileBarcode.push_back((*p)->barcode());
     

     if((abs(charge)>0)&&((*p)->status()==1))
       {	      	    	    
	 double etag=(*p)->momentum().eta();
		
	 t2obj.GeneratorStableCHPartileEta.push_back(etag);
	 double phig=(*p)->momentum().phi();
	 t2obj.GeneratorStableCHPartilePhi.push_back(phig);
	 barcodestableChp.insert(pair<int,int>((*p)->barcode(),etag));
	 
       }

     //  std::cout<<"HA"<<std::endl;
     // HepMC::GenVertex* inVertex = (*p)->production_vertex();
    
     //vector< GenParticle * > parents;
      
       //if (id == 11 || id == -11)  
	
	  // single primary electrons or electrons from Zs or Ws
	  HepMC::GenParticle* mother = 0;
	  bool present=false;
	  int ancestorsize=0;

	  if ( (*p)->production_vertex() )  {
	    if ( (*p)->production_vertex()->particles_begin(HepMC::parents) != (*p)->production_vertex()->particles_end(HepMC::parents)){
	      mother = *((*p)->production_vertex()->particles_begin(HepMC::parents));
	      present=true;
	      for ( HepMC::GenVertex::particle_iterator moth = (*p)->production_vertex()->particles_begin(HepMC::parents);moth != (*p)->production_vertex()->particles_end(HepMC::parents);  ++moth ){
		ancestorsize++;
	      }
	    }
	  }
	
	

	int OlderMotherParticle_PId=0;
	
	if(present){
	  //( (*p)->production_vertex()->particles_begin(HepMC::parents) - (*p)->production_vertex()->particles_end(HepMC::parents));
	  OlderMotherParticle_PId=mother->pdg_id();
	}



     // for(GenVertex::particles_in_const_iterator iter = inVertex->particles_in_const_begin();iter != inVertex->particles_in_const_end();iter++)
       
     //std::cout<<(*iter)->pdg_id()<<std::endl;
     
     //std::cout<<"HC!"<<parents.size()<<std::endl;
     

     /*
     std::cout<<"HA"<<std::endl;
     if ( (*p)->production_vertex() )
       std::cout<<"HAc!"<<std::endl;
     */
       
       //HepMC::GenVertex::particle_iterator ancestor;
       // ancestor = (*p)->production_vertex()->particles_begin(HepMC::parents);
    
     //std::cout<<"HB"<<std::endl;
     /*
     if ( (*p)->production_vertex() )
       for ( HepMC::GenVertex::particle_iterator mother = (*p)->production_vertex()->particles_begin(HepMC::parents);mother != (*p)->production_vertex()->particles_end(HepMC::parents);  ++mother ) {
	 std::cout << "\t";
	 (*mother)->print();
       }
     */
     
     /*
     for(;ancestor != (*p)->production_vertex()->particles_end(HepMC::parents); ++ancestor) {
       // std::cout<<"  "<<((*ancestor)->pdg_id());
       //motherbarcode to return is saved.
       //int motherbarcode=(*ancestor)->barcode();
       OlderMotherParticle_PId=(*ancestor)->pdg_id();
       ancestorsize++;       
     } 
     */

     
     t2obj.GeneratorAllPartileAncestorSize.push_back(ancestorsize);
     t2obj.GeneratorAllPartileAncestorPID.push_back(OlderMotherParticle_PId);
     

     /*
      for(;ancestor != (*p)->production_vertex()->particles_end(HepMC::parents); ++ancestor) {            
       OlderMotherParticle_Id=(*ancestor)->pdg_id();
       ancestorsize++;       
       } 
     */

    
   }
 }

 T1T2TrackCollection::const_iterator TrkCit;
 T2GeometryUtil::T2DetInfo planeinfo;
 T2GeometryUtil conv;


 
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

//TO DECOMMENT   		   
 /*
     for(its= StripDigiptr->begin(); its!=StripDigiptr->end(); ++its)
       {
	 T2DetId mydet=(*its).first; 
	 int detnumb = mydet.arm()*20 + mydet.halfTelescope()*10 + mydet.plane()*2 + mydet.planeSide();

	 for(std::vector<T2StripDigi>::const_iterator itstrip =(*its).second.first; itstrip !=(*its).second.second; itstrip++)
	   {
	     t2obj.Strip_row.push_back((*itstrip).getRow());
	     t2obj.Strip_col.push_back((*itstrip).getCol());
	     t2obj.Strip_det.push_back(detnumb);
	   }
       }

 */

 
  	  	  
	 
	 
	 
	 
     int rectrackindex=0;

     for(TrkCit=trackCollection->begin(); TrkCit!=trackCollection->end(); TrkCit++){
       //std::vector<T2Hit> hitvector;	
       unsigned int class1Hitcounter=0;
       unsigned int numh0=0;
       unsigned int numh1=0;
       unsigned int numh2=0;
       unsigned int numh3=0;
       double theEnergy=-1.;double thePt=-1.;
       int barcode=-1;
       double GenEta=0.; double GeantEta=0.;
      
       PrimaryTrkenergy_ForTheRecoHits_Bis.clear();
       trkinfo_bis.clear();
       ThePrimaryTrkIdBarcode.clear();
       ThePrimaryTrkIdFromRecHits.clear();
       TheSecondaryTrkIdPDGID.clear();
       TheSecondaryTrkIdFromRecHits.clear();
       PrimaryTrkpt_ForTheRecoHits_Bis.clear();

       unsigned int nature=0;
       int Trk_NumHitPrimary_Tuput=-1;
       if(LoadPrimaryGeantInformation){   
	 
	 ThePrimaryTrkIdBarcode.clear();
	 const GenEvent* evt = EvtHandle->GetEvent();
	 /*
	 for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){    
	   std::cout<<"Barcode:"<<((*p)->barcode())<<" PdgID:"<<(*p)->pdg_id()<<" status_code:"<<(*p)->status()<<std::endl;    
	   }
*/

	 trkinfo_bis=RecotracksInfoFromSimu_BIS(iEvent,theSimHits,theSimTracks,theSimVertices,evt,(*TrkCit),PrimaryTrkenergy_ForTheRecoHits_Bis,PrimaryTrkpt_ForTheRecoHits_Bis,ThePrimaryTrkIdBarcode,ThePrimaryTrkIdFromRecHits,TheSecondaryTrkIdPDGID,TheSecondaryTrkIdFromRecHits);

	 
	 Trk_NumHitPrimary_Tuput=trkinfo_bis[1];
	 if(trkinfo_bis[2]>trkinfo_bis[1])
	   {
	     nature=2;
	     /*
	     if(TheSecondaryTrkIdFromRecHits.size()>=3){
	       T1T2Track GeantTrk; 
	       bool GtrackReco=CalculateGeantEtaFromTrkId(ThePrimaryTrkIdFromRecHits,theSimHits,theSimTracks,GeantTrk);
	       double thetaHitAverage=0.;     
	       double hemis=0.; 
	       unsigned int trknumhit=GeantTrk.GetHitEntries();
	       for(unsigned int r=0;r<trknumhit;r++)
		 {
		   thetaHitAverage += fabs(GeantTrk.GetHitT2(r).GetHitR()/GeantTrk.GetHitT2(r).GetHitZ());
		   hemis=fabs(GeantTrk.GetHitT2(r).GetHitZ())/GeantTrk.GetHitT2(r).GetHitZ();
		   //std::cout<<t2trackVectorBothSide.at(u).GetHitT2(r).GetHitZ()<<" "<<t2trackVectorBothSide.at(u).GetHitT2(r).GetHitR()<<" "<<t2trackVectorBothSide.at(u).GetHitT2(r).GetHitPhi()<<"  Pad-StripCLS:"<<t2trackVectorBothSide.at(u).GetHitT2(r).GetHitNumPad()<<" "<<t2trackVectorBothSide.at(u).GetHitT2(r).GetHitNumStrip()<<"    |    ";
		 }
	       thetaHitAverage=thetaHitAverage/trknumhit;  
	       GeantEta= (hemis*(-1.0)*log(thetaHitAverage/2.0)); //low angle approximation
	       }
	     */
	   }
	 else{
	   if(PrimaryTrkenergy_ForTheRecoHits_Bis.size()>0)
	     {
	       nature=1;
	       theEnergy=MostProbableTrkEnergy(PrimaryTrkenergy_ForTheRecoHits_Bis);
	       thePt=MostProbableTrkEnergy(PrimaryTrkpt_ForTheRecoHits_Bis); //MostProbableTrkEnergy can be used as it is also for pt
	       GenEta=CalculateGeneratorEtaFromTrkId(ThePrimaryTrkIdBarcode,evt);	       
	       barcode=ThePrimaryTrkIdBarcode.at(0);

	       bool foundbar=false;	       
	       for(unsigned int hh=0;hh<ThePrimaryTrkIdBarcode.size();hh++)
		 if(barcodestableChp.find(ThePrimaryTrkIdBarcode.at(hh))!=barcodestableChp.end()){
		   barcode=ThePrimaryTrkIdBarcode.at(hh);
		   foundbar=true;
		 }
	       
	       if(foundbar==false){
		 std::cout<<"This primary Trk does not have barcode in HepMC. Barc:"<<barcode<<" etaGen:"<<GenEta<<" energy:"<<theEnergy<<std::endl;
	       }

	       T1T2Track GeantTrk; 
	       bool GtrackReco=CalculateGeantEtaFromTrkId(ThePrimaryTrkIdFromRecHits,theSimHits,theSimTracks,GeantTrk);
	       
	       if(GtrackReco){
		 double thetaHitAverage=0.;     
		 double hemis=0.; 
		 unsigned int trknumhit=GeantTrk.GetHitEntries();
		 for(unsigned int r=0;r<trknumhit;r++)
		   {
		     thetaHitAverage += fabs(GeantTrk.GetHitT2(r).GetHitR()/GeantTrk.GetHitT2(r).GetHitZ());
		     hemis=fabs(GeantTrk.GetHitT2(r).GetHitZ())/GeantTrk.GetHitT2(r).GetHitZ();
		     //std::cout<<t2trackVectorBothSide.at(u).GetHitT2(r).GetHitZ()<<" "<<t2trackVectorBothSide.at(u).GetHitT2(r).GetHitR()<<" "<<t2trackVectorBothSide.at(u).GetHitT2(r).GetHitPhi()<<"  Pad-StripCLS:"<<t2trackVectorBothSide.at(u).GetHitT2(r).GetHitNumPad()<<" "<<t2trackVectorBothSide.at(u).GetHitT2(r).GetHitNumStrip()<<"    |    ";
		   }
		 thetaHitAverage=thetaHitAverage/trknumhit;  
		 GeantEta= (hemis*(-1.0)*log(thetaHitAverage/2.0)); //low angle approximation
	       }
	     }
	   else{
	     nature=0;
	   }
	 }
	 	 
       }else
	 nature=9;





       double C0=((TrkCit->GetHitT2(0).GetHitX()*TrkCit->GetHitT2(0).GetHitX()+TrkCit->GetHitT2(0).GetHitY()*TrkCit->GetHitT2(0).GetHitY())/(TrkCit->GetTx()*TrkCit->GetHitT2(0).GetHitX()+TrkCit->GetTy()*TrkCit->GetHitT2(0).GetHitY()));
       double Z0impact=TrkCit->GetHitT2(0).GetHitZ()-C0;
       
       if(TheSecondaryTrkIdPDGID.size()>3){
	 int pdgsec=-1;
	 pdgsec=CaluclateCommonPDGSEC(TheSecondaryTrkIdPDGID);
	 t2obj.TrkSecondaryRecoIndex.push_back(rectrackindex);
	 t2obj.TrkSecondaryPDG.push_back(pdgsec);
	 int pdgsecmoth=-1;
	 const HepMC::GenEvent* evt = EvtHandle->GetEvent(); 
	 int motherbarcode=0;
	 pdgsecmoth=CaluclateCommonPDGSECMoth(TheSecondaryTrkIdFromRecHits,theSimVertices,theSimTracks, evt,false,motherbarcode);
	 t2obj.TrkSecondaryMotherPDG.push_back(pdgsecmoth);
	 
	 if((pdgsecmoth!=-1)&&(pdgsecmoth!=0))//It is a product of a decay particle
	   nature=fabs(pdgsecmoth);//Was 3
	 else
	   std::cout<<"WarningLOSS: TheSecondaryTrkIdPDGID>3 but pdgsecmoth="<<pdgsecmoth<<std::endl;

	 if((pdgsecmoth==22)||(pdgsecmoth==111)){//Phoj desn't have 111	    
	   nature=4; 
	   int rectrkindex=rectrackindex;
	   CaluclateCommonPDGSECMoth(TheSecondaryTrkIdFromRecHits,theSimVertices,theSimTracks, evt,false,motherbarcode);
	   std::cout<<"-> Isert Pi0Map with barcode"<<motherbarcode<<std::endl;
	   //   mymap.insert(pair<int,int>('a',100) );
	   //TracksInfoAssociationDEBUG(iEvent,theSimHits,theSimTracks,theSimVertices,evt,(*TrkCit),PrimaryTrkenergy_ForTheRecoHits_Bis,ThePrimaryTrkIdBarcode,ThePrimaryTrkIdFromRecHits,TheSecondaryTrkIdPDGID,TheSecondaryTrkIdFromRecHits);
	   map<int,std::vector<int> >::iterator it; 	  
	   it=Pi0Map.find(motherbarcode);
	   if(it==Pi0Map.end()){
	     std::vector<int> aa;
	     aa.push_back(TheSecondaryTrkIdPDGID.at(0));
	     std::vector<double> bb;
	     bb.push_back(Z0impact);
	     std::vector<int> cc;	    
	     cc.push_back(rectrkindex);

	     Pi0MapZImp.insert(pair<int,std::vector<double> >(motherbarcode,bb));
	     Pi0Map.insert(pair<int,std::vector<int> >(motherbarcode,aa));	 
	     Pi0MapTrkIndex.insert(pair<int,std::vector<int> >(motherbarcode,cc));
	   }else{
	     Pi0Map[motherbarcode].push_back(TheSecondaryTrkIdPDGID.at(0));
	     Pi0MapZImp[motherbarcode].push_back(Z0impact);
	     Pi0MapTrkIndex[motherbarcode].push_back(rectrkindex);
	   }
	    
	 }
	 if(pdgsecmoth==310)
	    nature=5;

	 if(fabs(pdgsecmoth)==3122)
	    nature=6;

	 if((fabs(pdgsecmoth)==3112)||(fabs(pdgsecmoth)==3222))
	   nature=7;

	 if((fabs(pdgsecmoth)==3312)||(fabs(pdgsecmoth)==3322))
	   nature=8;
	 //Look around line 1280
       }else
	 if(nature==2){
	  
	   if(fabs(Z0impact)<1500){
	      //const HepMC::GenEvent* evt = 
		EvtHandle->GetEvent();
	      //std::cout<<"Small Z0impact="<<Z0impact<<" secondary tracks not traced backword. #HitAssociated:"<<TheSecondaryTrkIdPDGID.size()<<" SecVsPrim:"<<trkinfo_bis[2]<<" "<<trkinfo_bis[1]<<std::endl;
	     // TracksInfoAssociationDEBUG(iEvent,theSimHits,theSimTracks,theSimVertices,evt,(*TrkCit),PrimaryTrkenergy_ForTheRecoHits_Bis,ThePrimaryTrkIdPDGID,ThePrimaryTrkIdFromRecHits,TheSecondaryTrkIdPDGID,TheSecondaryTrkIdFromRecHits);
	   }
	 }
       
       t2obj.Trk_NumHitPrimary.push_back(Trk_NumHitPrimary_Tuput);
       t2obj.TrkNature.push_back(nature);  
       
       t2obj.TrkPrimaryEnergy.push_back(theEnergy);
       t2obj.TrkPrimaryPt.push_back(thePt);

       t2obj.TrkPrimaryBarcode.push_back(barcode);
       t2obj.TrkPrimaryGeneratorEta.push_back(GenEta);

       t2obj.TrkPrimaryGeantEta.push_back(GeantEta);

	 for(unsigned int i=0; i<(*TrkCit).GetHitEntries();i++)
	   {
	     if((*TrkCit).GetHitT2(i).GetHitClass()==1){
	       class1Hitcounter++;
	        

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
	      
	       
	     }
	     //  hitvector.push_back((*TrkCit).GetHitT2(i));
	   }
	 


	  t2obj.TrkNumHitInH0.push_back(numh0);
	  t2obj.TrkNumHitInH1.push_back(numh1);
	  t2obj.TrkNumHitInH2.push_back(numh2);
	  t2obj.TrkNumHitInH3.push_back(numh3);


	  t2obj.TrkEta_XY.push_back(TrkCit->Eta());            
	  
	  t2obj.TrkZmin_XY.push_back(TrkCit->Z_at_Rmin());
	 
	  t2obj.TrkRmin_XY.push_back(TrkCit->Rmin());  
	 
	  double eta2=/*T2CutsUtil.*/EtaFromAveragesVtxCorr((*TrkCit),0.,0.,0.); //11250.
	  t2obj.TrkEta2.push_back(eta2);
	  

	  t2obj.TrkHitCounter.push_back((*TrkCit).GetHitEntries());
	  t2obj.TrkClass1HitCounter.push_back(class1Hitcounter);
	 
	  double RTrkEntry=sqrt(TrkCit->GetHitT2(0).GetHitX()*TrkCit->GetHitT2(0).GetHitX()+TrkCit->GetHitT2(0).GetHitY()*TrkCit->GetHitT2(0).GetHitY());
	  //double hitRecoEta=
	log((RTrkEntry/(TrkCit->GetHitT2(0).GetHitZ()))/2.0);
	  // std::cout<<"RecoEta2:"<<eta2<<"   Geant-Eta:"<<GeantEta<<" FirstHitRecoEta"<<hitRecoEta<<std::endl;
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



	 //R0xy=sqrt(((*TC_iter).bx_)*((*TC_iter).bx_)+((*TC_iter).by_)*((*TC_iter).by_));
	 //R0rz=sqrt(((*TC_iter).bx_)*((*TC_iter).bx_)+((*TC_iter).by_)*((*TC_iter).by_)); 
	 //Z0xy=(*TC_iter).Z_at_Rmin();
	 // Z0rz=(*TC_iter).Z_at_Rmin();
	 double ProbChi2X_xy=TMath::Prob((*TrkCit).ChiSquaredX(),((*TrkCit).GetHitEntries()-2));
	 double ProbChi2Y_xy=TMath::Prob((*TrkCit).ChiSquaredY(),((*TrkCit).GetHitEntries()-2));
	 //ProbChi2R_rz=TMath::Prob((*TC_iter).ChiSquaredR() , ((*TC_iter).GetHitEntries()-2));
	 //ProbChi2Phi_rz=TMath::Prob((*TC_iter).ChiSquaredPhi(),((*TC_iter).GetHitEntries()-1));
	 
	 t2obj.TrkChi2YProb.push_back(ProbChi2Y_xy);
	 t2obj.TrkChi2XProb.push_back(ProbChi2X_xy);
	   
	 //To be uncommented for the old (T2TrackProducer2) algorithm.
	 /*
	 T1T2Track trk2rz=T2CutsUtil.TrackFromHits(false,hitvector);//RZFIT

	 t2obj.TrkZ0_RZFit.push_back(trk2rz.Z_at_Rmin());
	 t2obj.TrkPhi_RZFit.push_back(trk2rz.Phi()*180/3.14159); 
	 t2obj.TrkEta_RZFit.push_back(trk2rz.Eta());
	 t2obj.TrkThetaR_RZFit.push_back(trk2rz.GetTy());
	 t2obj.TrkBX_RZFit.push_back(trk2rz.bx_);
	 t2obj.TrkBY_RZFit.push_back(trk2rz.by_);
	 */

	 unsigned int numhitall=(*TrkCit)._numHit_t2;
	 unsigned int numpadonly=(*TrkCit)._numPadOnlyHit_t2;
	 unsigned int numstriponly=(*TrkCit)._numStripOnlyHit_t2;  	 
	 unsigned int numhitacl1=(*TrkCit)._numCl1HitHit_t2;  
	 double ProbChi2_XY=TMath::Prob((*TrkCit).ChiSquared(),(numhitacl1*2-4));
	 /*
	 if(ProbChi2_XY<0.01){
	   std::cout<<"Warning! Lowchiprob:"<<std::endl;
	   for(unsigned int i=0; i<(*TrkCit).GetHitEntries();i++)
	     {
	       // if((*TrkCit).GetHitT2(i).GetHitClass()==1){
	       std::cout<<"X-Y: "<<(*TrkCit).GetHitT2(i).GetHitX()<<" "<<(*TrkCit).GetHitT2(i).GetHitY()<<" class:"<<(*TrkCit).GetHitT2(i).GetHitClass()<<"NumPad:"<<(*TrkCit).GetHitT2(i).GetHitNumPad()<<"NumStrip:"<<(*TrkCit).GetHitT2(i).GetHitNumStrip()std::endl;
		 //}
	     }
	 }
	 */
	 t2obj.TrkChiProb.push_back(ProbChi2_XY);
	 t2obj.Trknumpadonly.push_back(numpadonly);
	 t2obj.Trknumstriponly.push_back(numstriponly);
	 t2obj.Trknumhitall.push_back(numhitall);
	 t2obj.Trknumhitacl1.push_back(numhitacl1);

	 t2obj.TrkZ0_RZFit.push_back((*TrkCit)._Z_at_RminRZ);
	 t2obj.TrkPhi_RZFit.push_back((*TrkCit)._phiRZ*180/3.14159);
	 t2obj.TrkEta_RZFit.push_back((*TrkCit)._etaRZ);
	 t2obj.TrkThetaR_RZFit.push_back((*TrkCit)._TanthetaRZ);
	 t2obj.TrkBX_RZFit.push_back((*TrkCit)._bxRZ);
	 t2obj.TrkBY_RZFit.push_back((*TrkCit)._byRZ);
	 
	 //TrkEtaALL->Fill(trk2rz.Eta());

	 rectrackindex++;
     }
 
     //TO DECOMNENT!!
     /*
 for(T2HitCollection::const_iterator ithit = t2hitcoll->begin(); ithit != t2hitcoll->end(); ithit++)
    {
      
      uint32_t rawiddet=(*ithit).GetHitDetRawId();
      unsigned int type=0;
      //planeinfo=conv.GetT2Info(rawiddet);
      if(((*ithit).GetHitNumStrip()==0)&&((*ithit).GetHitNumPad()>0))
	type=0;
      if(((*ithit).GetHitNumStrip()>0)&&((*ithit).GetHitNumPad()==0))
	type=1;
      if(((*ithit).GetHitNumStrip()>0)&&((*ithit).GetHitNumPad()>0))	
	type=2;
	
      t2obj.HitR.push_back((*ithit).GetHitR());
      t2obj.HitPhi.push_back((*ithit).GetHitPhi());   
      t2obj.HitType.push_back(type);  
      t2obj.HitNumPad.push_back((*ithit).GetHitNumPad());  
      t2obj.HitNumStrip.push_back((*ithit).GetHitNumStrip()); 
    
    }
     */

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
      delete detID;
 }

  
  t2obj.NumPadCluH0=totHitAllH0; 
  t2obj.NumPadCluH1=totHitAllH1; 
  t2obj.NumPadCluH2=totHitAllH2; 
  t2obj.NumPadCluH3=totHitAllH3;

   
     //Fill the T2 ntpl.
 


  for(map<int,std::vector<int> >::iterator it=Pi0Map.begin();it!=Pi0Map.end();it++)
    {
      std::cout<<"Pi0 Barcode: "<<it->first<<"  size:"<<it->second.size()<<std::endl;
    }


  
 for(map<int,std::vector<double> >::iterator it2=Pi0MapZImp.begin();it2!=Pi0MapZImp.end();it2++)
    {
      t2obj.GammaVtxSizeZImp.push_back(10000.);//Vtx separator
      t2obj.GammaVtxSizeTrkIndex.push_back(10000);
      std::vector<double> theZ0s=it2->second;
      std::vector<int> theTrkIndex=Pi0MapTrkIndex[it2->first];

       t2obj.GammaVtxSize.push_back(theZ0s.size());
      //      if(theZ0s.size()>1)
      for(unsigned int u=0;u<theZ0s.size();u++){
	std::cout<<"Pi0 Barcode: "<<it2->first<<"  ZImpact ("<<u<<"):"<<theZ0s.at(u)<<std::endl;
	//if(theZ0s.size()>2)
	  t2obj.GammaVtxSizeZImp.push_back(theZ0s.at(u));
	  t2obj.GammaVtxSizeTrkIndex.push_back(theTrkIndex.at(u));
      }
      
    }
  //std::cout<<"End Event"<<std::endl;
 //t2tree->Fill();
}






     





// ------------ method called once each job just before starting event loop  ------------
void T2Ntuplizer::CreateBranches(const edm::EventSetup&, TTree *tree)
{





 
 knowPart.push_back(11); knowPart.push_back(12);  knowPart.push_back(13);  knowPart.push_back(14); knowPart.push_back(15); 
  knowPart.push_back(16); knowPart.push_back(17); knowPart.push_back(18); knowPart.push_back(22); knowPart.push_back(23);
  knowPart.push_back(24); 

  knowPart.push_back(130); knowPart.push_back(310); knowPart.push_back(311); 
  
  knowPart.push_back(111); knowPart.push_back(211); knowPart.push_back(113); knowPart.push_back(213); knowPart.push_back(215); 
  knowPart.push_back(117); knowPart.push_back(217); knowPart.push_back(119); knowPart.push_back(219); knowPart.push_back(221); 
  knowPart.push_back(331); knowPart.push_back(223); knowPart.push_back(333); knowPart.push_back(225); knowPart.push_back(335); 
  knowPart.push_back(227); 
  knowPart.push_back(337); knowPart.push_back(229);  knowPart.push_back(411); knowPart.push_back(421); knowPart.push_back(413); 
  knowPart.push_back(423); knowPart.push_back(415); knowPart.push_back(425); knowPart.push_back(431); knowPart.push_back(433); 
  knowPart.push_back(435); knowPart.push_back(511); knowPart.push_back(521); knowPart.push_back(513); knowPart.push_back(523); 
  knowPart.push_back(515); knowPart.push_back(525); knowPart.push_back(531); knowPart.push_back(533); knowPart.push_back(535); 
  knowPart.push_back(541); knowPart.push_back(543); knowPart.push_back(545); knowPart.push_back(441); knowPart.push_back(443); 
  knowPart.push_back(445); knowPart.push_back(551); knowPart.push_back(553); knowPart.push_back(555); knowPart.push_back(557); 
  knowPart.push_back(2212); knowPart.push_back(2112); knowPart.push_back(2224); knowPart.push_back(2214); knowPart.push_back(2114); 
  knowPart.push_back(1114); knowPart.push_back(3122); knowPart.push_back(3222); knowPart.push_back(3212); 
  knowPart.push_back(3112); knowPart.push_back(3224); knowPart.push_back(3214); knowPart.push_back(3114); knowPart.push_back(3322); 
  knowPart.push_back(3312); knowPart.push_back(3324); knowPart.push_back(3314); knowPart.push_back(3334); knowPart.push_back(4122); 
  knowPart.push_back(4222); knowPart.push_back(4212); knowPart.push_back(4112); knowPart.push_back(4224); knowPart.push_back(4214); 
  knowPart.push_back(4114); knowPart.push_back(4232); knowPart.push_back(4132); knowPart.push_back(4322); knowPart.push_back(4312); 
  knowPart.push_back(4324); knowPart.push_back(4314); knowPart.push_back(4332); knowPart.push_back(4334); knowPart.push_back(4412); 
  knowPart.push_back(4422); knowPart.push_back(4414); knowPart.push_back(4424); knowPart.push_back(4432); knowPart.push_back(4434); 
  knowPart.push_back(4444); knowPart.push_back(5122); knowPart.push_back(5112); knowPart.push_back(5212); knowPart.push_back(5222); 
  knowPart.push_back(5114); knowPart.push_back(5214); knowPart.push_back(5224); knowPart.push_back(5132); knowPart.push_back(5232); 
  knowPart.push_back(5312); knowPart.push_back(5322); knowPart.push_back(5314); knowPart.push_back(5324); knowPart.push_back(5332); 
  knowPart.push_back(5334); knowPart.push_back(5142); knowPart.push_back(5242); knowPart.push_back(5412); knowPart.push_back(5422); 
  knowPart.push_back(5414); knowPart.push_back(5424); knowPart.push_back(5342); knowPart.push_back(5432); knowPart.push_back(5434); 
  knowPart.push_back(5442); knowPart.push_back(5444); knowPart.push_back(5512); knowPart.push_back(5522); knowPart.push_back(5514); 
  knowPart.push_back(5524); knowPart.push_back(5532); knowPart.push_back(5534); knowPart.push_back(5542); knowPart.push_back(5544); 
  knowPart.push_back(5554);



  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  
  
  //AllowedDecayToCount.push_back(11); AllowedDecayToCount.push_back(12); electr and it neutrino Excluded (not decayng particle)
  //AllowedDecayToCount.push_back(13);  AllowedDecayToCount.push_back(14); muon ans its neutrino Excluded (lifetime too big ~10-6 to decay before T2)

  //WARNING: 66% of Tau decay involve hadrons in the final states.
  AllowedDecayToCount.push_back(15); //10-12s 10-13s (Tau), will decay before T2.
  //AllowedDecayToCount.push_back(16); AllowedDecayToCount.push_back(17); AllowedDecayToCount.push_back(18); Tau and its neutrino excluded

  AllowedDecayToCount.push_back(23); //W 10-25s
  AllowedDecayToCount.push_back(24); //Z 

  AllowedDecayToCount.push_back(130); AllowedDecayToCount.push_back(310); AllowedDecayToCount.push_back(311); 
  //K0L, K0S,K0. 10-8 and 10-11 sec can be seen in T2. Have pi0 and charged decay mode: To review
   //K+ lifetime 10-8
  
  //AllowedDecayToCount.push_back(111); Pi0 Excluded
  
  
  AllowedDecayToCount.push_back(211); //Pi+- lifetime is 10-8 but gamma~10 for E= 1GeV Can be seen in T2.
  
  AllowedDecayToCount.push_back(113); //pi+ + pi- decay
  
  //WARNING: to rewievw: pi+ + pi- and 2pi0 decay
  AllowedDecayToCount.push_back(213);  //rho+
  
  AllowedDecayToCount.push_back(215); //a2 meson ?


 
  //Other quite rare I=1 mesons
  AllowedDecayToCount.push_back(117); AllowedDecayToCount.push_back(217); AllowedDecayToCount.push_back(119); AllowedDecayToCount.push_back(219); 

  //WARNING: Neutral decay mode is 72%. At the moment I drop it. Very short lifetime
  //AllowedDecayToCount.push_back(221); AllowedDecayToCount.push_back(331);//eta, eta prime
   

  AllowedDecayToCount.push_back(223); AllowedDecayToCount.push_back(333); //w phi
  AllowedDecayToCount.push_back(225); AllowedDecayToCount.push_back(335); //f2 f2'
  AllowedDecayToCount.push_back(227); //w3 phi3
  AllowedDecayToCount.push_back(337); AllowedDecayToCount.push_back(229);  
  AllowedDecayToCount.push_back(411); AllowedDecayToCount.push_back(421); AllowedDecayToCount.push_back(413); 
  AllowedDecayToCount.push_back(423); AllowedDecayToCount.push_back(415); AllowedDecayToCount.push_back(425); AllowedDecayToCount.push_back(431); AllowedDecayToCount.push_back(433); 
  AllowedDecayToCount.push_back(435); AllowedDecayToCount.push_back(511); AllowedDecayToCount.push_back(521); AllowedDecayToCount.push_back(513); AllowedDecayToCount.push_back(523); 
  AllowedDecayToCount.push_back(515); AllowedDecayToCount.push_back(525); AllowedDecayToCount.push_back(531); AllowedDecayToCount.push_back(533); AllowedDecayToCount.push_back(535); 
  AllowedDecayToCount.push_back(541); AllowedDecayToCount.push_back(543); AllowedDecayToCount.push_back(545); AllowedDecayToCount.push_back(441); AllowedDecayToCount.push_back(443); 
  AllowedDecayToCount.push_back(445); AllowedDecayToCount.push_back(551); AllowedDecayToCount.push_back(553); AllowedDecayToCount.push_back(555); AllowedDecayToCount.push_back(557); 
  
  AllowedDecayToCount.push_back(2224); AllowedDecayToCount.push_back(2214); AllowedDecayToCount.push_back(2114); 
  AllowedDecayToCount.push_back(1114); AllowedDecayToCount.push_back(3122); AllowedDecayToCount.push_back(3222); AllowedDecayToCount.push_back(3212); 
  AllowedDecayToCount.push_back(3112); AllowedDecayToCount.push_back(3224); AllowedDecayToCount.push_back(3214); AllowedDecayToCount.push_back(3114); AllowedDecayToCount.push_back(3322); 
  AllowedDecayToCount.push_back(3312); AllowedDecayToCount.push_back(3324); AllowedDecayToCount.push_back(3314); AllowedDecayToCount.push_back(3334); AllowedDecayToCount.push_back(4122); 
  AllowedDecayToCount.push_back(4222); AllowedDecayToCount.push_back(4212); AllowedDecayToCount.push_back(4112); AllowedDecayToCount.push_back(4224); AllowedDecayToCount.push_back(4214); 
  AllowedDecayToCount.push_back(4114); AllowedDecayToCount.push_back(4232); AllowedDecayToCount.push_back(4132); AllowedDecayToCount.push_back(4322); AllowedDecayToCount.push_back(4312); 
  AllowedDecayToCount.push_back(4324); AllowedDecayToCount.push_back(4314); AllowedDecayToCount.push_back(4332); AllowedDecayToCount.push_back(4334); AllowedDecayToCount.push_back(4412); 
  AllowedDecayToCount.push_back(4422); AllowedDecayToCount.push_back(4414); AllowedDecayToCount.push_back(4424); AllowedDecayToCount.push_back(4432); AllowedDecayToCount.push_back(4434); 
  AllowedDecayToCount.push_back(4444); AllowedDecayToCount.push_back(5122); AllowedDecayToCount.push_back(5112); AllowedDecayToCount.push_back(5212); AllowedDecayToCount.push_back(5222); 
  AllowedDecayToCount.push_back(5114); AllowedDecayToCount.push_back(5214); AllowedDecayToCount.push_back(5224); AllowedDecayToCount.push_back(5132); AllowedDecayToCount.push_back(5232); 
  AllowedDecayToCount.push_back(5312); AllowedDecayToCount.push_back(5322); AllowedDecayToCount.push_back(5314); AllowedDecayToCount.push_back(5324); AllowedDecayToCount.push_back(5332); 
  AllowedDecayToCount.push_back(5334); AllowedDecayToCount.push_back(5142); AllowedDecayToCount.push_back(5242); AllowedDecayToCount.push_back(5412); AllowedDecayToCount.push_back(5422); 
  AllowedDecayToCount.push_back(5414); AllowedDecayToCount.push_back(5424); AllowedDecayToCount.push_back(5342); AllowedDecayToCount.push_back(5432); AllowedDecayToCount.push_back(5434); 
  AllowedDecayToCount.push_back(5442); AllowedDecayToCount.push_back(5444); AllowedDecayToCount.push_back(5512); AllowedDecayToCount.push_back(5522); AllowedDecayToCount.push_back(5514); 
  AllowedDecayToCount.push_back(5524); AllowedDecayToCount.push_back(5532); AllowedDecayToCount.push_back(5534); AllowedDecayToCount.push_back(5542); AllowedDecayToCount.push_back(5544); 
  AllowedDecayToCount.push_back(5554);
































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

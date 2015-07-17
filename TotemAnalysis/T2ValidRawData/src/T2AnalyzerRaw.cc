//4_4_5_version

#include "TotemAnalysis/T2ValidRawData/interface/T2AnalyzerRaw.h"
#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include "CLHEP/Random/RandGaussQ.h"
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
#include "DataFormats/T2Cluster/interface/cluster_entry.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "FWCore/Framework/interface/MakerMacros.h"



T2AnalyzerRaw::T2AnalyzerRaw(const edm::ParameterSet& iConfig){

  t2PadDigiCollectionLabel = iConfig.getParameter<edm::InputTag>("T2PadDigiCollectionLabel");
  t2StripDigiCollectionLabel = iConfig.getParameter<edm::InputTag>("T2StripDigiCollectionLabel");
  t2DigiVfatCollectionLabel = iConfig.getParameter<edm::InputTag>("T2DigiVfatCollectionLabel");
  rawEventLabel = iConfig.getParameter<edm::InputTag>("RawEventLabel");

  simufile=iConfig.getParameter<bool>("simufile");
  produceVfatEffiFile=iConfig.getParameter<bool>("produceVfatEffiFile");
  

  Chicut = iConfig.getParameter<double>("Chicut");
  PhiChiProbCut = iConfig.getParameter<double>("PhiChiProbCut");
  RChiProbCut = iConfig.getParameter<double>("RChiProbCut");  
  chiRredCut = iConfig.getParameter<double>("chiRredCut");
  chiPhiredCut = iConfig.getParameter<double>("chiPhiredCut");
  AllowedDRTrackDistance=iConfig.getParameter<double>("AllowedDRTrackDistance");
  
  energy=iConfig.getParameter<double>("FitdzEnergy");
  DZScale=iConfig.getParameter<double>("DZScale");
  tracketamin=iConfig.getParameter<double>("TrkEtamin");
  tracketamax=iConfig.getParameter<double>("TrkEtaMAX");


  // numhittrkeff=iConfig.getParameter<int>("Numhittrkeff");
  // Testedcamera=iConfig.getParameter<std::vector<unsigned int> >("Testedcamera"); 

  MaxPad=iConfig.getParameter<unsigned int>("MaxPad");
  MaxStrip=iConfig.getParameter<unsigned int>("MaxStrip");
  MaxDphi=iConfig.getParameter<double>("MaxDphi");
  NumHitGood=iConfig.getParameter<unsigned int>("NumHitGood");
  maxdrhit=iConfig.getParameter<double>("maxdrhit");
  maxdphihit=iConfig.getParameter<double>("maxdphihit");
  Effgoodhitnumber=iConfig.getParameter<unsigned int>("Effgoodhitnumber");
  ExcludeNoisyplane=iConfig.getParameter<bool>("ExcludeNoisyplane");
  CommonNoiseClSize=iConfig.getParameter<unsigned int>("CommonNoiseClSize");

  outputFileName = iConfig.getUntrackedParameter<std::string>("OutputFile");
  CluLabel = iConfig.getParameter<std::string>("CluLabel");
  HitLabel = iConfig.getParameter<std::string>("HitLabel");
  RoadLabel = iConfig.getParameter<std::string>("RoadLabel");
  TrackLabel= iConfig.getParameter<std::string>("TrackLabel");

  xmlfilenameUsed_NotDead= iConfig.getParameter<std::string>("xmlfilenameUsed_NotDead");
  xmlfilenameFull= iConfig.getParameter<std::string>("xmlfilenameFull");
  DeadSectFileName= iConfig.getParameter<std::string>("DeadSectFileName");

  MaxEvents= iConfig.getParameter<unsigned int>("MaxEvents");

  Effmaxdphihit=iConfig.getParameter<double>("Effmaxdphihit");
  Effmaxdrhit=iConfig.getParameter<double>("Effmaxdrhit");
  EffMaxPad=iConfig.getParameter<unsigned int>("EffMaxPad");
  EffMaxStrip=iConfig.getParameter<unsigned int>("EffMaxStrip");

  NoiseDphiMAX=iConfig.getParameter<double>("NoiseDphiMAX");
  NoiseDrMAX=iConfig.getParameter<double>("NoiseDrMAX");


  HitNumb4Align=iConfig.getParameter<unsigned int>("HitNumb4Align");
  MeasuredXYResol=iConfig.getParameter<double>("MeasuredXYResol");
  SHIFTprescale=iConfig.getParameter<double>("SHIFTprescale");
  MaxStepalignstep=iConfig.getParameter<unsigned int>("MaxStepalignstep");
  Idreferencedet=iConfig.getParameter<unsigned int>("Idreferencedet");
  AlignmentHitRMax=iConfig.getParameter<double>("AlignmentHitRMax");
  UseJointProb=iConfig.getParameter<unsigned int>("UseJointProb");
  FitgravcenterZ=iConfig.getParameter<double>("FitgravcenterZ");
  

  DetForNoiseStudies=iConfig.getParameter<unsigned int>("DetForNoiseStudies");
  PhiMinForNoiseStudies=iConfig.getParameter<double>("PhiMinForNoiseStudies");
  PhiMaxForNoiseStudies=iConfig.getParameter<double>("PhiMaxForNoiseStudies");
  useRZforResol=iConfig.getParameter<unsigned int>("useRZforResol");
  SelectedHalf=iConfig.getParameter<unsigned int>("SelectedHalf");
  verbosity=iConfig.getParameter<bool>("verbosity");

  MinTrkInQuarter=iConfig.getParameter<unsigned int>("MinTrkInQuarter");
  MaxTrkInQuarter=iConfig.getParameter<unsigned int>("MaxTrkInQuarter");
  MaxPadAllowedInQuarter=iConfig.getParameter<unsigned int>("MaxPadAllowedInQuarter");//=40

  skipSelectedEvents=iConfig.getParameter<bool>("skipSelectedEvents");
  skipEventFileName= iConfig.getParameter<std::string>("skipEventFileName");
  LookToRawEvent= iConfig.getParameter<bool>("LookToRawEvent");
  UseUncorrupetdEventMap= iConfig.getParameter<bool>("UseUncorrupetdEventMap");
  requiregoodChi= iConfig.getParameter<bool>("requiregoodChi");
  DispVtx=iConfig.getParameter<bool>("DispVtx");


  OnlycorruptionAnalysis=iConfig.getParameter<bool>("OnlycorruptionAnalysis");
  OnlyClusterAnalysis=iConfig.getParameter<bool>("OnlyClusterAnalysis");
  VFATMonitoring=iConfig.getParameter<bool>("VFATMonitoring");
  //MakeVfatCumulativePlots=iConfig.getParameter<bool>("MakeVfatCumulativePlots");
 bunchesToAnalyse=iConfig.getParameter<std::vector<unsigned int>  >("bunchesToAnalyse");
}


T2AnalyzerRaw::~T2AnalyzerRaw()
{
}


void T2AnalyzerRaw::ConvertCluRPhiInEffiGeoCell(double Rcl, double Phicl, uint32_t thedet, int &rCell, int &phiCell)
{
  
  //  144 and 41 mm correspond to 5.3 and 6.5 134 5.35
  int rsect=0;int phisect=0;
  int plane0_39=0;
 
  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo=conv.GetT2Info(thedet);
  plane0_39=planeinfo.symb;
  double zpiano=14000.;
  zpiano=planeinfo.Zdet;

  //((expR1<46)||(expR1>132))//Cut used in trk selection

  double etacluster=-log((Rcl/fabs(zpiano))/2.0);
  
  // 
 rsect=(int)(((6.4-etacluster)/(6.4-5.35))*numRsectEffi);
if(etacluster<5.35)
    rsect=numRsectEffi;//0;
  if((etacluster>6.4))
    rsect=-1;//numRsectEffi-1;


 /*
if(Rcl<41.)
    rsect=0;
  if((Rcl>144.)||(rsect==numRsectEffi))
    rsect=numRsectEffi-1;

 */
/*
rsect=(int)(((6.5-etacluster)/(6.5-5.3))*numRsectEffi);
  if(Rcl<41.)
    rsect=0;
  if((Rcl>144.)||(rsect==numRsectEffi))
    rsect=numRsectEffi-1;
*/
  




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
      // phisect=(int)((Phicl/192.)*numPhisectEffi);
      if(Phicl<4.)
	phisect=0;
      if(Phicl>188)
	phisect=numPhisectEffi-1;
      
      if((Phicl>4.)&&(Phicl<188))
	phisect=1+(int)(((Phicl-4)/(192.0-4.0-4.0))*(numPhisectEffi-2));

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
      //phisect=(int)((Phicl/192.)*numPhisectEffi);
      if(Phicl<4.)
	phisect=0;
      if(Phicl>188)
	phisect=numPhisectEffi-1;
      
      if((Phicl>4.)&&(Phicl<188))
	phisect=1+(int)(((Phicl-4)/(192.0-4.0-4.0))*(numPhisectEffi-2));
    }

  }
  
  if(phisect==numPhisectEffi)
    phisect=numPhisectEffi-1;


 
  
  rCell=rsect;
  phiCell=phisect;
}


/*
{
  
  //  144 and 41 mm correspond to 5.3 and 6.5
  int rsect=0;int phisect=0;
  int plane0_39=0;
 
  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo=conv.GetT2Info(thedet);
  plane0_39=planeinfo.symb;
  double zpiano=14000.;
  zpiano=planeinfo.Zdet;



  double etacluster=-log((Rcl/fabs(zpiano))/2.0);
  
  rsect=(int)(((6.5-etacluster)/(6.5-5.3))*numRsectEffi);
 

  if(Rcl<41.)
    rsect=0;
  if((Rcl>144.)||(rsect==numRsectEffi))
    rsect=numRsectEffi-1;
  


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

    double Phiclbeg=Phicl;
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
*/







bool T2AnalyzerRaw::HitInDeadSector(unsigned int Rsector,unsigned int plane0_10){
  bool toret=false;

  //0-100-200-...-40000
  //unsigned int RsectorConverted=((SelectedHalf*10)+plane0_10)*100+Rsector;
  
  // if(plane0_10==1)
  //std::cout<<"looking for dead sector in Plane:"<<plane0_10<<" Sector:"<<Rsector<<std::endl;

  unsigned int planeconv=((SelectedHalf*10)+plane0_10);
  
  for(unsigned int j=0;j<VectDeadSect_Plane.size();j++){
    if(VectDeadSect_Plane.at(j)==planeconv) 
      if(VectDeadSect_Sector.at(j)==Rsector) 
	toret=true;
  }
  
  // if(toret==true)    
  //std::cout<<" Found dead sector for Plane:"<<plane0_10<<" Sector:"<<Rsector<<std::endl;

  return toret;
}


bool T2AnalyzerRaw::IsVfatMapped(int symbvfat)
{
  bool toret=false;
  map<unsigned int, VFATRegisters>::const_iterator dit = Map_NOTDEAD->readoutIdToRegisters.find(symbvfat);
	
  if(dit != Map_NOTDEAD->readoutIdToRegisters.end()) 
    {
      toret=true;
    }
  else 
    {
      toret=false;        
    }
  
    

  return toret;
}




//--------------------------------------------------------------------------------------------------------------------------------------
// BEGIN ALIGNMENT PART FUNCTIONS
//--------------------------------------------------------------------------------------------------------------------------------------
unsigned int RawtoSymb(uint32_t thedet)
{
  T2DetId converter;
  unsigned int pl=converter.plane(thedet);
  unsigned int pls=converter.planeSide(thedet);
  unsigned int ht=converter.halfTelescope(thedet);
  unsigned int arm=converter.arm(thedet);
  unsigned int symbolic=pl*2+pls+ht*10+20*arm;	  
 
  return symbolic;
}




std::vector<double> MyLinearfitX(std::vector<T2Hit> hitvec,unsigned int UseJointProb)
{  

double Sr=0.;
double Srz=0.;
double Szz_r=0.;
double Sz_r=0.; 
double S0_r=0.; 

unsigned int  sizeHitv=hitvec.size();
std::vector<double> r;
std::vector<double> z;

std::vector<double> er;
 std::vector<double> ez;

 double a_rz, b_rz;

  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
   
      r.push_back(hitvec[jj].GetHitX());
      

      z.push_back(hitvec[jj].GetHitZ());    
 
      double phirad=hitvec[jj].GetHitPhi()*3.14159/180.0;
      double sigmax=cos(phirad)*cos(phirad)*0.12*0.12+sin(phirad)*sin(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
      
      if(UseJointProb==1)
	sigmax=(0.12*hitvec[jj].GetHitR())*(0.12/**hitvec[jj].GetHitR()*/);
      
      er.push_back(sqrt(sigmax));
   
            
      Srz += r[jj]*z[jj]/er[jj]/er[jj];
      Szz_r += z[jj]*z[jj]/er[jj]/er[jj];
      Sz_r += z[jj]/er[jj]/er[jj];
      Sr += r[jj]/er[jj]/er[jj];
      S0_r += 1.0/er[jj]/er[jj];

  }


a_rz = (Srz*S0_r - Sz_r*Sr) / (Szz_r*S0_r - Sz_r*Sz_r);   // angular coefficient
b_rz = (Sr*Szz_r - Sz_r*Srz) / (Szz_r*S0_r - Sz_r*Sz_r);  // intercept   R=(a_rz)Z + b_rz

 double e_a_rz = sqrt( S0_r / (S0_r*Szz_r - Sz_r*Sz_r) );
 double e_b_rz = sqrt( Szz_r / (S0_r*Szz_r - Sz_r*Sz_r) );

 std::vector<double> vect;
 for (unsigned int i=0;i<5;i++)
   vect.push_back(0.0);

 vect[0]=a_rz;
 vect[1]=b_rz;
 vect[2]=e_a_rz;
 vect[3]=e_b_rz;

 double correl=(-1.0)*Sz_r*(1.0/(S0_r*Szz_r-(Sz_r*Sz_r)));
 vect[4]=correl;

 return vect;


}





std::vector<double> MyLinearfitY(std::vector<T2Hit> hitvec,unsigned int UseJointProb)
{  

double Sr=0.;
double Srz=0.;
double Szz_r=0.;
double Sz_r=0.; 
double S0_r=0.; 

unsigned int  sizeHitv=hitvec.size();
std::vector<double> r;
std::vector<double> z;

std::vector<double> er;
 std::vector<double> ez;

 double a_rz, b_rz;

  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      //   std::cout<<hitvec[jj].GetHitR()<<" "<<hitvec[jj].GetHitPhi()<<" "<<hitvec[jj].GetHitZ()<<std::endl;
      r.push_back(hitvec[jj].GetHitY());
      

      z.push_back(hitvec[jj].GetHitZ());    
 
      double phirad=hitvec[jj].GetHitPhi()*3.14159/180.0;
      double sigmay=sin(phirad)*sin(phirad)*0.12*0.12+cos(phirad)*cos(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
      
      if(UseJointProb==1)
	sigmay=(0.12*hitvec[jj].GetHitR())*(0.12/**hitvec[jj].GetHitR()*/);
	
      er.push_back(sqrt(sigmay));
      
      
      Srz += r[jj]*z[jj]/er[jj]/er[jj];
      Szz_r += z[jj]*z[jj]/er[jj]/er[jj];
      Sz_r += z[jj]/er[jj]/er[jj];
      Sr += r[jj]/er[jj]/er[jj];
      S0_r += 1.0/er[jj]/er[jj];

  }


a_rz = (Srz*S0_r - Sz_r*Sr) / (Szz_r*S0_r - Sz_r*Sz_r);   // angular coefficient
b_rz = (Sr*Szz_r - Sz_r*Srz) / (Szz_r*S0_r - Sz_r*Sz_r);  // intercept   R=(a_rz)Z + b_rz
 double e_a_rz = sqrt( S0_r / (S0_r*Szz_r - Sz_r*Sz_r) );
 double e_b_rz = sqrt( Szz_r / (S0_r*Szz_r - Sz_r*Sz_r) );

  

 std::vector<double> vect;
 for (unsigned int i=0;i<5;i++)
   vect.push_back(0.0);

 vect[0]=a_rz;
 vect[1]=b_rz;
 vect[2]=e_a_rz;
 vect[3]=e_b_rz;

 double correl=(-1.0)*Sz_r*(1.0/(S0_r*Szz_r-(Sz_r*Sz_r)));
 vect[4]=correl;


 return vect;


}




std::vector<double> chi2X(std::vector<T2Hit> hitvec,unsigned int UseJointProb)
{

  std::vector<double> vpar=MyLinearfitX(hitvec,UseJointProb);
  unsigned int  sizeHitv=hitvec.size();
  double chi2x_=0.;
  for(unsigned int m =0; m<sizeHitv; m++)
    {
      double phirad=hitvec[m].GetHitPhi()*3.14159/180.0;
      double sigmax=cos(phirad)*cos(phirad)*0.12*0.12+sin(phirad)*sin(phirad)*0.015*0.015*hitvec[m].GetHitR()*hitvec[m].GetHitR();  
      if(UseJointProb==1)
      	sigmax=(0.12*hitvec[m].GetHitR())*(0.12/**hitvec[m].GetHitR()*/);
      chi2x_=chi2x_+((vpar[0]*hitvec[m].GetHitZ()+vpar[1]-hitvec[m].GetHitX())/sigmax)*((vpar[0]*hitvec[m].GetHitZ()+vpar[1]-hitvec[m].GetHitX())/sigmax);
    }
  //  return chi2x_;
  double chi2prob;
  chi2prob=TMath::Prob(chi2x_,(sizeHitv-2));
  
  std::vector<double> retV;
  retV.push_back(chi2prob);
  retV.push_back(chi2x_);
  
   return retV;

}



std::vector<double> chi2Y(std::vector<T2Hit> hitvec,unsigned int UseJointProb)
{

  std::vector<double> vpar=MyLinearfitY(hitvec,UseJointProb);
  unsigned int  sizeHitv=hitvec.size();
  double chi2y_=0.;

  for(unsigned int m =0; m<sizeHitv; m++)
    {
      double phirad=hitvec[m].GetHitPhi()*3.14159/180.0;
      double sigmay=sin(phirad)*sin(phirad)*0.12*0.12+cos(phirad)*cos(phirad)*0.015*0.015*hitvec[m].GetHitR()*hitvec[m].GetHitR();
      if(UseJointProb==1)
	sigmay=(0.12*hitvec[m].GetHitR())*(0.12/**hitvec[m].GetHitR()*/);
      chi2y_=chi2y_+((vpar[0]*hitvec[m].GetHitZ()+vpar[1]-hitvec[m].GetHitY())/sigmay)*((vpar[0]*hitvec[m].GetHitZ()+vpar[1]-hitvec[m].GetHitY())/sigmay);
    }
  
  double chi2prob;
  chi2prob=TMath::Prob(chi2y_,(sizeHitv-2));
  
  std::vector<double> retV;
  retV.push_back(chi2prob);
  retV.push_back(chi2y_);
  
   return retV;
  //return chi2y_;

}

//--------------------------------------------------------------------------------------------------------------------------------------
// END ALIGNMENT PART FUNCTIONS
//--------------------------------------------------------------------------------------------------------------------------------------



std::vector<double> MyLinearfit(std::vector<T2Hit> hitvec,unsigned int UseJointProb)
{  

double Sr=0.;
double Srz=0.;
double Szz_r=0.;
double Sz_r=0.; 
double S0_r=0.; 

double Sphi=0.;
double Sphiz=0.;
double Szz_phi=0.;
double Sz_phi=0.; 
double S0_phi=0.; 

unsigned int  sizeHitv=hitvec.size();
 std::vector<double> r;   
std::vector<double> z;
std::vector<double> phi;
std::vector<double> er;
 std::vector<double> ez;
std::vector<double> ephi;
 double a_rz, b_rz, a_phiz, b_phiz;

  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      //   std::cout<<hitvec[jj].GetHitX()<<" "<<hitvec[jj].GetHitPhi()<<" "<<hitvec[jj].GetHitZ()<<std::endl;
      //r->x   phi->y
      r.push_back(hitvec[jj].GetHitX());
      
      phi.push_back(hitvec[jj].GetHitY());
      z.push_back(hitvec[jj].GetHitZ());    
 
      double phirad=hitvec[jj].GetHitPhi()*3.14159/180.0;
      
      //er.push_back(0.1);
      double sigmax=cos(phirad)*cos(phirad)*0.12*0.12+sin(phirad)*sin(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
      //dalla Formula di JointProb viene sigma2x=r*sigma2R
      if(UseJointProb==1)
      	sigmax=(0.12*hitvec[jj].GetHitR())*(0.12/**hitvec[jj].GetHitR()*/);

      er.push_back(sqrt(sigmax));
      

      double sigmay=sin(phirad)*sin(phirad)*0.12*0.12+cos(phirad)*cos(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
      if(UseJointProb==1)
	sigmay=(0.12*hitvec[jj].GetHitR())*(0.12/**hitvec[jj].GetHitR()*/);

      ephi.push_back(sqrt(sigmay));
      
     
      //ez.push_back(hitvec[jj].GetHitDZ());
            
      Srz += r[jj]*z[jj]/er[jj]/er[jj];
      Szz_r += z[jj]*z[jj]/er[jj]/er[jj];
      Sz_r += z[jj]/er[jj]/er[jj];
      Sr += r[jj]/er[jj]/er[jj];
      S0_r += 1.0/er[jj]/er[jj];



      Sphiz += phi[jj]*z[jj]/ephi[jj]/ephi[jj];
      Szz_phi += z[jj]*z[jj]/ephi[jj]/ephi[jj];
      Sz_phi += z[jj]/ephi[jj]/ephi[jj];
      Sphi += phi[jj]/ephi[jj]/ephi[jj];
      S0_phi += 1.0/ephi[jj]/ephi[jj];




  }


a_rz = (Srz*S0_r - Sz_r*Sr) / (Szz_r*S0_r - Sz_r*Sz_r);   // angular coefficient
b_rz = (Sr*Szz_r - Sz_r*Srz) / (Szz_r*S0_r - Sz_r*Sz_r);  // intercept   X=(a_rz)Z + b_rz



a_phiz = (Sphiz*S0_phi - Sz_phi*Sphi) / (Szz_phi*S0_phi - Sz_phi*Sz_phi);   // angular coefficient
b_phiz = (Sphi*Szz_phi - Sz_phi*Sphiz) / (Szz_phi*S0_phi - Sz_phi*Sz_phi);  // intercept   


 double e_a_rz = sqrt( S0_r / (S0_r*Szz_r - Sz_r*Sz_r) );
 double e_b_rz = sqrt( Szz_r / (S0_r*Szz_r - Sz_r*Sz_r) );

 double e_a_phiz = sqrt( S0_phi / (S0_phi*Szz_phi - Sz_phi*Sz_phi) );
 double e_b_phiz = sqrt( Szz_phi / (S0_phi*Szz_phi - Sz_phi*Sz_phi) );

 std::vector<double> vect;
 for (unsigned int i=0;i<10;i++)
   vect.push_back(0.0);

 vect[0]=a_rz;
 vect[1]=b_rz;
 vect[2]=a_phiz;
 vect[3]=b_phiz;

 vect[4]=e_a_rz;
 vect[5]=e_b_rz;
 vect[6]=e_a_phiz;
 vect[7]=e_b_phiz;

 double correl=(-1.0)*Sz_r*(1.0/(S0_r*Szz_r-(Sz_r*Sz_r)));
 vect[8]=correl;
 correl=(-1.0)*Sz_phi*(1.0/(S0_phi*Szz_phi-(Sz_phi*Sz_phi)));
 vect[9]=correl;

 //std::cout<<vect[0]<<" "<<vect[1]<<" "<<vect[2]<<" "<<vect[3]<<" "<<std::endl;
 return vect;

}





std::vector<double> MyLinearfitCorr(std::vector<T2Hit> hitvec,TMatrixD &par_covariance_matrix,double &chi2_)
{
  
  //std::cout<<" Inside MyLinearfitCorr .. init"<<std::endl;

  double sigmaR=0.12;
  double sigmaPhi=0.015;
  unsigned int  sizeHitv=hitvec.size();

  TMatrixD ParCov_matr(4,4); //matrice di covarianza dei parametri trovati;
  unsigned int sizeArighe=sizeHitv*2;
  TMatrixD A(sizeArighe,4);
  TMatrixD At(4,sizeArighe);
  //A ??la matrice per cui Mis=A(Param)
  TMatrixD Vy(sizeArighe,sizeArighe); //matrice di covarianza delle misure (una per ogni xy, quindi ??diag a blocchi);
  TMatrixD Vym1(sizeArighe,sizeArighe); 

  TMatrixD Ls(4,4);//((A^T)(Vy^-1)A)
  // TMatrixD Cs(4,4);//(A^T)(Vy^-1)
  TMatrixD Cs(4,sizeArighe);//(A^T)(Vy^-1)

  TVectorD FittedParam(4);
 
  TVectorD Yvect(sizeArighe);

  std::vector<TVectorD> allXY; //vettore degli (xi,yi);
  
  for(unsigned int k =0; k<sizeHitv; k++)
    {
      TVectorD insert(2);
      insert[0]=hitvec[k].GetHitX();
      insert[1]=hitvec[k].GetHitY();
      allXY.push_back(insert);
    }

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

  //std::cout<<" MyLinearfitCorr Start computation  .. "<<std::endl;

for(unsigned int k =0; k<sizeHitv; k++)
 {
   TMatrixD OneVy(2,2); 
   OneVy.Zero();
   double phirad=hitvec[k].GetHitPhi()*3.14159/180.0;
   double r=hitvec[k].GetHitR();
   OneVy(0,0)=sigmaR*sigmaR*cos(phirad)*cos(phirad)+r*r*sin(phirad)*sin(phirad)*sigmaPhi*sigmaPhi;
   OneVy(0,1)=cos(phirad)*sin(phirad)*(sigmaR*sigmaR-r*r*sigmaPhi*sigmaPhi);  
   OneVy(1,0)=OneVy(0,1);
   OneVy(1,1)=sigmaR*sigmaR*sin(phirad)*sin(phirad)+ r*r*cos(phirad)*cos(phirad)*sigmaPhi*sigmaPhi;
   
   //Invert Vy matrix
   TMatrixD OneVym1(2,2); 
   OneVym1=OneVy;
   Double_t deti;	
   OneVym1.Invert(&deti);
   
   //OneVy.Print();
   if(fabs(deti)<0.001)     
     std::cout<<" Possible Vy Zero Determinant!!!"<<std::endl;
    
   Vym1(2*k,2*k)=OneVym1(0,0);
   Vym1(2*k,2*k+1)= OneVym1(0,1); 
   Vym1(2*k+1,2*k)=OneVym1(1,0);
   Vym1(2*k+1,2*k+1)=OneVym1(1,1);
   

   Vy(2*k,2*k)=OneVy(0,0);
   Vy(2*k,2*k+1)= OneVy(0,1); 
   Vy(2*k+1,2*k)=OneVy(1,0);
   Vy(2*k+1,2*k+1)=OneVy(1,1);

 

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


 }

//std::cout<<" MyLinearfitCorr Parameter calculation  .. "<<std::endl;

//Nota che qui Vy ??in realt??Vy^-1 

 Ls(0,0)=Ls00(0,0);  Ls(0,1)=Ls00(0,1);   Ls(0,2)=Ls01(0,0);   Ls(0,3)=Ls01(0,1);
 Ls(1,0)=Ls00(1,0);  Ls(1,1)=Ls00(1,1);   Ls(1,2)=Ls01(1,0);   Ls(1,3)=Ls01(1,1); 

 Ls(2,0)=Ls10(0,0);  Ls(2,1)=Ls10(0,1);   Ls(2,2)=Ls11(0,0);   Ls(2,3)=Ls11(0,1);
 Ls(3,0)=Ls10(1,0);  Ls(3,1)=Ls10(1,1);   Ls(3,2)=Ls11(1,0);   Ls(3,3)=Ls11(1,1);

 
 Double_t determ;	
 Ls.Invert(&determ);
 if(fabs(determ)<0.001)
     std::cout<<" Possible Zero LS  Determinant!!!"<<std::endl;

 At.Zero();
 At.Transpose(A);
 Cs.Zero();
 Cs=At*Vym1;


 
 FittedParam= Ls*Cs*Yvect;
 
 //std::cout<<" MyLinearfitCorr FittedParam: "<<std::endl;
 //FittedParam.Print();

 ParCov_matr=At*Vym1*A;


 ParCov_matr.Invert(&determ);
 //std::cout<<" MyLinearfitCorr Error Matrix: "<<std::endl;
 //ParCov_matr.Print();

 if(fabs(determ)<0.001)
   std::cout<<" Possible Zero ParCov_matr  Determinant!!!"<<std::endl;  

 //Calcolo chi2
 double chi2=0.;
 TVectorD Residui(sizeArighe);
 Residui=(Yvect-A*FittedParam);
 
 //TVectorD ResiduiT;
 //ResiduiT.Transpose(Residui);
 TVectorD VdotRes(sizeArighe);
  VdotRes=(Vym1*Residui);
 chi2=Residui*VdotRes;

 //----------------------
 // Save the output
 //----------------------
 chi2_=chi2;
 
 if(chi2<0)
   {
     std::cout<<" WARNING: MyLinearfitCorr Chi2: "<<chi2<<std::endl;
     std::cout<<" Residui"<<std::endl;
     Residui.Print();
     std::cout<<" Cov Matrix Vy"<<std::endl;
     Vy.Print();
     std::cout<<" Vym1*Residui"<<std::endl;
     VdotRes.Print();

   }
 par_covariance_matrix=ParCov_matr;
 
 //std::cout<<" MyLinearfitCorr Error Matrix: "<<std::endl;
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

 return vect;

}




bool lookforactivepad(edm::Handle<T2HitCollection> t2hitcoll,T1T2Track TrkCit,unsigned int m)
{
 bool flag=false;
 return flag;
}
bool lookforactivestrip(edm::Handle<T2HitCollection> t2hitcoll,T1T2Track TrkCit,unsigned int m)
{
  bool flag=false;
  return flag;
}




std::vector<double>  ResiduiRPhi(std::vector<double> vpar,double x, double y,double z) //vpar[4]
{
  //double retdrdphi[2];
   std::vector<double> retdrdphi;
 for (unsigned int i=0;i<2;i++)
   retdrdphi.push_back(0.0);

  retdrdphi[0]=(vpar[0]*z+vpar[1]-x);
  retdrdphi[1]=(vpar[2]*z+vpar[3]-y);  
  return retdrdphi;
}
/*
double T2AnalyzerRaw::hitrealdphi(double phi1,double phi2)
{
  double dphi=0.;
  int q1=phi1/90;
  int q2=phi2/90;
  
  if(q1==q2)
    dphi=
  
  return dphi;
}
*/

bool T2AnalyzerRaw::PadIsNeighbour(std::vector<T2Hit>refhitv,T2Hit clusterchecked)
{
  bool CloseCluster=false;
  int row1=0; int row2=0;
 int col1=0; int col2=0;
  for(unsigned int h=0;h<refhitv.size();h++)
    {
      if(clusterchecked.GetHitNumPad()>0)
	if(refhitv.at(h).GetHitNumPad()>0)
	  for(unsigned int i=0;i<clusterchecked.ClusterPad_entries.size();i++)
	    {
	      row1=clusterchecked.ClusterPad_entries[i].rad_coord;
	      col1=clusterchecked.ClusterPad_entries[i].ang_coord;
	
     
	      for(unsigned int j=0;j<refhitv.at(h).ClusterPad_entries.size();j++)
		{
	    
		  row2=refhitv.at(h).ClusterPad_entries[j].rad_coord;//bug found on 3 jan 010. i instead of j
		  col2=refhitv.at(h).ClusterPad_entries[j].ang_coord;

	  
		  int odd_evencheck=RawtoSymb(refhitv.at(h).GetHitDetRawId())-RawtoSymb(clusterchecked.GetHitDetRawId());//cluplane
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

std::vector<double> T2AnalyzerRaw::ResiduiForPad(std::vector<T2Hit>refhitv, T2Hit hit,double &expPadXErr,double &expPadYErr)
{

  T2Hit extrhit=GiveExtrapolatedHit((RawtoSymb(hit.GetHitDetRawId())%10),refhitv);
  //Extrapolation done with x-y fit

  double totDx=0; double totDy=0;
  for(unsigned int i=0;i<refhitv.size();i++)
    {
      totDx+=refhitv.at(i).GetHitDX();
      totDy+=refhitv.at(i).GetHitDY();
    }
  
  totDx=totDx/refhitv.size();
  totDy=totDy/refhitv.size();


  expPadYErr=totDy;
  expPadXErr=totDx;

  /*
  T1T2Track thetrk= cutobj.TrackFromHits(true, refhitv);
  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo=conv.GetT2Info(referenceHits.at(0).GetHitDetRawId());
  int simbol0=planeinfo.symb;
  int hitsimbol=(planeinfo.symb/10)*10 + plane0_10;
  planeinfo=conv.GetT2Info(hitsimbol);  
  double expX=thetrk.X0()+thetrk.GetTx()*planeinfo.Zdet;
  double expY=thetrk.Y0()+thetrk.GetTy()*planeinfo.Zdet;
  */


  std::vector<double> retdrdphi;
  for (unsigned int i=0;i<2;i++)
    retdrdphi.push_back(0.0);
  

  retdrdphi[0]=extrhit.GetHitX()- hit.GetHitX();
  retdrdphi[1]=extrhit.GetHitY()- hit.GetHitY();

  /*
  retdrdphi[0]=extrhit.GetHitR()- hit.GetHitR();
  int refsize=refhitv.size();

  double dphi1=fabs(refhitv.at(refsize-1).GetHitPhi()- hit.GetHitPhi());
  double dphi2=fabs(refhitv.at(0).GetHitPhi()- hit.GetHitPhi());
  double dphi=10.;

  if(dphi1<dphi2)
    dphi=dphi1;
  else
    dphi=dphi2;

  //std::cout<<" dphi="<<dphi<<std::endl;

  retdrdphi[1]=dphi;
  */

  return retdrdphi;
  
}


std::vector<double> T2AnalyzerRaw::ResiduiForStrip(std::vector<T2Hit>refhitv, T2Hit hit)
{
  
  
  for(unsigned int oo=0;oo<refhitv.size();oo++)
    {
      if(refhitv[oo].GetHitNumPad()>1500)
	std::cout<<" Problem In HIT cluster_entri" <<refhitv[oo].ClusterStrip_entries.size()<<" Num strip: "<<refhitv[oo].GetHitNumStrip()<<" Num Pad:"<<refhitv[oo].GetHitNumPad()<<std::endl;       				 
    }
  
  T2Hit extrhit=GiveExtrapolatedHit((RawtoSymb(hit.GetHitDetRawId())%10),refhitv);
  std::vector<double> retdrdphi;
  for (unsigned int i=0;i<2;i++)
    retdrdphi.push_back(0.0);
  
  retdrdphi[0]=extrhit.GetHitR()- hit.GetHitR();
 
  //std::cout<<" DR= "<<retdrdphi[0]<<std::endl;
  
  unsigned int refsize=refhitv.size();
  
  int symbplane=(RawtoSymb(hit.GetHitDetRawId())%10);
  //std::cout<<" Here symbplane="<<symbplane<<std::endl;
  bool ispari;
  if((symbplane%2)==0)
    ispari=true;
  else
    ispari=false;
  
  //cerco un pari o un dispari.
  T2Hit usethis;
  bool sameparityhitfound=false;
  for(unsigned int i=0;i<refsize;i++)
    {
    

      int actsymb=RawtoSymb(refhitv.at(i).GetHitDetRawId())%10;
      //  std::cout<<" act-symb "<<actsymb<<" Num Strip:"<<refhitv.at(i).GetHitNumStrip()<<" Num Pad:"<<refhitv.at(i).GetHitNumPad()<<std::endl;

      if(((actsymb%2)==0)&&(ispari))
	{
	  sameparityhitfound=true;
	  usethis=refhitv.at(i);
	  //  std::cout<<"  Num Strip:"<<usethis.GetHitNumStrip()<<" Num Pad:"<<usethis.GetHitNumPad()<<std::endl;
	}
      if(((actsymb%2)!=0)&&(ispari==false))
	{
	  sameparityhitfound=true;
	  usethis=refhitv.at(i);
	  //std::cout<<"  Num Strip:"<<usethis.GetHitNumStrip()<<" Num Pad:"<<usethis.GetHitNumPad()<<std::endl;
	}
  
    }
  //std::cout<<" Here2"<<std::endl;


  if(sameparityhitfound==true)
    {
      //std::cout<<" Here3"<<std::endl;
      //VERIFICA CHE LA COLONNA SIA UGUALE between hit and usethis.
      std::vector<cluster_entry> entriesstripcl1= hit.ClusterStrip_entries;
      std::vector<cluster_entry> entriesstripcl2= usethis.ClusterStrip_entries;

      if((entriesstripcl1.size()==0)||(entriesstripcl2.size()==0))
	{
	    std::cout<<" Problem In HIT cluster_entries: fix the bug in the code!!  Hit-refhit strip size:"<<hit.ClusterStrip_entries.size()<<" "<<usethis.ClusterStrip_entries.size()<<" |||  Num strip: "<<usethis.GetHitNumStrip()<<" Num Pad:"<<usethis.GetHitNumPad()<<std::endl;
 	}
      
      if(entriesstripcl1.at(0).ang_coord==entriesstripcl2.at(0).ang_coord)
	retdrdphi[1]=0.1;//It only means that it Is small
      else
	retdrdphi[1]=970.1;//90.

    }
  else
    {
      //std::cout<<" Here4"<<std::endl;
      for(unsigned int i=0;i<refsize;i++)
	{
	  int actsymb=RawtoSymb(refhitv.at(i).GetHitDetRawId())%10;
	  if((symbplane%2)!=(actsymb%2))
	    {usethis=refhitv.at(i);
	       //VERIFICA CHE LA COLONNA SIA *DIVERSA* between hit and usethis.
	      std::vector<cluster_entry> entriesstripcl1= hit.ClusterStrip_entries;
	      std::vector<cluster_entry> entriesstripcl2= usethis.ClusterStrip_entries;
	      if(entriesstripcl1.at(0).ang_coord!=entriesstripcl2.at(0).ang_coord)
		retdrdphi[1]=0.1;//It only means that it Is small
	      else
		retdrdphi[1]=970.1;//90.
	    }
	}
    }


  


//  (RawtoSymb(hit.GetHitDetRawId())%10)
  /*
  if(fabs(hit.GetHitPhi()-refhitv.at(0).GetHitPhi())>180.)
    {
      if(hit.GetHitPhi())
    }

  if(fabs(refhitv.at(refsize-1).GetHitPhi()-refhitv.at(0).GetHitPhi())<10.)
    retdrdphi[1]=(refhitv.at(refsize-1).GetHitPhi()+refhitv.at(0).GetHitPhi())/2.0 - hit.GetHitPhi();
  else
    retdrdphi[1]=std::min(fabs(refhitv.at(refsize-1).GetHitPhi()- hit.GetHitPhi()),fabs(refhitv.at(0).GetHitPhi()- hit.GetHitPhi()));
  */
  return retdrdphi;
}


std::vector<double> T2AnalyzerRaw::ResiduiForC1HIT(std::vector<T2Hit>refhitv, T2Hit hit)
{
 T2Hit extrhit=GiveExtrapolatedHit((RawtoSymb(hit.GetHitDetRawId())%10),refhitv);
 std::vector<double> retdrdphi;
  for (unsigned int i=0;i<2;i++)
    retdrdphi.push_back(0.0);
  
  retdrdphi[0]=extrhit.GetHitR()- hit.GetHitR();

//  int refsize=refhitv.size();
//  double dphi1=fabs(refhitv.at(refsize-1).GetHitPhi()- hit.GetHitPhi());
//  double dphi2=fabs(refhitv.at(0).GetHitPhi()- hit.GetHitPhi());
//  double dphi=10.;
//  if(dphi1<dphi2)
//    dphi=dphi1;
//  else
//    dphi=dphi2;

  return retdrdphi;
}


bool  T2AnalyzerRaw::HitIsInTrackColl(T2Hit trackinghit,T1T2TrackCollection trackColl) 
{
  bool HitisAlreadyInTrk=false;
  T1T2TrackCollection::const_iterator TrkCit;
  for(TrkCit=trackColl.begin(); TrkCit!=trackColl.end(); TrkCit++){
    
    for(unsigned int u=0;u<(*TrkCit).GetHitEntries();u++)

    if (fabs(trackinghit.GetHitPhi()-(*TrkCit).GetHitT2(u).GetHitPhi())<0.001)
      if (fabs(trackinghit.GetHitR()-(*TrkCit).GetHitT2(u).GetHitR())<0.001)
	if (fabs(trackinghit.GetHitZ()-(*TrkCit).GetHitT2(u).GetHitZ())<0.001)
	  {
	    HitisAlreadyInTrk=true;
	  }
  
  }
  
  return HitisAlreadyInTrk;
}



T2Hit T2AnalyzerRaw::GiveExtrapolatedHit(int plane0_10, std::vector<T2Hit> referenceHits)
{
  T2SelectionCutUtils cutobj;
  T1T2Track thetrk= cutobj.TrackFromHits(true, referenceHits);
  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo=conv.GetT2Info(referenceHits.at(0).GetHitDetRawId());
  int hitsimbol=(planeinfo.symb/10)*10 + plane0_10;
  planeinfo=conv.GetT2Info(hitsimbol);
  
  double expX=thetrk.X0()+thetrk.GetTx()*planeinfo.Zdet;
  double expY=thetrk.Y0()+thetrk.GetTy()*planeinfo.Zdet;
  // std::cout<<" ^^^^^^^^^ ExtrXYZ:"<<expX<<" "<<expY<<" "<<planeinfo.Zdet<<std::endl;

  T2Hit thehit;
  thehit=T2Hit(expX,expY,planeinfo.Zdet,planeinfo.cmsswid);
  
  return thehit;
}


bool T2AnalyzerRaw::HitInQuarterAnalyzed(T2Hit &aHit){

  bool toret=false;

  if((SelectedHalf==0)||(SelectedHalf==2))
    if((aHit.GetHitPhi()<96)||(aHit.GetHitPhi()>264))
      toret=true;

  if((SelectedHalf==1)||(SelectedHalf==3))
    if((aHit.GetHitPhi()>96)&&(aHit.GetHitPhi()<276))
      toret=true;

  if(aHit.GetHitPhi()==246.964)
    if(toret)
      std::cout<<"HitInQuarterAnalyzed FAIL for 246.964"<<std::endl;
    else
      std::cout<<"HitInQuarterAnalyzed OK for 246.964"<<std::endl;

  return toret;

}



std::vector<bool> T2AnalyzerRaw::ExludeThisNoisyplaneForEffMeas(T2StripClusterCollection::const_iterator itstrip, T2PadClusterCollection::const_iterator itpad)
{
  bool FoundBIGNoisySTRIPCluster=false;
  bool FoundBIGNoisyPADCluster=false;
  
  vector<T2Cluster> stripClv = itstrip->second;
  for(unsigned int k=0;k<stripClv.size();k++){
    if(stripClv[k].GetNoOfEntries()>CommonNoiseClSize)
      FoundBIGNoisySTRIPCluster=true;
  }
    

  vector<T2Cluster> padClv = itpad->second;
  for(unsigned int k=0;k<padClv.size();k++){
    if(padClv[k].GetNoOfEntries()>CommonNoiseClSize)
      FoundBIGNoisyPADCluster=true;	  
  }
  
  std::vector<bool> toret;
  toret.push_back(FoundBIGNoisySTRIPCluster);
  toret.push_back(FoundBIGNoisyPADCluster);

  return toret;
}



bool T2AnalyzerRaw::HitIndetector(std::vector<double> vpar, unsigned int symbdetid)
{
  bool flag=false;
  double z=0.;
  double expY;
  double expX;
  double expRaggio;

  T2GeometryUtil conv;

  
  T2GeometryUtil::T2DetInfo detInfo2;
  detInfo2=conv.GetT2Info((symbdetid+SelectedHalf*10));
  z=detInfo2.Zdet;


  
  //if((mapsymbtoZ.find(symbdetid))!= mapsymbtoZ.end())
    {
      //  std::map<unsigned int,double>::iterator Itr;
      //Itr=mapsymbtoZ.find(symbdetid);      
      //z=Itr->second;
      expX = (vpar[0]*z+vpar[1]); 
      expY=(vpar[2]*z+vpar[3]);
      expRaggio=sqrt(expX*expX+expY*expY);
      
      //if((expX<140.)&&(expX>(-140.))) //144, 42
	{
	  //  if((fabs(expY)<140.)&&(fabs(expY)>(-140.)))
	    {
	      if((expRaggio>44.)&&(expRaggio<144.))
		{
		  flag=true;
		}
	      else
		if(verbosity) {
		  std::cout<<" Raggio Aspettato Fuori det: "<<expRaggio<<" z="<<z<<std::endl;
		  std::cout<<"X_ab:"<<vpar[0]<<" "<<vpar[1]<<"   Y_ab:"<<vpar[2]<<" "<<vpar[3]<<std::endl;
		}
	    }
	    //else
	    //if(verbosity) std::cout<<" Y Aspettato Fuori det: "<<expY<<std::endl;
	}
	// else
	//if(verbosity) std::cout<<" X Aspettato Fuori det:"<<expX<<std::endl;
    }

    /*
  else
    {
      //L'associazione non c'e' ancora.Metto flag=false per ignorare l'evento
      std::cout<<" HitIndetector: DETECTOR not yet ASSIGNED "<<std::endl;
      flag=false;
    }
    */

  //if(flag==false)
  //std::cout<<" Expected hit out of detector; X="<<expr<<"Y="<<expphi<<" R="<<expRaggio  <<" Z="<<z<<"  Symb="<<symbdetid<<std::endl;
 
  return flag;
}


bool T2AnalyzerRaw::PlaneInRightHalf(uint32_t cmsswid)
{
  bool ret=false;
  if((RawtoSymb(cmsswid)/10)==SelectedHalf)
    {
      ret=true;
      //  if(verbosity)
      //std::cout<<" PlaneInRightHalf "<<RawtoSymb(cmsswid)<<" Accepted"<<std::endl;
    }
  /*
  else
    {
      if(verbosity)
	std::cout<<" SYMB PLANE="<<RawtoSymb(cmsswid)<<" WHILE SELECTED-HALF="<<SelectedHalf<<std::endl;
    }
  */
  return ret; 
}


bool T2AnalyzerRaw::MettiTracciaBuona(std::vector<TrackInfo>* matricetracce,T1T2TrackCollection::const_iterator  trk, unsigned int numevento)
{
  bool flag=false;
  
  unsigned int counthitgood=0;
  std::vector<T2Hit> hitv;
  std::vector<double> vpar;

  
  //14057.1;
  double signz=1.0;
  if((*trk).GetHitT2(0).GetHitZ()<0)
     signz=-1.0;

  double centrodigravitaZ=FitgravcenterZ*signz;

   for(unsigned int a=0;a<(*trk).GetHitEntries(); a++)
      {
	//	if((*trk).GetHitT2(a).GetHitPlane()>4)
	  //  std::cout<<" Warning RECO!! Plane: "<<(*trk).GetHitT2(a).GetHitPlane()<<std::endl;
	if((((*trk).GetHitT2(a).GetHitNumStrip()==0)||((*trk).GetHitT2(a).GetHitNumPad()==0))&&((*trk).GetHitT2(a).GetHitClass()==1))
	  std::cout<<" MettitracciaBuona: Here THERE IS A BUG IN THE HIT!!"<<std::endl;
		  
	if(((*trk).GetHitT2(a).GetHitNumStrip()>513)||((*trk).GetHitT2(a).GetHitNumPad()>1600))
	  std::cout<<" MettitracciaBuona: Here THERE IS A BUG IN THE HIT!!"<<std::endl;

	if(PlaneInRightHalf((*trk).GetHitT2(a).GetHitDetRawId()))
	if((*trk).GetHitT2(a).GetHitClass()==1)
	  if((*trk).GetHitT2(a).GetHitNumStrip()<=MaxStrip)
	    if((*trk).GetHitT2(a).GetHitNumPad()<=MaxPad)
	      if((*trk).GetHitT2(a).GetHitDPhi()<=MaxDphi)
		{
		  hitv.push_back((*trk).GetHitT2(a));
		  hitv.at(counthitgood).SetHitZ((*trk).GetHitT2(a).GetHitZ()-centrodigravitaZ);
		  counthitgood++;
		}

	//	 if((RawtoSymb((*trk).GetHitT2(a).GetHitDetRawId())/10)==0)
	// if(*trk).GetHitT2(a).GetHitPhi()
      }

   
   double expX1=0.;
   double expY1=0.;
   double expR1=0.;
   double zz=0;
   for(unsigned int i=13700;i<=14200;i=i+100)
     {
       zz=(double)i;
       zz=zz*signz;
       expX1=(*trk).X0()+(*trk).GetTx()*zz;
       expY1=(*trk).Y0()+(*trk).GetTy()*zz; 
       expR1=sqrt(expX1*expX1+expY1*expY1);
       if((expR1<46)||(expR1>134))//Avoid to include the border effect in the effiestimation.
	 return false;
     
     }
   

   
   //  bool chicutpassed=false; 
   bool chicutpassed=true; //removed on 26 july 010
   //T2SelectionCutUtils cutobj;
   //T1T2Track thetrk= cutobj.TrackFromHits(true, referenceHits);
   /*
   if(((chi2X(hitv,UseJointProb).at(0))>0.001)&&((chi2Y(hitv,UseJointProb).at(0))>0.001))
     chicutpassed=true;
   else
     {
     if(verbosity)
       std::cout<<" Chi2-prob cut failed in Mettitracciabuona"<<std::endl;
     
     std::cout<<" Chi2-prob cut failed in Mettitracciabuona"<<std::endl;
     }
   */
   if((counthitgood>=NumHitGood)&&(chicutpassed))
      {
	//std::cout<<" counthitgood1 passes and chicutpassed"<<std::endl;
	TMatrixD covmat(4,4);
	covmat.Zero();
	double chi2corr;
	vpar=MyLinearfitCorr(hitv,covmat,chi2corr);
	
	if(useRZforResol==2)
	  vpar=MyLinearfit(hitv,UseJointProb);
	
	    //double *drdphi;//drdphi[2];
	std::vector<double>  drdphi;
	unsigned int symb;
	unsigned int countgoodhitaligned=0;
    
	//std::cout<<"  -  "<<vpar[0]<<" -  "<<vpar[1]<<" -  "<<vpar[2]<<" -  "<<vpar[3]<<std::endl;
	TrackInfo onetrack;
	std::vector<TrackInfo::IdandHit> tempidactive;
    
	for(unsigned int u=0;u<10;u++)
	  {
	    TrackInfo::IdandHit TrackInMatrixRow;
	    TrackInMatrixRow.iddet=0;
	    tempidactive.push_back(TrackInMatrixRow);	
	  }
    
	std::vector<T2Hit> hitv2;
	std::vector<double> vpar2; 
	std::vector<double> OLDvpar2;
	//double* vpar2; //4

	//	for(unsigned int a=0;a<(*trk).GetHitEntries(); a++)
	for(unsigned int a=0;a<hitv.size();a++)
	  {
	    // drdphi=ResiduiRPhi(vpar,(*trk).GetHitT2(a).GetHitX(),(*trk).GetHitT2(a).GetHitY(),(*trk).GetHitT2(a).GetHitZ());
	    //bug discovered on 22-05-09
	    drdphi=ResiduiRPhi(vpar,hitv.at(a).GetHitX(),hitv.at(a).GetHitY(),hitv.at(a).GetHitZ());
	    if(drdphi[0]<maxdrhit)
	      if(drdphi[1]<maxdphihit)
		{		
		  hitv2.push_back(hitv.at(a));
		  //hitv2.push_back((*trk).GetHitT2(a));	
		  countgoodhitaligned++;
		}
	  }
    
	
	//std::cout<<" B"<<std::endl;
	 if(verbosity) std::cout<<" countgoodhitaligned="<<countgoodhitaligned<<"  (NumHitGood= "<<NumHitGood<<")"<<std::endl;
	 //double C0=((*trk).GetHitT2(0).GetHitX()*(*trk).GetHitT2(0).GetHitX()+(*trk).GetHitT2(0).GetHitY()*(*trk).GetHitT2(0).GetHitY())/((*trk).GetTx()*(*trk).GetHitT2(0).GetHitX()+(*trk).GetTy()*(*trk).GetHitT2(0).GetHitY());

	(*trk).GetHitT2(0).GetHitZ();

	if(countgoodhitaligned>=NumHitGood)
	  {	       
	    onetrack.eventnumber=numevento;		
	    onetrack.goodhitnumber=countgoodhitaligned;
	

	    covmat.Zero();
	    double chi2corr;
	    vpar2=MyLinearfitCorr(hitv2,covmat,chi2corr);

	    if(useRZforResol==2)
	      vpar2=MyLinearfit(hitv,UseJointProb);
	   
	    //std::cout<<"  -  "<<vpar2[0]<<" -  "<<vpar2[1]<<" -  "<<vpar2[2]<<" -  "<<vpar2[3]<<std::endl;
	    for(unsigned int a=0;a<hitv2.size(); a++)
	      {	   
		symb=hitv2.at(a).GetHitPlane()*2+hitv2.at(a).GetHitPlaneSide();
		if((symb>9)||(hitv2.at(a).GetHitPlane()>4))
		  std::cout<<" WARNING"<<hitv2.at(a).GetHitPlane()<<hitv2.at(a).GetHitPlaneSide()<<std::endl;
		
		tempidactive.at(symb).iddet=1;		
		tempidactive.at(symb).thehit=hitv2.at(a);
		if(tempidactive.at(symb).iddet==1)
		  {
		    drdphi=ResiduiRPhi(vpar2,hitv2.at(a).GetHitX(),hitv2.at(a).GetHitY(),hitv2.at(a).GetHitZ());
		    tempidactive.at(symb).dr=drdphi[0];//dx
		    tempidactive.at(symb).dphi=drdphi[1];//dy
		  }
	
	      }
	
	
	    onetrack.ar=vpar2[0];
	    onetrack.br= vpar2[1];
	    onetrack.aphi=vpar2[2];
	    onetrack.bphi= vpar2[3];
	    /*
	    onetrack.OLDar=OLDvpar2[0];
	    onetrack.OLDbr= OLDvpar2[1];
	    onetrack.OLDaphi=OLDvpar2[2];
	    onetrack.OLDbphi= OLDvpar2[3];
	    */

	    onetrack.idactive=tempidactive;	   	   

	    (*matricetracce).push_back(onetrack);
	    flag=true;
	  }		
      }
   else
     {
       //std::cout<<" counthitgood1 NOT passes and chicutpassed"<<std::endl;
     }
    
    if(flag==true)
      {
	if(verbosity) std::cout<<" MettiTracciaBuona Flag-True"<<std::endl;
      }
    else
      if(verbosity) std::cout<<" MettiTracciaBuona Flag-False"<<std::endl; 
    
    return flag;
}


bool  T2AnalyzerRaw::CloseTrksInEvt(T1T2TrackCollection trackColl)
{
  bool toret=false;

  T1T2TrackCollection::const_iterator TrkCit;
  T1T2TrackCollection::const_iterator TrkCit2;

  for(TrkCit=trackColl.begin(); TrkCit!=trackColl.end(); TrkCit++){
    //std::cout<<" R | Phi:  "<<(*TrkCit).GetHitT2(0).GetHitR()<<"|"<< (*TrkCit).GetHitT2(0).GetHitPhi()<<"            ";
    if(PlaneInRightHalf((*TrkCit).GetHitT2(0).GetHitDetRawId()))
      for(TrkCit2=TrkCit; TrkCit2!=trackColl.end(); TrkCit2++){
     
	if(PlaneInRightHalf((*TrkCit2).GetHitT2(0).GetHitDetRawId()))
	  if(TrkCit2!=TrkCit)
	    if(fabs((*TrkCit).GetHitT2(0).GetHitR()-(*TrkCit2).GetHitT2(0).GetHitR())<30.)//Was 50  
	      if(fabs((*TrkCit).GetHitT2(0).GetHitPhi()-(*TrkCit2).GetHitT2(0).GetHitPhi())<15.)
		{
		  toret=true;
		  //std::cout<<" Close tracks! Evt discarded! "<<((*TrkCit).GetHitT2(0).GetHitR()-(*TrkCit2).GetHitT2(0).GetHitR())<<" "<<((*TrkCit).Phi()*180/3.14159265-(*TrkCit2).Phi()*180/3.14159265)<<std::endl;
		}
      }

  }

  // if(toret)
    //std::cout<<" Close tracks! Evt discarded! "<<std::endl;
  return toret;
}

//bool  T2AnalyzerRaw::HitIsInTrackColl(T2Hit trackinghit,T1T2TrackCollection trackColl) 

//--------------------------------------------------------------------------------------------------------------------------------------
// ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******
//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
// ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******
//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
// ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******  ANLYZER BEGIN ******
//--------------------------------------------------------------------------------------------------------------------------------------
 

void T2AnalyzerRaw::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace HepMC;

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  /* LOADING OF ALL THE RECORDS FROM THE EVENT */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  if(verbosity)
    {
      std::cout<<" -----------------------------------------------------------------------------"<<std::endl;
      std::cout<<" numevent: "<<numevent<<std::endl;
      std::cout<<" -----------------------------------------------------------------------------"<<std::endl;
    }

    /* :::::::::::::TakeDigi::::::::::::*/
  edm::Handle<T2PadDigiCollection> t2paddigicoll;
  iEvent.getByLabel(t2PadDigiCollectionLabel, t2paddigicoll);
  const T2PadDigiCollection* PadDigiptr;
  PadDigiptr= t2paddigicoll.product();
  edm::Handle<T2StripDigiCollection> t2stripdigicoll;
  iEvent.getByLabel(t2StripDigiCollectionLabel, t2stripdigicoll);
  const T2StripDigiCollection* StripDigiptr;
  StripDigiptr= t2stripdigicoll.product(); 
  DigiContainerIterator<T2DetId, T2PadDigi> itp;
  DigiContainerIterator<T2DetId, T2StripDigi> its;


     /*:::::::::::::TAKE VFAT DATA or SIMU:::::::::::::::*/
  
  edm::Handle<T2DigiVfatCollection> t2digivfatcoll;
  iEvent.getByLabel(t2DigiVfatCollectionLabel, t2digivfatcoll);
  const T2DigiVfatCollection* DigiVFatptr;
  DigiVFatptr= t2digivfatcoll.product();
  DigiContainerIterator<T2DetId, T2DigiVfat> itv;

  
  
  //std::cout<<"EventoInit:"<<std::endl;
  
  edm::Handle<T2VfatInformation> t2vfatinfo;
  const T2VfatInformation* t2vfatinfoptr;
  
  edm::Handle<Totem::RawEvent> inputRaw;

  /*
  if(((LookToRawEvent)||(UseUncorrupetdEventMap))&&(OnlycorruptionAnalysis==false))   //Added 1 august 2010 in order to work also with digi
    { 
      iEvent.getByType(inputRaw);
      iEvent.getByType(t2vfatinfo);
      
      t2vfatinfoptr= t2vfatinfo.product();     
      std::cout<<"Evento:"<<inputRaw->triggerData.bunch_num<<" "<<inputRaw->triggerData.input_status_bits<<std::endl;
      //unsigned int event_num, bunch_num, src_id;
      //unsigned int orbit_num;
      //unsigned char revision_num;
      //unsigned int run_num, trigger_num, inhibited_triggers_num, input_status_bits;
      BunchStatusBit->Fill(inputRaw->triggerData.bunch_num,inputRaw->triggerData.input_status_bits);
    }
  */
  
  
  /* :::::::::::::Take The Clusters::::::::::::*/
  Handle<T2StripClusterCollection> t2strclcoll;
  Handle<T2PadClusterCollection> t2padclcoll;
  /*::::::Take  T2  Hits::::::*/
  Handle<T2HitCollection> t2hitcoll;
    /*::::::Take  T2  Roads::::::*/
  Handle<T2RoadCollection> t2roadcoll;

 
  /*:::::: Take T2 tracks ::::::*/
  Handle<T1T2TrackCollection> trackCollection;
  bool GoodEventForTheQuarterExcludingDead=true;
  



  /*
  unsigned int thisbunch=inputRaw->triggerData.bunch_num;
  
  
  
  if((std::find(bunchesToAnalyse.begin(),bunchesToAnalyse.end(),thisbunch) == bunchesToAnalyse.end())&&(bunchesToAnalyse.size()!=0))
    return;
  */




  if((UseUncorrupetdEventMap)&&(OnlycorruptionAnalysis==false))
    {
  
      map<unsigned int, VFATRegisters>::const_iterator mapIterator;
      
      for( mapIterator = Map_FullVfats->readoutIdToRegisters.begin(); mapIterator != Map_FullVfats->readoutIdToRegisters.end(); mapIterator++ ) {
	int symbvfid=mapIterator->first;
	T2VfatInformation::const_iterator itvf;
	itvf= t2vfatinfoptr->find(symbvfid);

	if(itvf!=t2vfatinfo->end())
	  {
	    if((*itvf).second==5)
	      {
		if(unsigned((symbvfid/100)/10)==SelectedHalf)
		  {
		    if(verbosity)
		      std::cout<<" vfat "<<symbvfid%100<<"  in plane "<<symbvfid/100<<"  Not found in data"<<std::endl;
		  }
	      }
	    
	    if((*itvf).second==4)
	      {
		if(unsigned((symbvfid/100)/10)==SelectedHalf)
		  {
		    if(verbosity)
		      std::cout<<" vfat "<<symbvfid%100<<"  in plane "<<symbvfid/100<<"  Has been not included in the map"<<std::endl;
		  }
	      }
	    if((*itvf).second==3)
	      {
		if(unsigned((symbvfid/100)/10)==SelectedHalf)
		  {
		    //std::cout<<" vfat "<<symbvfid%100<<"  in plane "<<symbvfid/100<<"  CRC error"<<std::endl;
		  }
	      }
	    
	    if((*itvf).second==2)
	      {
		if(unsigned((symbvfid/100)/10)==SelectedHalf)
		  {
		    //std::cout<<" vfat "<<symbvfid%100<<"  in plane "<<symbvfid/100<<"  FootPrint error"<<std::endl;
		  }
	      }
	    if((*itvf).second==1)
	      {
		if(unsigned((symbvfid/100)/10)==SelectedHalf)
		  {
			std::cout<<" vfat "<<symbvfid%100<<"  in plane "<<symbvfid/100<<"  Data Mapping inconsistence Error"<<std::endl;
		  }
	      }

	    //Histogram
	    
	    if(IsVfatMapped(symbvfid)==false)//Is not in the alive list (IS DEad)
	      if((*itvf).second==0)          //Send Correct frames 
		{
		  T2vfatinfoErrorMap->Fill(symbvfid/100,symbvfid%100);//You can have dead vfat with correct frames
		  T2vfatinfoErrorMap_DEADWithCorrectFrame->Fill(symbvfid/100,symbvfid%100);
		}

	    if((*itvf).second<=5)
	      {
		if((*itvf).second>0)
		  {
		    T2vfatinfoErrorMap->Fill(symbvfid/100,symbvfid%100);

		    if(IsVfatMapped(symbvfid)==true){
		      T2vfatinfoErrorMapNOTDEAD->Fill(symbvfid/100,symbvfid%100);
		      if(unsigned((symbvfid/100)/10)==SelectedHalf)
			GoodEventForTheQuarterExcludingDead=false;
		    
		    }
		    //geometry histo
		    unsigned int quarter=((symbvfid/100)/10);
		    int plane=((symbvfid/100)%10);
		    int vfatiid=(symbvfid%100);
		    Double_t radialC=-1.;Double_t azimccord=-1.;
		    if(quarter==SelectedHalf){

		    if((plane%2)==0)
		      {
			if((vfatiid==0))
			  radialC=0;
			if((vfatiid==16))
			  radialC=16;
			if((vfatiid==1))
			  radialC=1;
			if((vfatiid==15))
			  radialC=15;
			
			if((vfatiid>1)&&(vfatiid<15))
			  azimccord=vfatiid-2;			 			
		      }
		    else
		      {
			if((vfatiid==0))
			  radialC=16;
			if((vfatiid==16))
			  radialC=0;
			if((vfatiid==1))
			  radialC=15;
			if((vfatiid==15))
			  radialC=1;

			if((vfatiid>1)&&(vfatiid<15))
			  azimccord=13-(vfatiid-2);
		      }
		    //    std::cout<<radialC<<" "<<azimccord<<std::endl;

		    if(azimccord!=(-1.))
		      {
			
			    T2vfatinfoErrorMapGeoQ1->Fill(azimccord,0.);
			    T2vfatinfoErrorMapGeoQ1->Fill(azimccord,1.);
			  
		      }

		    if(radialC!=(-1.))
		      {
		       
			if(radialC==0.)
			  {
			    T2vfatinfoErrorMapGeoQ1->Fill(2.,0.);T2vfatinfoErrorMapGeoQ1->Fill(3.,0.);
			    T2vfatinfoErrorMapGeoQ1->Fill(4.,0.);T2vfatinfoErrorMapGeoQ1->Fill(5.,0.);
			    T2vfatinfoErrorMapGeoQ1->Fill(6.,0.);T2vfatinfoErrorMapGeoQ1->Fill(7.,0.);
			  }
			if(radialC==1.)
			  {			
			    T2vfatinfoErrorMapGeoQ1->Fill(2.,1.);T2vfatinfoErrorMapGeoQ1->Fill(3.,1.);
			    T2vfatinfoErrorMapGeoQ1->Fill(4.,1.);T2vfatinfoErrorMapGeoQ1->Fill(5.,1.);
			    T2vfatinfoErrorMapGeoQ1->Fill(6.,1.);T2vfatinfoErrorMapGeoQ1->Fill(7.,1.);
			  }
			if(radialC==15.)
			  {
			    T2vfatinfoErrorMapGeoQ1->Fill(9.,1.);T2vfatinfoErrorMapGeoQ1->Fill(10.,1.);
			    T2vfatinfoErrorMapGeoQ1->Fill(11.,1.);T2vfatinfoErrorMapGeoQ1->Fill(12.,1.);
			    T2vfatinfoErrorMapGeoQ1->Fill(13.,1.);T2vfatinfoErrorMapGeoQ1->Fill(14.,1.);
			  }
			if(radialC==16.)
			  {			
			    T2vfatinfoErrorMapGeoQ1->Fill(9.,0.);T2vfatinfoErrorMapGeoQ1->Fill(10.,0.);
			    T2vfatinfoErrorMapGeoQ1->Fill(11.,0.);T2vfatinfoErrorMapGeoQ1->Fill(12.,0.);
			    T2vfatinfoErrorMapGeoQ1->Fill(13.,0.);T2vfatinfoErrorMapGeoQ1->Fill(14.,0.);
			  }

		      }
		    }
		    
		  }

	
		if((*itvf).second==2)
		  {
		    T2vfatinfoErrorMapFootPr->Fill(symbvfid/100,symbvfid%100);
		  }
		if((*itvf).second==3)
		  {
		    T2vfatinfoErrorMapCRC->Fill(symbvfid/100,symbvfid%100);
		  }
		

		if((*itvf).second==5)
		  {
		    T2vfatinfoErrorNotInData->Fill(symbvfid/100,symbvfid%100);
		  }

		if((*itvf).second==4)
		  {
		    T2vfatinfoErrorMapNotIncluded->Fill(symbvfid/100,symbvfid%100);
		  }
		if((*itvf).second==1)
		  {
		    T2vfatinfoErrorErrorMapping->Fill(symbvfid/100,symbvfid%100);
		  }
	      }
	    else
	      std::cout<<" Error in t2vfatinfo Reco analysis"<<std::endl;
	   
	  }
	  else
	    {
	      std::cout<<" VFat in FullXMLMAP not found boolean map.. Non Dovrebbe accadere perch??la boleana dovrebbe avere  17*40 vfat"<<"FullMapiSymbolicId="<<symbvfid<<std::endl;
	    }
	//vfatsSymbId.push_back(symbvfid);
      }
      //  std::cout<<" Found "<<vfatsSymbId.size()<<" vfats in the map utilized for efficiency caluclation"<<std::endl;

    }



   


  if(OnlycorruptionAnalysis==false)
    {
      

      iEvent.getByLabel(CluLabel,"T2StripClusters",t2strclcoll);
      
      iEvent.getByLabel(CluLabel,"T2PadClusters",t2padclcoll);  

      if(OnlyClusterAnalysis==false)
	{
	  iEvent.getByLabel(HitLabel,"T2Hits",t2hitcoll);
	  
	  iEvent.getByLabel(RoadLabel,"T2RoadColl",t2roadcoll);
	  
	  iEvent.getByLabel(TrackLabel,"T2TrackColl",trackCollection);
	}
      //std::cout<<"Tracce Prese!!"<<std::endl;
    }
  else
    return;




 


  
  bool alreadyprint=false;

  T2PadDigiCollection::DigiRangeIterator detUnitItP;
  //std::vector<int> howmanypads;

  /*
  for (detUnitItP=t2paddigicoll->begin(); detUnitItP != t2paddigicoll->end(); ++detUnitItP) {

    const T2PadDigiCollection::Range& range = (*detUnitItP).second;
    unsigned int padCounterP=0;

    T2DetId *detID =new T2DetId((*detUnitItP).first);
    //const T2DetId& detID = (*detUnitItP).first;
    uint32_t cmsswdId= detID->calculateRawId(detID->arm(),detID->halfTelescope(),detID->plane(),detID->planeSide());
    unsigned int symbol=RawtoSymb(cmsswdId);

    for (T2PadDigiCollection::const_iterator digiItP = range.first; digiItP != range.second; ++digiItP) {
      //USA:
      //listOfTrackingplanesinWantedQuarter
      //if(digiItP->getRow()<18)
	padCounterP++;
      //if(numevent==263)
	// if(symbol<10)
	//std::cout<<digiItP->getRow()<<" "<<digiItP->getCol()<<" plane: "<<symbol<<std::endl;      
    }
    
    if(std::find(listOfTrackingplanesinWantedQuarter.begin(), listOfTrackingplanesinWantedQuarter.end(), symbol)!=listOfTrackingplanesinWantedQuarter.end())
      {
	if((symbol/10)==SelectedHalf)
	  if(padCounterP>10)
	    Planelowmultiplicityflag=false;
	
	//listOfTrackingplanesinWantedQuarter

	numpadperplane[symbol]=padCounterP;
    
	if((symbol/10)==SelectedHalf){
	  totalPadInQuarter+=padCounterP;   
	  numplaneActive++;
	}
      }

  }

  bool toomanypad=false;
  bool toomanystrip=false;

  if((totalPadInQuarter>MaxPadAllowedInQuarter)||(Planelowmultiplicityflag==false))
      {
	toomanypad=true;
	if(verbosity)
	  std::cout<<" Numevent "<<numevent<<" with too many pad : "<<totalPadInQuarter<<" ON of the wanted Q.( "<<numplaneActive<<" different planes active)"<<std::endl;
      }
    else
      if(verbosity)
	std::cout<<" Low multiplicty event: "<<numevent<<" with "<<totalPadInQuarter<<" pad ON of the wanted Q.("<<numplaneActive<<" different planes active)"<<std::endl;

    
  
   
    
    
    T2StripDigiCollection::DigiRangeIterator detUnitItS;
    for (detUnitItS=t2stripdigicoll->begin(); detUnitItS != t2stripdigicoll->end(); ++detUnitItS) {

      const T2StripDigiCollection::Range& range = (*detUnitItS).second;
 
      //T2DetId *detID =new T2DetId((*detUnitItS).first);
      //uint32_t cmsswdId= detID->calculateRawId(detID->arm(),detID->halfTelescope(),detID->plane(),detID->planeSide());

      T2DetId detID;
      detID=(*detUnitItS).first;
      uint32_t cmsswdId= detID.calculateRawId(detID.arm(),detID.halfTelescope(),detID.plane(),detID.planeSide());

      
      unsigned int symbol=RawtoSymb(cmsswdId);
      unsigned int oneplanecount=0;
    
      if(std::find(listOfTrackingplanesinWantedQuarter.begin(), listOfTrackingplanesinWantedQuarter.end(), symbol)!=listOfTrackingplanesinWantedQuarter.end())
	for (T2StripDigiCollection::const_iterator digiItS = range.first; digiItS != range.second; ++digiItS) {	 
	  //std::cout<<" ev300 Strip at R "<<digiItS->getRow()<<" plane: "<<symbol<<std::endl;
	  oneplanecount++;
	} 
      
      if((symbol/10)==SelectedHalf)
	if(oneplanecount>1024)
	  toomanystrip=true;
	
    }
*/

    /*
    if(toomanystrip==true)
      if(verbosity)
	std::cout<<" Event with too many Strips"<<std::endl;
    
    if(toomanypad==true)
      if(verbosity)
	std::cout<<" Event with too many Pads"<<std::endl;
    */
   


 std::vector<unsigned int> idexplored;
    double rr=0.;double ff=0.;double xx=0.;double yy=0.;

 for(T2StripClusterCollection::const_iterator itstrip = t2strclcoll->begin(); itstrip != t2strclcoll->end(); itstrip++){
    

   vector<T2Cluster> stripClv = itstrip->second;
    T2DetId *detID =new T2DetId(itstrip->first);
    
    uint32_t cmsswdId= detID->calculateRawId(detID->arm(),detID->halfTelescope(),detID->plane(),detID->planeSide());
    unsigned int symbol=RawtoSymb(cmsswdId);

    unsigned int quarter=symbol/10;
    unsigned int plane=symbol%10;

    /*  
	NumPadCluVsPlaneAll3H0_CuttedMinus->Write();
	NumPadCluVsPlaneAll3H1_CuttedMinus->Write();
	NumPadCluVsPlaneAll3H2_CuttedMinus->Write();
	NumPadCluVsPlaneAll3H3_CuttedMinus->Write();
	NumStripCluVsPlaneAll3H0_CuttedMinus->Write();
	NumStripCluVsPlaneAll3H1_CuttedMinus->Write();
	NumStripCluVsPlaneAll3H2_CuttedMinus->Write();
	NumStripCluVsPlaneAll3H3_CuttedMinus->Write();
	
    */ 
    
    unsigned int countCutted0=0;
    unsigned int countCutted1=0;
    unsigned int countCutted2=0;
    unsigned int countCutted3=0;
    unsigned int countCutted0Plus=0;unsigned int countCutted0Minus=0;
    unsigned int countCutted1Plus=0;unsigned int countCutted1Minus=0;
    unsigned int countCutted2Plus=0;unsigned int countCutted2Minus=0;
    unsigned int countCutted3Plus=0;unsigned int countCutted3Minus=0;
    idexplored.push_back(symbol);
    /*
      rr=padClv[k].GetClusterR();
      ff=padClv[k].GetClusterPhi()*3.14159/180.;    
      xx=rr*cos(ff);
      yy=rr*sin(ff);
    */
    for(unsigned int k=0;k<stripClv.size();k++){
      
      rr=stripClv[k].GetClusterR();
      ff=stripClv[k].GetClusterPhi()*3.14159/180.;    
      xx=rr*cos(ff);
      yy=rr*sin(ff);

      if(quarter==0)
	{
	  StripCluSizeVsPlaneAll3H0->Fill(plane,stripClv[k].GetNoOfEntries());
	  CumulativeStripCluSizeAll3H0->Fill(stripClv[k].GetNoOfEntries());
	  if(stripClv[k].GetNoOfEntries()<=5){
	    countCutted0++;
	    if(yy>0)
	      countCutted0Plus++;
	    else
	      countCutted0Minus++;
	  }
	}
      
      if(quarter==1)
	{
	  StripCluSizeVsPlaneAll3H1->Fill(plane,stripClv[k].GetNoOfEntries());
	  if(stripClv[k].GetNoOfEntries()<=5){
	    countCutted1++;
	    if(yy>0)
	      countCutted1Plus++;
	    else
	      countCutted1Minus++;
	  }
	}

      if(quarter==2)
	{
	  StripCluSizeVsPlaneAll3H2->Fill(plane,stripClv[k].GetNoOfEntries());
	  if(stripClv[k].GetNoOfEntries()<=5){
	    countCutted2++;
	    if(yy>0)
	      countCutted2Plus++;
	    else
	      countCutted2Minus++;
	  }
	}

      if(quarter==3)
	{
	  StripCluSizeVsPlaneAll3H3->Fill(plane,stripClv[k].GetNoOfEntries());
	  if(stripClv[k].GetNoOfEntries()<=5){
	    countCutted3++;
	  
	    if(yy>0)
	      countCutted3Plus++;
	    else
	      countCutted3Minus++;
	  }
	}

      //if(symbol<10)
      //if(verbosity)
      //  std::cout<<" ev300Cluster Strip at "<<stripClv[k].GetClusterR()<<" plane: "<<symbol<<std::endl;
    }


    if(quarter==0){
      NumStripCluVsPlaneAll3H0->Fill(plane,stripClv.size());
      NumStripCluVsPlaneAll3H0_Cutted->Fill(plane,countCutted0);
      NumStripCluVsPlaneAll3H0_CuttedMinus->Fill(plane,countCutted0Minus);
      NumStripCluVsPlaneAll3H0_CuttedPlus->Fill(plane,countCutted0Plus);	
      CumulativeNumStripCluAll3H0->Fill(stripClv.size());
    }
    
    if(quarter==1){
      NumStripCluVsPlaneAll3H1->Fill(plane,stripClv.size());
      NumStripCluVsPlaneAll3H1_Cutted->Fill(plane,countCutted1);
      NumStripCluVsPlaneAll3H1_CuttedMinus->Fill(plane,countCutted1Minus);
      NumStripCluVsPlaneAll3H1_CuttedPlus->Fill(plane,countCutted1Plus);
    }

    if(quarter==2){
      NumStripCluVsPlaneAll3H2->Fill(plane,stripClv.size());
      NumStripCluVsPlaneAll3H2_Cutted->Fill(plane,countCutted2);
      NumStripCluVsPlaneAll3H2_CuttedMinus->Fill(plane,countCutted2Minus);
      NumStripCluVsPlaneAll3H2_CuttedPlus->Fill(plane,countCutted2Plus);
    }

    if(quarter==3){
      NumStripCluVsPlaneAll3H3->Fill(plane,stripClv.size());
      NumStripCluVsPlaneAll3H3_Cutted->Fill(plane,countCutted3);
      NumStripCluVsPlaneAll3H3_CuttedMinus->Fill(plane,countCutted3Minus);
      NumStripCluVsPlaneAll3H3_CuttedPlus->Fill(plane,countCutted3Plus);
    }


  }

 for(unsigned int i=0; i<40;i++)
   {
     if (std::find(idexplored.begin(),idexplored.end(),i)==idexplored.end())
       {
	 unsigned int qq=i/10;
	 unsigned int plane=i%10;
	 if(qq==0){
	   NumStripCluVsPlaneAll3H0->Fill(plane,0);
	   
	 }
	 if(qq==1){
	   NumStripCluVsPlaneAll3H1->Fill(plane,0);
	 }
	 if(qq==2){
	   NumStripCluVsPlaneAll3H2->Fill(plane,0);
	 }
	 if(qq==3){
	   NumStripCluVsPlaneAll3H3->Fill(plane,0);              
	 }
       }
   }
 
 idexplored.clear();

 for(T2PadClusterCollection::const_iterator itpad= t2padclcoll->begin(); itpad != t2padclcoll->end(); itpad++){
    vector<T2Cluster> padClv = itpad->second;
    T2DetId *detID =new T2DetId(itpad->first);
    
    uint32_t cmsswdId= detID->calculateRawId(detID->arm(),detID->halfTelescope(),detID->plane(),detID->planeSide());
    unsigned int symbol=RawtoSymb(cmsswdId);
    unsigned int quarter=symbol/10;
    unsigned int plane=symbol%10;

    unsigned int countCutted0=0;
    unsigned int countCutted1=0;
    unsigned int countCutted2=0;
    unsigned int countCutted3=0;
  unsigned int countCutted0Plus=0;unsigned int countCutted0Minus=0;
    unsigned int countCutted1Plus=0;unsigned int countCutted1Minus=0;
    unsigned int countCutted2Plus=0;unsigned int countCutted2Minus=0;
    unsigned int countCutted3Plus=0;unsigned int countCutted3Minus=0;
    idexplored.push_back(symbol);

      for(unsigned int k=0;k<padClv.size();k++){
	rr=padClv[k].GetClusterR();
	ff=padClv[k].GetClusterPhi()*3.14159/180.;    
	xx=rr*cos(ff);
	yy=rr*sin(ff);

	
	PadClusterR_AllvsPlane[symbol]->Fill(padClv[k].GetClusterR());
	PadClusterSize_AllvsPlane[symbol]->Fill(padClv[k].GetNoOfEntries());

	if(quarter==0)
	  {
	    PadCluSizeVsPlaneAll3H0->Fill(plane,padClv[k].GetNoOfEntries());
	    CumulativePadCluSizeAll3H0->Fill(padClv[k].GetNoOfEntries());
	    if(padClv[k].GetNoOfEntries()<=5){
	      countCutted0++;
	      if(yy>0)
		countCutted0Plus++;
	      else
		countCutted0Minus++;
	    }
	  }
      
	if(quarter==1)
	  {
	    PadCluSizeVsPlaneAll3H1->Fill(plane,padClv[k].GetNoOfEntries());
	    if(padClv[k].GetNoOfEntries()<=5){
	      countCutted1++;
	       if(yy>0)
		countCutted1Plus++;
	      else
		countCutted1Minus++;
	    }
	  }

	if(quarter==2)
	  {
	    PadCluSizeVsPlaneAll3H2->Fill(plane,padClv[k].GetNoOfEntries());
	    if(padClv[k].GetNoOfEntries()<=5){
	      countCutted2++;
	       if(yy>0)
		countCutted2Plus++;
	       else
		countCutted2Minus++;
	    }
	  }
	
	if(quarter==3)
	  {
	    PadCluSizeVsPlaneAll3H3->Fill(plane,padClv[k].GetNoOfEntries());
	    if(padClv[k].GetNoOfEntries()<=5){
	      countCutted3++;
	      if(yy>0)
		countCutted3Plus++;
	      else
		countCutted3Minus++;
	    }
	  }
		
    }


      if(quarter==0)
	{
	  NumPadCluVsPlaneAll3H0->Fill(plane,padClv.size());
	  CumulativeNumPadCluAll3H0->Fill(padClv.size());
	  NumPadCluVsPlaneAll3H0_Cutted->Fill(plane,countCutted0);
	  NumPadCluVsPlaneAll3H0_CuttedMinus->Fill(plane,countCutted0Minus);
	  NumPadCluVsPlaneAll3H0_CuttedPlus->Fill(plane,countCutted0Plus);
	  
	  if(countCutted0<5)
	    NumPadCluVsPlaneAll3H0_Cutted_LOWMultipl->Fill(plane,countCutted0);
	}
      
      if(quarter==1)
	{
	  NumPadCluVsPlaneAll3H1->Fill(plane,padClv.size());
	  NumPadCluVsPlaneAll3H1_Cutted->Fill(plane,countCutted1);
	  NumPadCluVsPlaneAll3H1_CuttedMinus->Fill(plane,countCutted1Minus);
	  NumPadCluVsPlaneAll3H1_CuttedPlus->Fill(plane,countCutted1Plus);

	  if(countCutted1<5)
	    NumPadCluVsPlaneAll3H1_Cutted_LOWMultipl->Fill(plane,countCutted1);
	}

      if(quarter==2)
	{
	  NumPadCluVsPlaneAll3H2->Fill(plane,padClv.size());
	  NumPadCluVsPlaneAll3H2_Cutted->Fill(plane,countCutted2);
	  NumPadCluVsPlaneAll3H2_CuttedMinus->Fill(plane,countCutted2Minus);
	  NumPadCluVsPlaneAll3H2_CuttedPlus->Fill(plane,countCutted2Plus);

	  if(countCutted2<5)
	    NumPadCluVsPlaneAll3H2_Cutted_LOWMultipl->Fill(plane,countCutted2);
	}

      if(quarter==3)
	{
	  NumPadCluVsPlaneAll3H3->Fill(plane,padClv.size());   
	  NumPadCluVsPlaneAll3H3_Cutted->Fill(plane,countCutted3);
	   NumPadCluVsPlaneAll3H3_CuttedMinus->Fill(plane,countCutted3Minus);
	  NumPadCluVsPlaneAll3H3_CuttedPlus->Fill(plane,countCutted3Plus);
	  if(countCutted3<5)
	    NumPadCluVsPlaneAll3H3_Cutted_LOWMultipl->Fill(plane,countCutted3);
	}

  }

 for(unsigned int i=0; i<40;i++)
   {
     if (std::find(idexplored.begin(),idexplored.end(),i)==idexplored.end())
       {
	 unsigned int qq=i/10;
	 unsigned int plane=i%10;
	 if(qq==0)
	   NumPadCluVsPlaneAll3H0->Fill(plane,0);
	 if(qq==1)
	   NumPadCluVsPlaneAll3H1->Fill(plane,0);
	 if(qq==2)
	   NumPadCluVsPlaneAll3H2->Fill(plane,0);
	 if(qq==3)
	   NumPadCluVsPlaneAll3H3->Fill(plane,0);              
       }
   }







 
 if(OnlyClusterAnalysis==true){
   numevent++;
   return;
 }



//--*****************************        Reco TRACK        ********************************--/
  std::vector<unsigned int> listOfTrackingplanesinWantedQuarter;

  T1T2TrackCollection::const_iterator TrkCit;
  T1T2TrackCollection  t2trackVectorALL;
   //I'm looking if the event is too schifoso
  unsigned int trkcountergd=0;
  unsigned int trkcounterallInQuarterAnalyzed=0;
  unsigned int trkcounterallEVERYQUARTER=0;

  T2SelectionCutUtils T2CutsUtil_ForChi,T2CutsUtil;

  unsigned int count_t2trackVectorALL_H0=0; 
  unsigned int count_t2trackVectorALL_H1=0; 
  unsigned int count_t2trackVectorALL_H2=0; 
  unsigned int count_t2trackVectorALL_H3=0;
  std::vector<int> dummy_0; dummy_0.push_back(0);	  
  std::vector<int> dummy_1; dummy_1.push_back(1);
  std::vector<int> dummy_2; dummy_2.push_back(2);
  std::vector<int> dummy_3; dummy_3.push_back(3);

  if(verbosity)
    std::cout<<"Numero Tracce  "<<  trackCollection->end()- trackCollection->begin() <<std::endl;



  unsigned int trkInplus=0; unsigned int trkInminus=0;

  for(TrkCit=trackCollection->begin(); TrkCit!=trackCollection->end(); TrkCit++){     

    //     std::cout<<"  ?????????????????????????????????????????????? $$$$$$$$$$$$$$$$ "<<(*TrkCit).GetHitT2(0).GetHitHalftele()<<std::endl;
    t2trackVectorALL.push_back((*TrkCit));
    //  if(((*TrkCit).GetHitT2(0).GetHitHalftele())==1)
    //  std::cout<<"  ?????????????????????????????????????????????? $$$$$$$$$$$$$$$$ "<<(*TrkCit).GetHitT2(0).GetHitHalftele()<<std::endl;
    unsigned int symb=RawtoSymb(((*TrkCit).GetHitT2(0).GetHitDetRawId()));

    if((*TrkCit).Eta()>0){
      trkInplus++;
      TrkphiRZPlus->Fill((*TrkCit)._phiRZ*180.0/3.14159);
      if(T2CutsUtil.TrkAlsoInQuarter((*TrkCit),dummy_0))
	count_t2trackVectorALL_H0++;
      if(T2CutsUtil.TrkAlsoInQuarter((*TrkCit),dummy_1))  
	count_t2trackVectorALL_H1++;	      
    }
    else{
      trkInminus++;
      TrkphiRZMinus->Fill((*TrkCit)._phiRZ*180.0/3.14159);
      if(T2CutsUtil.TrkAlsoInQuarter((*TrkCit),dummy_2))
	count_t2trackVectorALL_H2++;
      if(T2CutsUtil.TrkAlsoInQuarter((*TrkCit),dummy_3))  
	count_t2trackVectorALL_H3++;  
    }
	 
    
    double trkphi=0.;
    double trketa=0.;
    trketa= (*TrkCit).Eta();
    TrketaALL->Fill(trketa);
    //   std::cout<<(*TrkCit).Phi()<<std::endl;
    trkphi=(*TrkCit).GetHitT2(0).GetHitPhi();//(*TrkCit).Phi()*180/3.14159265; 
    //  std::cout<<trkphi<<std::endl;
    
    TrkphiALL->Fill(trkphi);
    NumhitinTrackALL->Fill((*TrkCit).GetHitEntries());
    
    TrkQuarterIdALL->Fill(symb);
    unsigned int hindex=symb/10;
    unsigned int hindex2=0;
    bool trkinrighhalf=true;
    
    for(unsigned int p=0;p<(*TrkCit).GetHitEntries();p++){
      
      if(PlaneInRightHalf((*TrkCit).GetHitT2(p).GetHitDetRawId())==false)
	trkinrighhalf=false;
      else{
	listOfTrackingplanesinWantedQuarter.push_back(symb);
      }
    }

    if(trkinrighhalf)
      {
	trkcounterallInQuarterAnalyzed++;
	if(T2CutsUtil_ForChi.ChiCutCond((*TrkCit)))
	  for(unsigned int p=0;p<(*TrkCit).GetHitEntries();p++)
	    {
	      hindex2=RawtoSymb(((*TrkCit).GetHitT2(p).GetHitDetRawId()));	    
	      TrkHitR_vsplane[(hindex2%10)]->Fill((*TrkCit).GetHitT2(p).GetHitR());
	      
	    }


	if(((*TrkCit).GetHitEntries()>=NumHitGood)&&((*TrkCit).GetHitEntries()<=10))	
	  trkcountergd++;
	
      }
 
    HalfTeleTrkRadiographyXY[hindex]->Fill((*TrkCit).GetHitT2(0).GetHitX(),(*TrkCit).GetHitT2(0).GetHitY());
    
    trkcounterallEVERYQUARTER++;

  }
  
  // if(verbosity)
  //std::cout<<"#TrkReco: "<<t2trackVectorALL.size()<<" "<<std::endl;
 

  NumTrackALL_ONEQuarter->Fill(trkcounterallInQuarterAnalyzed);
  NumTrackALL_EveryQuarter->Fill(trkcounterallEVERYQUARTER);

  Count_t2trackVectorALL_H0->Fill(count_t2trackVectorALL_H0); 
  Count_t2trackVectorALL_H1->Fill(count_t2trackVectorALL_H1); 
  Count_t2trackVectorALL_H2->Fill(count_t2trackVectorALL_H2); 
  Count_t2trackVectorALL_H3->Fill(count_t2trackVectorALL_H3);


  if((trkInminus==0)&&(trkInplus==0))
    EventHemisphere->Fill(0);
  if((trkInminus==0)&&(trkInplus>0))
    EventHemisphere->Fill(1);
  if((trkInminus>0)&&(trkInplus==0))
    EventHemisphere->Fill(-1);
  if((trkInminus>0)&&(trkInplus>0))
    EventHemisphere->Fill(2);


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



  if((trkInplus>0)&&(trkInminus>0))
  for(unsigned int i=0;i<10;i++){
    PadCluMultiplVsPlane->Fill((i+30),ActivePlH3[i]);
    PadCluMultiplVsPlane->Fill((i+20),ActivePlH2[i]);
    PadCluMultiplVsPlane->Fill((i+10),ActivePlH1[i]);
    PadCluMultiplVsPlane->Fill((i),ActivePlH0[i]);    
  }
















 /*
if(numevent==263)
    {

			  
 for(T2HitCollection::const_iterator ithit = t2hitcoll->begin(); ithit != t2hitcoll->end(); ithit++)
    {
      
      uint32_t rawiddet=(*ithit).GetHitDetRawId();
      unsigned int  symb=RawtoSymb(rawiddet);  

      if(symb<10)
	if(verbosity)
	  std::cout<<" ev263 Reco Hit With numstrip "<<(*ithit).GetHitNumStrip()<<"R="<<(*ithit).GetHitR()<<" Phi="<<(*ithit).GetHitPhi()<<" In plane "<<symb<< "NumStrip-Pad: "<<(*ithit).GetHitNumStrip()<<"|"<<(*ithit).GetHitNumPad()<<std::endl;
   }
    }
 */


  if(LookToRawEvent)
    {
       Handle< RawEvent >  input;
       iEvent.getByLabel(rawEventLabel, input);
       if (!input.isValid()) return;

       for (VFATFrameCollection::Iterator fr(input->frames); !fr.IsEnd(); fr.Next()) {

	 const VFATFrame &frame = *fr.Data();
	 bool GoodVfat=true;
	 // skip unmapped S-Link positions
	 map<FramePosition, unsigned int>::const_iterator it = Map_NOTDEAD->readoutPositionToId.find(fr.Position());
	 if (it != Map_NOTDEAD->readoutPositionToId.end())
	   {	       
	     // get IDs
	     unsigned short symId = it->second;

	     
	     map<unsigned int, VFATRegisters>::const_iterator dit = Map_NOTDEAD->readoutIdToRegisters.find(symId);
	     
	     unsigned short dataId = frame.getChipID();
	     unsigned short confId = 0;
	     confId = dit->second.GetDataChipID();
	   

	     // check 12Bit IDs
	     if (confId != dataId && dit != Map_NOTDEAD->readoutIdToRegisters.end()) {
	       if(verbosity)
		 std::cout<<" Error type 1: ChibId 12Bit expected in this XML position does not match what find in the data"<<std::endl;
	       GoodVfat=false;
	     }
	 
	     // check the frame
	     if (!frame.checkFootprint()) {
	       if(verbosity)
		 std::cout<<" Error type 2: Vfat FootPrintCheck in data fails"<<std::endl;
	       GoodVfat=false;
	     }

	     if (!frame.checkCRC()) {
	       if(verbosity)
		 std::cout<<" Error type 3 Vfat CRC in data fails"<<std::endl;
	       GoodVfat=false;
	     }	     
	     
	     int plane=symId/100;
	     if(GoodVfat==false)
	       if(plane<40)
		 {
		   VfatsCorruptionFrequency[plane]->Fill(symId%100);
		   VfatsCorruptionFrequency2D->Fill(symId/100,symId%100);
		 }
	   }
       }
    }
  
  if(GoodEventForTheQuarterExcludingDead) 
    countegood++;

    

  bool acceptSimuTracks=false;
  unsigned int count0simu=0;
  unsigned int count1simu=0;

  
  if(simufile)
    {
       //::::::Take Geant local hit on T2::::::

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
      
       for(map<int, PSimHitContainer>::const_iterator hitMapItr = hitMap.begin(); hitMapItr != hitMap.end(); ++hitMapItr){

          const PSimHitContainer & planeSimHits = hitMapItr->second;
          T2DetId *theT2DetId =new T2DetId(hitMapItr->first);
	  if(theT2DetId->halfTelescope()==0)
	    {
	      count0simu++;
	      // std::cout<<" Geant Hit in Half 0"<<std::endl;
	    }
	  if(theT2DetId->halfTelescope()==1)
	    {
	      count1simu++;
	      //std::cout<<" Geant Hit in Half 1"<<std::endl;
	    }
	  if(RawtoSymb(theT2DetId->calculateRawId(theT2DetId->arm(),theT2DetId->halfTelescope(),theT2DetId->plane(),theT2DetId->planeSide()))==0)
	    {
	      //std::cout<<" qui"<<std::endl;
	      // std::cout<<"  Energy Loss:  "<<std::endl; 
	       for (unsigned int l=0; l<planeSimHits.size(); l++){

		 PSimHit myplanehits= planeSimHits[l];
		 double lmyx= myplanehits.localPosition().x();
		 double lmyy= myplanehits.localPosition().y();
		 
		 // std::cout<<myplanehits.energyLoss()<<" "<<std::endl;
		  if(numevent<=200)
		    {
		      Plside0PSimHitXD0->Fill(lmyx);
		      Plside0PSimHitYD0->Fill(lmyy);			  
		    }
	       }
	    }
	  if(RawtoSymb(theT2DetId->calculateRawId(theT2DetId->arm(),theT2DetId->halfTelescope(),theT2DetId->plane(),theT2DetId->planeSide()))==1)	
	    {
	      for (unsigned int l=0; l<planeSimHits.size(); l++){

		 PSimHit myplanehits= planeSimHits[l];
		 double lmyx= myplanehits.localPosition().x();
		 double lmyy= myplanehits.localPosition().y();
		
			
		 Plside1PSimHitXD0->Fill(lmyx);
		 Plside1PSimHitYD0->Fill(lmyy);	
			
		 
		 if(numevent<=200)
		   {
		     PSimHitXD1->Fill(lmyx);
		     PSimHitYD1->Fill(lmyy);
		   }
	      }	    
	    }
       }
       

       acceptSimuTracks=true; 
       if((count0simu==10)&&(count1simu==0))
	{
	  acceptSimuTracks=true;  //The track is in full acceptance
	 //std::cout<<" Geant Track accepted"<<std::endl;
	}
      // std::cout<<" Geant det"<<hitMap.size()<<std::endl;
    }
  
  
  if(verbosity) std::cout<<" numevent"<<numevent<<std::endl;

 
  for (unsigned int i=0;i<10;i++)
    {
      TotHITNoiseEnt.at(i)=0;
      TotSTRIPNoiseEnt.at(i)=0;
      TotPADNoiseEnt.at(i)=0;
      STRIPNoiseEnt.at(i)=0;
      PADNoiseEnt.at(i)=0;
    }
 
  

 

 
  
  bool EvtWithCloseTrk=false;
  EvtWithCloseTrk=CloseTrksInEvt(t2trackVectorALL);
  
  if(verbosity)
    std::cout<<" TrkCountergood: "<<trkcountergd<<" | EvtWithCloseTrk: "<<EvtWithCloseTrk<<std::endl;
  
  
  bool skipthisCorruptedevent=false;
  
  /*
  if(skipSelectedEvents)
    if(std::find(AllCorruptedEvents.begin(),AllCorruptedEvents.end(),numevent) != AllCorruptedEvents.end())    
      skipthisCorruptedevent=true;

  if(((skipthisCorruptedevent==false)&&(GoodEventGlobal))==false)
    if(verbosity) 
      {
	std::cout<<" FALSE CONDITION (skipthisCorruptedevent==false)&&(GoodEventGlobal)"<<std::endl;
	if(skipthisCorruptedevent==false)
	  std::cout<<" FALSE skipthisCorruptedevent"<<std::endl;
	if(GoodEventGlobal==false)
	  std::cout<<" FALSE GoodEventGlobal"<<std::endl;
      }
  
  if(((trkcountergd<=MaxTrkInProcess)&&(trkcountergd>0)&&(trkcounterall<=3))==false)
    if(verbosity)
      { 
	std::cout<<" FALSE CONDITION ((trkcountergd<=MaxTrkInProcess)&&(trkcountergd>0)&&(trkcounterall<=3))"<<std::endl;
	std::cout<<" trkcountergd="<<trkcountergd<<"  MaxTrkInProcess=="<<MaxTrkInProcess<<std::endl;
      }
  */

  bool EventUtilizedForEffi=false;
  T2DetId myT2Det;
  int myhalftele;
  int myplaneside;
  int myplane;
 

  if(EvtWithCloseTrk==true)
    {
      if(verbosity)
	std::cout<<" Event with close Trk: discarded"<<std::endl;
    }

  if(EvtWithCloseTrk==false)
  if((skipthisCorruptedevent==false)&&(GoodEventForTheQuarterExcludingDead))
    if((trkcountergd>0)&&(trkcounterallInQuarterAnalyzed<=MaxTrkInQuarter)&&(trkcounterallInQuarterAnalyzed>=MinTrkInQuarter))/*trkcountergd<=MaxTrkInProcess)&&*/ /*&&(trkcounterall<=3)*/
      if (numevent<MaxEvents) // carica la matrice
	{
   
	  //std::cout<<" Event "<<numevent<<std::endl;
	  //if(verbosity) 
      
	  EventCandidate->Fill(1);

	  for(TrkCit=trackCollection->begin(); TrkCit!=trackCollection->end(); TrkCit++)
	    {

	      std::vector<T2Hit> hitvector;	  
	      for(unsigned int mm=0;mm<(*TrkCit).GetHitEntries();mm++)
		hitvector.push_back((*TrkCit).GetHitT2(mm));	  
	      
	      T2SelectionCutUtils T2CutsUtil;
	      T1T2Track trk2rz=T2CutsUtil.TrackFromHits(false,hitvector);//RZFIT  

	      TMath::Prob(trk2rz.ChiSquaredR(),(trk2rz.GetHitEntries()-2));
	      TMath::Prob(trk2rz.ChiSquaredPhi(),(trk2rz.GetHitEntries()-1));  
	  


	  if(verbosity)
	    std::cout<<"Working with Event "<<numevent<<" with "<<trkcountergd<<" tracks"<<std::endl;
	  /*
	  if(PlaneInRightHalf((*TrkCit).GetHitT2(0).GetHitDetRawId())==false)
	    std::cout<<"Trk not in quarter"<<std::endl;
	  */


	  bool thistrkinrighhalf=true;
	  for(unsigned int p=0;p<(*TrkCit).GetHitEntries();p++)
	    if(PlaneInRightHalf((*TrkCit).GetHitT2(p).GetHitDetRawId())==false)
	      thistrkinrighhalf=false;



	  // if(TrkfromRawVtx==true)	   
	  if(thistrkinrighhalf)
	      {
		bool inserted=false;
		
		
	
		if(simufile)     //controllo per simulazione non voglio tracce che non sarebbero triggerabili!
		  {
		    if(acceptSimuTracks)
		      {		
			inserted=MettiTracciaBuona(&MatriceTracce,TrkCit,numevent);
			NumhitinTrack->Fill((*TrkCit).GetHitEntries());
		      }
		  }
		else   //real data case		
		  {
		    inserted=MettiTracciaBuona(&MatriceTracce,TrkCit,numevent);
		    if(PlaneInRightHalf((*TrkCit).GetHitT2(0).GetHitDetRawId()))
		      NumhitinTrack->Fill((*TrkCit).GetHitEntries());
		  }

		if(verbosity)
		  std::cout<<" Before matrix inserted"<<std::endl;
		
		if(inserted)//analisi
		  {
		    if(verbosity)
		      std::cout<<" After matrix inserted"<<std::endl; 
		    /*if(verbosity) */
		    double trkphi=0.;
		    double trketa=0.;
		    trketa= (*TrkCit).Eta();
		    Trketa->Fill(trketa);
		    //   std::cout<<(*TrkCit).Phi()<<std::endl;
		    trkphi=(*TrkCit).Phi()*180/3.14159265; 
		    //  std::cout<<trkphi<<std::endl;
		    Trkphi->Fill(trkphi);
		    NumhitinTrackGood->Fill((*TrkCit).GetHitEntries());
		    
		    TrackInfo TrackInMatrixRow;
		    
		    //TrackInfo::IdandHit TrackInMatrixRow2;
		    TrackInMatrixRow=MatriceTracce.at(matrixentries);
		    
		    matrixentries++;
		    
		    std::vector<TrackInfo::IdandHit> theidactive=TrackInMatrixRow.idactive;
		    
		    
		    //std::cout<<TrackInMatrixRow.aphi<<" - "<<TrackInMatrixRow.bphi<<"  Tot="<<TrackInMatrixRow.bphi+TrackInMatrixRow.aphi*14000.0<<std::endl;
		    //Have sense only with RZ reco so I don't include it by default
		    //Trkphigood->Fill(TrackInMatrixRow.bphi+TrackInMatrixRow.aphi*14000.0/*180.0/3.14159*/);
		    //PolarAngles->Fill(TrackInMatrixRow.ar);
		    
		    //  std::cout<<" Evento "<<TrackInMatrixRow.eventnumber<<":"<<std::endl;	    
		    // std::cout<<TrackInMatrixRow.aphi<<" - "<<TrackInMatrixRow.bphi<<std::endl;
		    if(verbosity) std::cout<<"  TrackInMatrixRow.goodhitnumber-Effgoodhitnumber:  "<<TrackInMatrixRow.goodhitnumber<<"-"<<Effgoodhitnumber<<std::endl;
		    
		    
	   
	   
		    if(TrackInMatrixRow.goodhitnumber>=Effgoodhitnumber)
		      {
			NumhitinTrackAligned->Fill(TrackInMatrixRow.goodhitnumber);
			
			
			//--------------------------------------------
			// Part for the alignment
			//--------------------------------------------
			
			
			// ALIGNMENT PART ....carico il vettore delle roads
			std::vector<T2Hit> alignmvect;
			for(unsigned int m=0;m<theidactive.size();m++)
			  {
			    if(theidactive.at(m).iddet==1)
			      {
				
				alignmvect.push_back(theidactive.at(m).thehit);
				uint32_t rawiddet=theidactive.at(m).thehit.GetHitDetRawId();
				AllDetId->Fill((RawtoSymb(rawiddet)%10));
				
				if(std::find(alldetid.begin(),alldetid.end(),rawiddet)==alldetid.end())
				  {
				    
				    alldetid.push_back(rawiddet);			  			  
				  }
			      }
			  }
			
			/*

			//alignmvect.size()>Effgoodhitnumber ??sempre soddisfatta
			ProbXHisto->Fill(chi2X(alignmvect,UseJointProb).at(0));
			ProbYHisto->Fill(chi2Y(alignmvect,UseJointProb).at(0));
			TrkChi2X->Fill(chi2X(alignmvect,UseJointProb).at(1));
			TrkChi2Y->Fill(chi2Y(alignmvect,UseJointProb).at(1));
			
			TMatrixD covmat(4,4);
			covmat.Zero();
			double chi2corr;		
			std::vector<double> trkparam=MyLinearfitCorr(alignmvect,covmat,chi2corr);
			

			if(useRZforResol==2)
			  trkparam=MyLinearfit(alignmvect,UseJointProb);
			
		 
			//std::vector<double> trkparam=MyLinearfit(alignmvect,UseJointProb);		
			
			AXError->Fill(trkparam.at(4));
			AXHisto->Fill(trkparam.at(0));
			AYHisto->Fill(trkparam.at(2));
			BXError->Fill(trkparam.at(5));
			AYError->Fill(trkparam.at(6));
			BYError->Fill(trkparam.at(7));
			
			
			if(alignmvect.at(0).GetHitR()<AlignmentHitRMax)
			  if(alignmvect.size()>=HitNumb4Align)
			    {
			      roadXfit.push_back(alignmvect);
			      roadYfit.push_back(alignmvect);
			    }
			
			std::vector<double> drphi4Alignment;
			
	
			T2Hit hpl[10];
			
			bool pres[10];
			pres[0]=false; pres[1]=false;  pres[2]=false;  pres[3]=false;  pres[4]=false; 
			pres[5]=false;  pres[6]=false;  pres[7]=false;  pres[8]=false;  pres[9]=false;
		 
			if(alignmvect.size()>=(HitNumb4Align-2))
			  for (unsigned int u=0;u<alignmvect.size();u++)
			    {
			      unsigned int  symb=alignmvect.at(u).GetHitPlane()*2+alignmvect.at(u).GetHitPlaneSide();
			      
			      if(alignmvect.at(u).GetHitNumStrip()<=4)
				if(alignmvect.at(u).GetHitNumPad()<=3)
				  {
				    if (symb==0)
				      {
					hpl[0]=alignmvect.at(u);
					pres[0]=true;
				      }
				    if (symb==1)
				      {
					hpl[1]=alignmvect.at(u);
					pres[1]=true;
				      }				    
				    if (symb==2)
				      {
					hpl[2]=alignmvect.at(u);
					pres[2]=true;
				      }	
				    if (symb==3)
				      {
					hpl[3]=alignmvect.at(u);
					pres[3]=true;
				      }
				    if (symb==4)
				      {
					hpl[4]=alignmvect.at(u);
					pres[4]=true;
				      }				    
				    if (symb==5)
				      {
					hpl[5]=alignmvect.at(u);
					pres[5]=true;
				      }	
				    if (symb==6)
				      {
					hpl[6]=alignmvect.at(u);
					pres[6]=true;	
				      }
				    if (symb==7)
				      {
					hpl[7]=alignmvect.at(u);
					pres[7]=true;
				      }				    
				    if (symb==8)
				      {
					hpl[8]=alignmvect.at(u);
					pres[8]=true;
				      }	
				    if (symb==9)
				      {
					hpl[9]=alignmvect.at(u);
					pres[9]=true;	
				      }
				    
				  }
			    }
			
			if(alignmvect.size()>=(HitNumb4Align-2))
			  for (unsigned int u=0;u<alignmvect.size();u++)
			    {
			      drphi4Alignment=ResiduiRPhi(trkparam,alignmvect.at(u).GetHitX(), alignmvect.at(u).GetHitY(),alignmvect.at(u).GetHitZ());
			      unsigned int  symb=alignmvect.at(u).GetHitPlane()*2+alignmvect.at(u).GetHitPlaneSide();
			      if(alignmvect.at(u).GetHitNumStrip()<=4)
				if(alignmvect.at(u).GetHitNumPad()<=3)
				  {
				    if(alignmvect.at(u).GetHitR()<AlignmentHitRMax)
				      if((fabs(trkparam[0]<0.03))&&(fabs(trkparam[2]<0.03)))
					{
					  //DXAlignDet[symb]->Fill(drphi4Alignment.at(0));
					  //DYAlignDet[symb]->Fill(drphi4Alignment.at(1));
					  if(symb>0)
					    if((pres[symb])&&(pres[symb-1]))
					      {
						DXResp0[symb]->Fill(hpl[symb].GetHitX()-hpl[symb-1].GetHitX());
						DYResp0[symb]->Fill(hpl[symb].GetHitY()-hpl[symb-1].GetHitY());
					      }
					  if(symb==0)
					    {
					      DXResp0[symb]->Fill(0.);
					      DYResp0[symb]->Fill(0.);
					    }
					  if(pres[9])
					    {
					      DXResp9[symb]->Fill(hpl[9].GetHitX()-alignmvect.at(u).GetHitX());
					      DYResp9[symb]->Fill(hpl[9].GetHitY()-alignmvect.at(u).GetHitY());
					    }
					  double expX=(trkparam[0]*alignmvect.at(u).GetHitZ()+trkparam[1]);
					  double expY=(trkparam[2]*alignmvect.at(u).GetHitZ()+trkparam[3]);  
					  double expPhi=atan(fabs(expY)/fabs(expX));
					  expPhi=expPhi*180.0/3.14159;
					  if(expY<0)
					    expPhi=360.0-expPhi;
					  double expR= sqrt(expX*expX+expY*expY);
					  
					  double measuredphi=alignmvect.at(u).GetHitPhi();
					  double mesuredR=alignmvect.at(u).GetHitR();  
					  
					  double measuredphi2=atan(fabs(alignmvect.at(u).GetHitY())/fabs(alignmvect.at(u).GetHitX()));
					  measuredphi2=measuredphi2*180.0/3.14159;
					  if(alignmvect.at(u).GetHitY()<0)
					    measuredphi2=360.0-measuredphi2;

					  double mesuredR2=sqrt(alignmvect.at(u).GetHitX()*alignmvect.at(u).GetHitX()+alignmvect.at(u).GetHitY()*alignmvect.at(u).GetHitY());
				   
					  
				
					  
					}
				  }
			    }
			

			
			
			*/


			
			if(verbosity) std::cout<<" theidactive.size()= "<<theidactive.size()<<std::endl;
			for(unsigned int l=0;l<theidactive.size();l++)
			  {
			    //std::cout<<" So qui1"<<std::endl;
			    unsigned int countref=0;
			    std::vector<T2Hit> refhitv;
			    //std::vector<T2Hit> allhits;

			    if(theidactive.at(l).iddet==1)
			      {
				
				HitMatrixR->Fill(theidactive.at(l).thehit.GetHitR());
				HitMatrixPhi->Fill(theidactive.at(l).thehit.GetHitPhi());
			      }
			    
			    for(unsigned int m=0;m<theidactive.size();m++)
			      {		      
				
				if(m!=l)
				  {
				    if(theidactive.at(m).iddet==1)
				      {
					countref++;
					refhitv.push_back(theidactive.at(m).thehit);
					
					if((theidactive.at(m).thehit.GetHitNumStrip()>513)||(theidactive.at(m).thehit.GetHitNumPad()>1600))
					  std::cout<<" ---->: Here THERE IS A BUG IN THE HIT!!"<<std::endl;
					if((theidactive.at(m).thehit.GetHitNumStrip()==0)||(theidactive.at(m).thehit.GetHitNumPad()==0))
					  std::cout<<" ---->: Here THERE IS A BUG IN THE HIT!!"<<std::endl;
				      }
				    
				  }
				//allhits.push_back(theidactive.at(m).thehit);
			      }
			  
			    
			    if(verbosity) 
			      std::cout<<" countref-Effgoodhitnumber"<<countref<<"-"<<Effgoodhitnumber<<std::endl;
			    if(countref>=Effgoodhitnumber)
			      {
				EventUtilizedForEffi=true;
				//std::cout<<" Debug-0"<<std::endl;
				// if(alreadyprint0==false)
				//{
				TrackPhi0R0->Fill((*TrkCit).GetHitT2(0).GetHitPhi(),(*TrkCit).GetHitT2(0).GetHitR());
				
				if(verbosity)
				  std::cout<<" @@@@@ >>>> One Track utilized with phi "<<(*TrkCit).GetHitT2(0).GetHitPhi()<<" and R0: "<<(*TrkCit).GetHitT2(0).GetHitR()<<" Utilized for Efficiency in plane "<<l<<std::endl;
				//  alreadyprint0=true;
				//}
				unsigned int originalTrkSize=TrkCit->GetHitEntries();
				unsigned int symb;
				for(unsigned mm=0;mm<originalTrkSize;mm++)
				  {			  
				    symb=TrkCit->GetHitT2(mm).GetHitDetRawId();
				    symb=RawtoSymb(symb);
				    RelEffi_ForTracking->Fill(symb%10);
				  }
				
				bool stripon=false;
				bool padon=false;
				bool hiton=false;
				bool stripongood=false;
				bool padongood=false;
				bool hitongood=false;
				std::vector<double> paramvect;
				std::vector<double> OLDparamvect;
				std::vector<double> paramvectForAlign;
				bool hitAlreadysaved=false;
				bool padAlreadysaved=false;
				bool stripAlreadysaved=false;
				
				T2Hit StripHitgood;
				T2Hit PadHitgood;
				T2Hit StripAndPadHitgood;
				
				TMatrixD covmat(4,4);
				covmat.Zero();
				double chi2corr;
				// paramvectForAlign=MyLinearfit(refhitv,UseJointProb); // exlude the tested chamber from getting parameter
				paramvectForAlign=MyLinearfitCorr(refhitv,covmat,chi2corr);
				if(useRZforResol==2)
				  paramvectForAlign=MyLinearfit(refhitv,UseJointProb);
		      
				
				//std::cout<<" Param V-Al:"<<std::endl;		    
				//std::cout<<paramvectForAlign.at(0)<<" - "<<paramvectForAlign.at(1)<<" - "<<paramvectForAlign.at(2)<<" - "<<paramvectForAlign.at(3)<<" - "<<std::endl;
				
				paramvect.push_back(TrackInMatrixRow.ar);
				paramvect.push_back(TrackInMatrixRow.br);
				paramvect.push_back(TrackInMatrixRow.aphi);
				paramvect.push_back(TrackInMatrixRow.bphi);
		      
				
				
				TracksplanestatN->Fill(l);
				bool lhitindet=false;
				lhitindet=HitIndetector(paramvect,l);
				
		      
				std::vector<T2Hit> firstHit;
				std::vector<T2Hit> secondHit;
				std::vector<unsigned int> planeindexfirst;
				std::vector<unsigned int> planeindexsecond;
		      
				bool NoisyplaneForPad=false;
				bool NoisyplaneForStrip=false;
				
				double RVal=(refhitv.at(0).GetHitR()+refhitv.at(refhitv.size()-1).GetHitR())/2.0;
				unsigned int Rsector=0;
				if(RVal<50.)
				  {
				    Rsector=0;
				  }
				else
				  {
				    if(RVal<92.5)
				      {
					Rsector=1;
				      }
				    else
				      {
					if(RVal<121.5)
					  Rsector=2;
					else
					  Rsector=3;
				      }
				  }
				
				bool RsectorDead=false;
				RsectorDead=HitInDeadSector(Rsector,l);

				if((lhitindet==true))				 
				  {
				    
				    if(alreadyprint==false)
				      {
					if(verbosity)
					  std::cout<<" ->->->->->->->->      CANDIDATE EVENT ("<<numevent<<") FOR EFF CALCULATION:"<<std::endl;	 
					alreadyprint=true;
				      }
				    
				    
				    
				    for(T2HitCollection::const_iterator ithit = t2hitcoll->begin(); ithit != t2hitcoll->end(); ithit++)
				      {
					
					//std::cout<<" -----"<<std::endl;	 
					unsigned int  symb=(*ithit).GetHitPlane()*2+(*ithit).GetHitPlaneSide();
					unsigned int quarter=(*ithit).GetHitArm()*2+(*ithit).GetHitHalftele();
					
					
					

					if((symb==l)&&(quarter==SelectedHalf))
					  {
					    
					    RPhiEvtHit[l]->Fill((*ithit).GetHitPhi(),(*ithit).GetHitR());
					    //std::cout<<" Debug-1.2"<<std::endl;
					    if((*ithit).GetHitNumStrip()>0)
					      {
						
						//if(numevent==300)
						//  std::cout<<" ((*ithit).GetHitNumStrip()>0)"<<std::endl;
						
						TotSTRIPNoiseEnt.at(l)++;  
						
						stripon=true;
						//std::vector<double> drphi=ResiduiRPhi(paramvect,(*ithit).GetHitX(), (*ithit).GetHitY(),(*ithit).GetHitZ());
						//std::cout<<" Debug-1.5"<<std::endl;
						std::vector<double> drphi=ResiduiForStrip(refhitv,(*ithit));
						//if(numevent==300)
				//std::cout<<" Debug-2"<<std::endl;
				//	std::cout<<" !!!!!!!!!!!!!! !!!!!!!! drphi.at(1)="<<drphi.at(1)<<" plane="<<l<<std::endl;
						
						if(fabs(drphi.at(1))<500.)//Was 5000
						  {
						    
						    
						    if(l==0)
						      {
							
							
							EFF_DRstripNoCutdet0->Fill(drphi.at(0));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSStripdet0NoCut->Fill((*ithit).GetHitNumStrip());
							
							if((*ithit).GetHitNumStrip()<=EffMaxStrip)
							  if(fabs(drphi.at(0))<Effmaxdrhit)
							    {					    
							      stripongood=true;
							      EFF_DRstripdet0->Fill(drphi.at(0));	
							      EFF_CLSStripdet0->Fill((*ithit).GetHitNumStrip());
							      // std::cout<<" OK! numevent: "<<numevent<<" CLS strip found for plane 0 at R="<<(*ithit).GetHitR()<<" Phi="<<(*ithit).GetHitPhi()<<"(expected act R="<<aaa.GetHitR()<<" Phi="<<aaa.GetHitPhi()<<")"<<std::endl;
							    }
							
							if(stripongood==false)
							  {					
							    //		std::cout<<" NO: numevent: "<<numevent<<" CLS strip found for plane 0 at R="<<(*ithit).GetHitR()<<" and Phi="<<(*ithit).GetHitPhi()<<" but expected at R="<<aaa.GetHitR()<<" Phi="<<aaa.GetHitPhi()<<std::endl;
							  }
							
						      }
						    if(l==2)
						      {
							EFF_DRstripNoCutdet2->Fill(drphi.at(0));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSStripdet2NoCut->Fill((*ithit).GetHitNumStrip());
							if((*ithit).GetHitNumStrip()<=EffMaxStrip)				      
							  if(fabs(drphi.at(0))<Effmaxdrhit)
							    { 
							      stripongood=true;	
							      EFF_DRstripdet2->Fill(drphi.at(0));	
							      EFF_CLSStripdet2->Fill((*ithit).GetHitNumStrip());
							    }
						      }
						    if(l==4)
						      {
							EFF_DRstripNoCutdet4->Fill(drphi.at(0));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSStripdet4NoCut->Fill((*ithit).GetHitNumStrip());
							if((*ithit).GetHitNumStrip()<=EffMaxStrip)				      
							  if(fabs(drphi.at(0))<Effmaxdrhit)
							    { 
							      stripongood=true;	 
							      EFF_DRstripdet4->Fill(drphi.at(0));	
							      EFF_CLSStripdet4->Fill((*ithit).GetHitNumStrip());
							    }
						      }
						    if(l==6)
						      {
							EFF_DRstripNoCutdet6->Fill(drphi.at(0));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSStripdet6NoCut->Fill((*ithit).GetHitNumStrip());
							if((*ithit).GetHitNumStrip()<=EffMaxStrip)
							  if(fabs(drphi.at(0))<Effmaxdrhit)
							    {
							      stripongood=true;
							      EFF_DRstripdet6->Fill(drphi.at(0));	
							      EFF_CLSStripdet6->Fill((*ithit).GetHitNumStrip());
							    }	  
						      }
						    if(l==8)
						      {
							EFF_DRstripNoCutdet8->Fill(drphi.at(0));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSStripdet8NoCut->Fill((*ithit).GetHitNumStrip());				    
							if((*ithit).GetHitNumStrip()<=EffMaxStrip)
							  if(fabs(drphi.at(0))<Effmaxdrhit)
							    { 
							      stripongood=true;	 
							      EFF_DRstripdet8->Fill(drphi.at(0));
							      EFF_CLSStripdet8->Fill((*ithit).GetHitNumStrip());
							    }

							  //else
							    //{
							      //std::cout<<drphi.at(0)<<" "<<Effmaxdrhit<<std::endl;
							    //}
						      }
						    
						    
						    if(l==1)
						      {
							
							T2Hit aaa=GiveExtrapolatedHit(1,refhitv);
							
							EFF_DRstripNoCutdet1->Fill(drphi.at(0));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSStripdet1NoCut->Fill((*ithit).GetHitNumStrip());
							
						
							if((*ithit).GetHitNumStrip()<=EffMaxStrip)	
							  if(fabs(drphi.at(0))<Effmaxdrhit)
							  {
							      stripongood=true;
							      EFF_DRstripdet1->Fill(drphi.at(0));	
							      EFF_CLSStripdet1->Fill((*ithit).GetHitNumStrip());
							      if(verbosity)
								std::cout<<" OK! numevent: "<<numevent<<" CLS strip found for plane 1 at R="<<(*ithit).GetHitR()<<" Phi="<<(*ithit).GetHitPhi()<<"(expected act R="<<aaa.GetHitR()<<" Phi="<<aaa.GetHitPhi()<<")"<<std::endl;
							    }
							if(stripongood==false)
							  {		
							    if(verbosity)			
							      std::cout<<" NO: numevent: "<<numevent<<" CLS strip found for plane 1 at R="<<(*ithit).GetHitR()<<" and Phi="<<(*ithit).GetHitPhi()<<" but expected at R="<<aaa.GetHitR()<<" Phi="<<aaa.GetHitPhi()<<std::endl;
							  }
							
						      }

						    if(l==3)
						      {
							
							EFF_DRstripNoCutdet3->Fill(drphi.at(0));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSStripdet3NoCut->Fill((*ithit).GetHitNumStrip());

							if((*ithit).GetHitNumStrip()<=EffMaxStrip)				      
							  if(fabs(drphi.at(0))<Effmaxdrhit)
							    { 
							      stripongood=true;	
							      EFF_DRstripdet3->Fill(drphi.at(0));	
							      EFF_CLSStripdet3->Fill((*ithit).GetHitNumStrip());
							    }
						      }
						    if(l==5)
						      {
							EFF_DRstripNoCutdet5->Fill(drphi.at(0));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSStripdet5NoCut->Fill((*ithit).GetHitNumStrip());
							
							if((*ithit).GetHitNumStrip()<=EffMaxStrip)				      
							  if(fabs(drphi.at(0))<Effmaxdrhit) 
							    { 
							      stripongood=true;	 
							      EFF_DRstripdet5->Fill(drphi.at(0));	
							      EFF_CLSStripdet5->Fill((*ithit).GetHitNumStrip());
							    }
						      }
						    if(l==7)
						      {
							EFF_DRstripNoCutdet7->Fill(drphi.at(0));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSStripdet7NoCut->Fill((*ithit).GetHitNumStrip());
							if((*ithit).GetHitNumStrip()<=EffMaxStrip)		  
							  if(fabs(drphi.at(0))<Effmaxdrhit)
							    {
							      stripongood=true;
							      EFF_DRstripdet7->Fill(drphi.at(0));	
							      EFF_CLSStripdet7->Fill((*ithit).GetHitNumStrip());
							    }	  
						      }
						    if(l==9)
						      {
							EFF_DRstripNoCutdet9->Fill(drphi.at(0));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSStripdet9NoCut->Fill((*ithit).GetHitNumStrip());
							
							if((*ithit).GetHitNumStrip()<=EffMaxStrip)				      				      				      if(fabs(drphi.at(0))<Effmaxdrhit)
							  { 
							    stripongood=true;	 
							    EFF_DRstripdet9->Fill(drphi.at(0));
							    EFF_CLSStripdet9->Fill((*ithit).GetHitNumStrip());
							  }
						      }
						  }//if fabs dphi<50.
						//	else 
						//std::cout<<" drphi.at(1):"<<drphi.at(1)<<std::endl;
						
						if(stripongood){
						  CumulativeStripCluSize_UsedInEffi_HX->Fill((*ithit).GetHitNumStrip()); 
						}
						
					
						
					      }//if((*ithit).GetHitNumStrip()>0)
					    
			    

					    if((*ithit).GetHitNumPad()>0)
					      {
						padon=true;
						//std::vector<double> drphi=ResiduiRPhi(paramvect,(*ithit).GetHitX(), (*ithit).GetHitY(),(*ithit).GetHitZ());
						//(std::vector<T2Hit>refhitv, T2Hit hit,double &expPadXErr,double &expPadYErr)
						double avgdx=0.; double avgdy=0.;
						std::vector<double> drphi=ResiduiForPad(refhitv,(*ithit),avgdx,avgdy);
						
						TotPADNoiseEnt.at(l)++;
						
						bool adjacentpad=false;
						adjacentpad=PadIsNeighbour(refhitv,(*ithit));
					
						if(((fabs(drphi.at(0))<3*avgdx)&&(fabs(drphi.at(1))<3*avgdy))||(adjacentpad))//Put at 2 as def
						  {
						    
						    if(l==0)
						      {
							// T2Hit aaa=GiveExtrapolatedHit(0,refhitv);
							EFF_DPhipadNoCutdet0->Fill(drphi.at(1));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSPaddet0NoCut->Fill((*ithit).GetHitNumPad());
							if((*ithit).GetHitNumPad()<=EffMaxPad)
							  //if(fabs(drphi.at(1))<Effmaxdphihit)
							  {
							    padongood=true;
							    EFF_DPhipaddet0->Fill(drphi.at(1));	
							    EFF_CLSPaddet0->Fill((*ithit).GetHitNumPad());
							    
							    //   std::cout<<" OK! numevent: "<<numevent<<" CLS PAD found for plane 0 at R="<<(*ithit).GetHitR()<<" Phi="<<(*ithit).GetHitPhi()<<"(expected act R="<<aaa.GetHitR()<<" Phi="<<aaa.GetHitPhi()<<")"<<std::endl;
							  }
							//
							// if(padongood==false)
							//   {					
							//	std::cout<<" NO: numevent: "<<numevent<<" CLS PAD found for plane 0 at R="<<(*ithit).GetHitR()<<" and Phi="<<(*ithit).GetHitPhi()<<" but expected at R="<<aaa.GetHitR()<<" Phi="<<aaa.GetHitPhi()<<std::endl;
							//  }
							
						      }
						    if(l==2)
						      {
							EFF_DPhipadNoCutdet2->Fill(drphi.at(1));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSPaddet2NoCut->Fill((*ithit).GetHitNumPad());
							if((*ithit).GetHitNumPad()<=EffMaxPad)
							  //if(fabs(drphi.at(1))<Effmaxdphihit)
							  {
							    padongood=true;
							    EFF_DPhipaddet2->Fill(drphi.at(1));	
							    EFF_CLSPaddet2->Fill((*ithit).GetHitNumPad());
							  }
						      }
						    if(l==4)
						      {
							EFF_DPhipadNoCutdet4->Fill(drphi.at(1));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSPaddet4NoCut->Fill((*ithit).GetHitNumPad());
							if((*ithit).GetHitNumPad()<=EffMaxPad)
							  //if(fabs(drphi.at(1))<Effmaxdphihit)
							  {
							    padongood=true;
							    EFF_DPhipaddet4->Fill(drphi.at(1));	
							    EFF_CLSPaddet4->Fill((*ithit).GetHitNumPad());
							  }
							
						      }
						    if(l==6)
						      {
							EFF_DPhipadNoCutdet6->Fill(drphi.at(1));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSPaddet6NoCut->Fill((*ithit).GetHitNumPad());
							if((*ithit).GetHitNumPad()<=EffMaxPad)
							  //if(fabs(drphi.at(1))<Effmaxdphihit)
							  {
							    padongood=true;
							    EFF_DPhipaddet6->Fill(drphi.at(1));
							    EFF_CLSPaddet6->Fill((*ithit).GetHitNumPad());
							  }
						      }
						    if(l==8)
						      {
							EFF_DPhipadNoCutdet8->Fill(drphi.at(1));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSPaddet8NoCut->Fill((*ithit).GetHitNumPad());
							if((*ithit).GetHitNumPad()<=EffMaxPad)
							  //if(fabs(drphi.at(1))<Effmaxdphihit)
							  {
							    padongood=true;
							    EFF_DPhipaddet8->Fill(drphi.at(1));	
							    EFF_CLSPaddet8->Fill((*ithit).GetHitNumPad());
							  }
						      }
						    
						    
						    
						    if(l==1)
						      {
							T2Hit aaa=GiveExtrapolatedHit(1,refhitv);
							EFF_DPhipadNoCutdet1->Fill(drphi.at(1));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSPaddet1NoCut->Fill((*ithit).GetHitNumPad());
							if((*ithit).GetHitNumPad()<=EffMaxPad)
							  //if(fabs(drphi.at(1))<Effmaxdphihit)
							  {
							    padongood=true;
							    EFF_DPhipaddet1->Fill(drphi.at(1));	
							    EFF_CLSPaddet1->Fill((*ithit).GetHitNumPad());
							    if(verbosity)
							      std::cout<<" OK! numevent: "<<numevent<<" CLS PAD found for plane 1 at R="<<(*ithit).GetHitR()<<" Phi="<<(*ithit).GetHitPhi()<<"(expected act R="<<aaa.GetHitR()<<" Phi="<<aaa.GetHitPhi()<<")   [drphi.at(1)= "<<(drphi.at(1)) <<"]"<<std::endl;
							  }
							
							if(padongood==false)
							  {	
							    if(verbosity)
							      std::cout<<" NO: numevent: "<<numevent<<" CLS PAD found for plane 1 at R="<<(*ithit).GetHitR()<<" and Phi="<<(*ithit).GetHitPhi()<<" but expected at R="<<aaa.GetHitR()<<" Phi="<<aaa.GetHitPhi()<<std::endl;
							  }
							
						      }
						    if(l==3)
						      {
							EFF_DPhipadNoCutdet3->Fill(drphi.at(1));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSPaddet3NoCut->Fill((*ithit).GetHitNumPad());
							if((*ithit).GetHitNumPad()<=EffMaxPad)
							  //if(fabs(drphi.at(1))<Effmaxdphihit)
							  {
							    padongood=true;
							    EFF_DPhipaddet3->Fill(drphi.at(1));	
							    EFF_CLSPaddet3->Fill((*ithit).GetHitNumPad());
							    
							  }
						      }
						    if(l==5)
						      {
							EFF_DPhipadNoCutdet5->Fill(drphi.at(1));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSPaddet5NoCut->Fill((*ithit).GetHitNumPad());
							if((*ithit).GetHitNumPad()<=EffMaxPad)
							  //if(fabs(drphi.at(1))<Effmaxdphihit)
							  {
							    padongood=true;
							    EFF_DPhipaddet5->Fill(drphi.at(1));	
							    EFF_CLSPaddet5->Fill((*ithit).GetHitNumPad());
							  }
							
						      }
						    if(l==7)
						      {
							EFF_DPhipadNoCutdet7->Fill(drphi.at(1));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSPaddet7NoCut->Fill((*ithit).GetHitNumPad());
							if((*ithit).GetHitNumPad()<=EffMaxPad)
							  //if(fabs(drphi.at(1))<Effmaxdphihit)
							  {
							    padongood=true;
							    EFF_DPhipaddet7->Fill(drphi.at(1));
							    EFF_CLSPaddet7->Fill((*ithit).GetHitNumPad());
							  }
						      }
						    if(l==9)
						      {
							EFF_DPhipadNoCutdet9->Fill(drphi.at(1));
							if(fabs(drphi.at(0))<(Effmaxdrhit*2))
							  EFF_CLSPaddet9NoCut->Fill((*ithit).GetHitNumPad());
							if((*ithit).GetHitNumPad()<=EffMaxPad)
							  //if(fabs(drphi.at(1))<Effmaxdphihit)
							  {
							    padongood=true;
							    EFF_DPhipaddet9->Fill(drphi.at(1));	
							    EFF_CLSPaddet9->Fill((*ithit).GetHitNumPad());
							  }
						      }
						    
						  }
						
						if(padongood)
						  CumulativePadCluSize_UsedInEffi_HX->Fill((*ithit).GetHitNumPad());
						else
						  if(verbosity)
						    std::cout<<"Pad Hit Fail. DX-DY:"<<drphi.at(0)<<" "<<drphi.at(1)<<" HIT xy:"<<(*ithit).GetHitX()<<" "<< (*ithit).GetHitY()<<" Avg xy:"<<avgdx<<" "<<avgdy<<" hit dx:"<<(*ithit).GetHitDX()<<std::endl;//
					      }
					    
					    
					    
					    if(((*ithit).GetHitNumPad()>0)&&((*ithit).GetHitNumStrip()>0))
					      {
						hiton=true;

						//std::vector<double> drphi=ResiduiRPhi(paramvect,(*ithit).GetHitX(), (*ithit).GetHitY(),(*ithit).GetHitZ());
						std::vector<double> drphi4Align=ResiduiRPhi(paramvectForAlign,(*ithit).GetHitX(), (*ithit).GetHitY(),(*ithit).GetHitZ());
						TotHITNoiseEnt.at(l)++;
						
						
						std::vector<double> drphi=ResiduiForC1HIT(refhitv,(*ithit));
						
						if(fabs(drphi.at(0))<Effmaxdrhit)
						  if(fabs(drphi.at(1))<Effmaxdphihit)
						    if((*ithit).GetHitNumPad()<=EffMaxPad)
						      if((*ithit).GetHitNumStrip()<=EffMaxStrip)
							{
							  CluPadentriesGoodHit->Fill((*ithit).GetHitNumPad());
							  CluStripentriesGoodHit->Fill((*ithit).GetHitNumStrip());
							}
						
						if(l==0)
						  {
						    EFF_DPhiHitNoCutdet0->Fill(drphi.at(1));
						    EFF_DRHitNoCutdet0->Fill(drphi.at(0));
						    
						    if((*ithit).GetHitNumStrip()<=EffMaxStrip)
						      if((*ithit).GetHitNumPad()<=EffMaxPad)
							{
							  AL_DPhiHitNoCutdet0->Fill(drphi4Align.at(1));
							  AL_DRHitNoCutdet0->Fill(drphi4Align.at(0));
							  
							  if(fabs(drphi.at(0))<Effmaxdrhit)
							    if(fabs(drphi.at(1))<Effmaxdphihit)
							      {
								hitongood=true;
								EFF_DPhiHitdet0->Fill(drphi.at(1));
								EFF_DRHitdet0->Fill(drphi.at(0));	
								EFF_C1CLSStripdet0->Fill((*ithit).GetHitNumStrip());
								EFF_C1CLSPaddet0->Fill((*ithit).GetHitNumPad());
							      }
							}
						    
						  }
						if(l==2)
						  {
						    EFF_DPhiHitNoCutdet2->Fill(drphi.at(1));
						    EFF_DRHitNoCutdet2->Fill(drphi.at(0));
						    if((*ithit).GetHitNumStrip()<=EffMaxStrip)
						      if((*ithit).GetHitNumPad()<=EffMaxPad)
							{
							  AL_DPhiHitNoCutdet2->Fill(drphi4Align.at(1));
							  AL_DRHitNoCutdet2->Fill(drphi4Align.at(0));
							  
							  if(fabs(drphi.at(0))<Effmaxdrhit)
							    if(fabs(drphi.at(1))<Effmaxdphihit)
							      {
								hitongood=true;
								EFF_DPhiHitdet2->Fill(drphi.at(1));
								EFF_DRHitdet2->Fill(drphi.at(0));	
								EFF_C1CLSStripdet2->Fill((*ithit).GetHitNumStrip());
								EFF_C1CLSPaddet2->Fill((*ithit).GetHitNumPad());
							      }
							}
						  }
						if(l==4)
						  {
						    EFF_DPhiHitNoCutdet4->Fill(drphi.at(1));
						    EFF_DRHitNoCutdet4->Fill(drphi.at(0));
						    if((*ithit).GetHitNumStrip()<=EffMaxStrip)
						      if((*ithit).GetHitNumPad()<=EffMaxPad)
							{
							  AL_DPhiHitNoCutdet4->Fill(drphi4Align.at(1));
							  AL_DRHitNoCutdet4->Fill(drphi4Align.at(0));
							  
							  if(fabs(drphi.at(0))<Effmaxdrhit)
							    if(fabs(drphi.at(1))<Effmaxdphihit)
							      {
								hitongood=true;
								EFF_DPhiHitdet4->Fill(drphi.at(1));
								EFF_DRHitdet4->Fill(drphi.at(0));	
								EFF_C1CLSStripdet4->Fill((*ithit).GetHitNumStrip());
								EFF_C1CLSPaddet4->Fill((*ithit).GetHitNumPad());
							      }
							}
						  }
						if(l==6)
						  {
						    EFF_DPhiHitNoCutdet6->Fill(drphi.at(1));
						    EFF_DRHitNoCutdet6->Fill(drphi.at(0));
						    if((*ithit).GetHitNumStrip()<=EffMaxStrip)
						      if((*ithit).GetHitNumPad()<=EffMaxPad)
							{
							  AL_DPhiHitNoCutdet6->Fill(drphi4Align.at(1));
							  AL_DRHitNoCutdet6->Fill(drphi4Align.at(0));
							  
							  if(fabs(drphi.at(0))<Effmaxdrhit)
							    if(fabs(drphi.at(1))<Effmaxdphihit)
							      {
								hitongood=true;
								EFF_DPhiHitdet6->Fill(drphi.at(1));
								EFF_DRHitdet6->Fill(drphi.at(0));		
								EFF_C1CLSStripdet6->Fill((*ithit).GetHitNumStrip());
								EFF_C1CLSPaddet6->Fill((*ithit).GetHitNumPad());
							      }
							}
						  }
						if(l==8)
						  {
						    EFF_DPhiHitNoCutdet8->Fill(drphi.at(1));
						    EFF_DRHitNoCutdet8->Fill(drphi.at(0));				    
						    if((*ithit).GetHitNumStrip()<=EffMaxStrip)
						      {
							if((*ithit).GetHitNumPad()<=EffMaxPad)
							  {
							    AL_DPhiHitNoCutdet8->Fill(drphi4Align.at(1));
							    AL_DRHitNoCutdet8->Fill(drphi4Align.at(0));
							    
							    if(fabs(drphi.at(0))<Effmaxdrhit)
							      {
								if(fabs(drphi.at(1))<Effmaxdphihit)
								  {
								    hitongood=true;
								    EFF_DPhiHitdet8->Fill(drphi.at(1));
								    EFF_DRHitdet8->Fill(drphi.at(0));	
								    EFF_C1CLSStripdet8->Fill((*ithit).GetHitNumStrip());
								    EFF_C1CLSPaddet8->Fill((*ithit).GetHitNumPad());				      
								  }
								else
								  if(verbosity)
								    std::cout<<" Look Plane8 EFF: class1 DPHI "<<drphi.at(1)<<std::endl;
								if(verbosity)
								  std::cout<<" Look Plane8 EFF: class1 DR Ok Hit cls: "<<(*ithit).GetHitNumStrip()<<(*ithit).GetHitNumPad()<</*" "<<(*ithit).GetHitArm()<<(*ithit).GetHitHalftele()<<*/" DR: "<<drphi.at(0)<<"  Hit R-Phi: "<<(*ithit).GetHitR()<<" "<<(*ithit).GetHitPhi()<<std::endl;
							      }
							    else
							      if(verbosity)
								std::cout<<" class1 dr "<<(*ithit).GetHitNumStrip()<<(*ithit).GetHitNumPad()<<" "<<(*ithit).GetHitArm()<<(*ithit).GetHitHalftele()<<" "<<drphi.at(0)<<"  "<<(*ithit).GetHitR()<<" "<<(*ithit).GetHitPhi()<<std::endl;
							  }
							//else
							  //std::cout<<(*ithit).GetHitNumPad()<<std::endl;
						      }
						   // else
						     // std::cout<<(*ithit).GetHitNumStrip()<<std::endl;
						  }	
						
						if(l==1)
						  {
						    EFF_DPhiHitNoCutdet1->Fill(drphi.at(1));
						    EFF_DRHitNoCutdet1->Fill(drphi.at(0));
						    
						    if((*ithit).GetHitNumStrip()<=EffMaxStrip)
						      if((*ithit).GetHitNumPad()<=EffMaxPad)
							{
							  AL_DPhiHitNoCutdet1->Fill(drphi4Align.at(1));
							  AL_DRHitNoCutdet1->Fill(drphi4Align.at(0));
							  
							  if(fabs(drphi.at(0))<Effmaxdrhit)
							    if(fabs(drphi.at(1))<Effmaxdphihit)
							      {
								hitongood=true;
								EFF_DPhiHitdet1->Fill(drphi.at(1));
								EFF_DRHitdet1->Fill(drphi.at(0));	
								EFF_C1CLSStripdet1->Fill((*ithit).GetHitNumStrip());
								EFF_C1CLSPaddet1->Fill((*ithit).GetHitNumPad());
								
							      }
							}
						    
						  }
						if(l==3)
						  {
						    EFF_DPhiHitNoCutdet3->Fill(drphi.at(1));
						    EFF_DRHitNoCutdet3->Fill(drphi.at(0));
						    if((*ithit).GetHitNumStrip()<=EffMaxStrip)
						      if((*ithit).GetHitNumPad()<=EffMaxPad)
							{
							  AL_DPhiHitNoCutdet3->Fill(drphi4Align.at(1));
							  AL_DRHitNoCutdet3->Fill(drphi4Align.at(0));
							  if(fabs(drphi.at(0))<Effmaxdrhit)
							    if(fabs(drphi.at(1))<Effmaxdphihit)
							      {
								hitongood=true;
								EFF_DPhiHitdet3->Fill(drphi.at(1));
								EFF_DRHitdet3->Fill(drphi.at(0));	
								EFF_C1CLSStripdet3->Fill((*ithit).GetHitNumStrip());
								EFF_C1CLSPaddet3->Fill((*ithit).GetHitNumPad());
							      }
							}
						  }
						if(l==5)
						  {
						    EFF_DPhiHitNoCutdet5->Fill(drphi.at(1));
						    EFF_DRHitNoCutdet5->Fill(drphi.at(0));
						    if((*ithit).GetHitNumStrip()<=EffMaxStrip)
						      if((*ithit).GetHitNumPad()<=EffMaxPad)
							{
							  AL_DPhiHitNoCutdet5->Fill(drphi4Align.at(1));
							  AL_DRHitNoCutdet5->Fill(drphi4Align.at(0));
							  if(fabs(drphi.at(0))<Effmaxdrhit)
							    if(fabs(drphi.at(1))<Effmaxdphihit)
							      {
								hitongood=true;
								EFF_DPhiHitdet5->Fill(drphi.at(1));
								EFF_DRHitdet5->Fill(drphi.at(0));	
								EFF_C1CLSStripdet5->Fill((*ithit).GetHitNumStrip());
								EFF_C1CLSPaddet5->Fill((*ithit).GetHitNumPad());
							      }
							}
						  }
						if(l==7)
						  {
						    EFF_DPhiHitNoCutdet7->Fill(drphi.at(1));
						    EFF_DRHitNoCutdet7->Fill(drphi.at(0));
						    if((*ithit).GetHitNumStrip()<=EffMaxStrip)
						      if((*ithit).GetHitNumPad()<=EffMaxPad)
							{
							  AL_DPhiHitNoCutdet7->Fill(drphi4Align.at(1));
							  AL_DRHitNoCutdet7->Fill(drphi4Align.at(0));
							  if(fabs(drphi.at(0))<Effmaxdrhit)
							    if(fabs(drphi.at(1))<Effmaxdphihit)
							      {
								hitongood=true;
								EFF_DPhiHitdet7->Fill(drphi.at(1));
								EFF_DRHitdet7->Fill(drphi.at(0));	
								EFF_C1CLSStripdet7->Fill((*ithit).GetHitNumStrip());
								EFF_C1CLSPaddet7->Fill((*ithit).GetHitNumPad());
							      }
							}
						  }
						if(l==9)
						  {
						    EFF_DPhiHitNoCutdet9->Fill(drphi.at(1));
						    EFF_DRHitNoCutdet9->Fill(drphi.at(0));
						    if((*ithit).GetHitNumStrip()<=EffMaxStrip)
						      if((*ithit).GetHitNumPad()<=EffMaxPad)
							{
							  AL_DPhiHitNoCutdet9->Fill(drphi4Align.at(1));
							  AL_DRHitNoCutdet9->Fill(drphi4Align.at(0));
							  if(fabs(drphi.at(0))<Effmaxdrhit)
							    if(fabs(drphi.at(1))<Effmaxdphihit)
							      {
								hitongood=true;
								EFF_DPhiHitdet9->Fill(drphi.at(1));
								EFF_DRHitdet9->Fill(drphi.at(0));	
								EFF_C1CLSStripdet9->Fill((*ithit).GetHitNumStrip());
								EFF_C1CLSPaddet9->Fill((*ithit).GetHitNumPad());
							      }
							}
						  }			      
						
					      }
					    
					    //Copy the good hit of the l plane
					    if(hitongood)
					      {
						if(hitAlreadysaved==false) //hitongood remain true after the hit successfully found!!
						  {
						    StripAndPadHitgood=(*ithit);
						    hitAlreadysaved=true;
						  }
					      }
					    else
					      {
						if(padongood)
						  if(padAlreadysaved==false) 
						    {
						      PadHitgood=(*ithit);
						      padAlreadysaved=true;
						    }
						
						if(stripongood)
						  {
						    if(stripAlreadysaved==false) 
						      {
							StripHitgood=(*ithit);
							stripAlreadysaved=true;
						      }
						  }
						//	std::cout<<" ||"<<hitongood<<"   "<<padongood<<"   "<<stripongood<<std::endl;
					      }

			    
			    
					    
					  }//if symb==l
			
					//if(((TotSTRIPNoiseEnt.at(l)==1)&&(TotPADNoiseEnt.at(l)==0))||((TotPADNoiseEnt.at(l)==1)&&(TotSTRIPNoiseEnt.at(l)==0)))
			

					if((TotHITNoiseEnt.at(l)==1)&&(hitongood))
					  {
					    firstHit.push_back((*ithit));		      
					    planeindexfirst.push_back(l);		
					  }
					//if(((TotSTRIPNoiseEnt.at(l)==2)&&(TotPADNoiseEnt.at(l)==0))||((TotPADNoiseEnt.at(l)==2)&&(TotSTRIPNoiseEnt.at(l)==0))||((TotPADNoiseEnt.at(l)==1)&&(TotSTRIPNoiseEnt.at(l)==1)))
					if(TotHITNoiseEnt.at(l)==2)  
					  {
					    secondHit.push_back((*ithit));
					    planeindexsecond.push_back(l);
					  }
					
					//std::cout<<" Debug-3"<<std::endl;
				      

				      }//end hit-collection Loop
				    
		      
				    int rCell=0;int phiCell=0;
				    T2Hit ExtrapolHit=GiveExtrapolatedHit(l,refhitv);

				    // (double Rcl, double Phicl, uint32_t thedet, int &rCell, int &phiCell)
				    ConvertCluRPhiInEffiGeoCell(ExtrapolHit.GetHitR(), ExtrapolHit.GetHitPhi(), ExtrapolHit.GetHitDetRawId(), rCell, phiCell);
				    int planenum=SelectedHalf*10+l;
				    //std::cout<<"Effi Geo Fill "<<planenum<<" "<<rCell<<" "<<phiCell<<" "<<ExtrapolHit.GetHitR()<<" "<<ExtrapolHit.GetHitPhi()<<std::endl;
				    
				    HGeometry_PadEfficiency_Den[planenum]->Fill(rCell,phiCell);
				    //std::cout<<"A "<<std::endl;
				    HGeometry_StripEfficiency_Den[planenum]->Fill(rCell,phiCell);    
				    // std::cout<<"b"<<std::endl;


				    Geometry_PadEfficiency_Den.at(planenum).at(rCell).at(phiCell)+=1.0;

				    //std::cout<<"c "<<std::endl;
				    Geometry_StripEfficiency_Den.at(planenum).at(rCell).at(phiCell)+=1.0;
				    //std::cout<<"d "<<std::endl;

				    if(padongood==true){				     
				      HGeometry_PadEfficiency_Num[planenum]->Fill(rCell,phiCell);				     
				      Geometry_PadEfficiency_Num.at(planenum).at(rCell).at(phiCell)+=1.0;
				      //  std::cout<<"e1 "<<std::endl;
				    }
				    
				    if(stripongood==true){
				      HGeometry_StripEfficiency_Num[planenum]->Fill(rCell,phiCell);
				      Geometry_StripEfficiency_Num.at(planenum).at(rCell).at(phiCell)+=1.0;
				      // std::cout<<"e2 "<<std::endl;
				    }
				    
				    //  std::cout<<"e3 "<<std::endl;
				    //std::cout<<" Debug-4"<<std::endl;
				    
				    bool stripvfatmapped=false;
				    bool padvfatmapped=false;
				    
				    // if((lhitindet==true)&&(RsectorDead==false))
				    //{

				    if(padon==true)
				      EFFpadVsIdNoCut->Fill(l,1.);
				    else
				      {
					EFFpadVsIdNoCut->Fill(l,0.);
				      }
				    
				    if(stripon==true)
				      EFFstripVsIdNoCut->Fill(l,1.);
				    else
				      EFFstripVsIdNoCut->Fill(l,0.);
				    
				    
				    T2Hit touseVfat;
				    if(padongood==true)
				      {
					touseVfat=PadHitgood;
					
					if(hitongood==true)
					  {
					    touseVfat=StripAndPadHitgood;
					  }
				      }

			  

			  if(padongood==true)
			    {
			      EFFpadVsId->Fill(l,1.);
			      T2GeometryUtil converter;
			     
			      //std::cout<<" Work with Real PadHitgood"<<touseVfat.GetHitR()<<"  "<<touseVfat.GetHitPhi()<<std::endl;
			      // std::vector<T2GeometryUtil::vfatid_channel> PadVFatVect = converter.PadVfatsIdsFromPadVect(PadHitgood);
			      T2GeometryUtil::vfatid_channel PadVFat=converter.PadVfatsIdFromRecoHit(touseVfat);
			      
			      if(PadVFat.vfatiid==-1)
				std::cout<<"Problem in Pad vfat association"<<PadHitgood.GetHitNumPad()<<" "<<touseVfat.GetHitNumPad()<<std::endl;

			      VFATEFF[l]->Fill(PadVFat.vfatiid,1);
			      VFATEFFNormalized[l]->Fill(PadVFat.vfatiid,1);
			      //std::cout<<" @@ real PadVFat.vfatiid="<<PadVFat.vfatiid<<std::endl;
			      /*
			      for(unsigned int h=0;h<PadVFatVect.size();h++)
			       {
				 VFATEFF[l]->Fill(PadVFatVect.at(h).vfatiid,1);
				 std::cout<<" @@"<<PadVFatVect.at(h).vfatiid<<std::endl;
			       }
			      */
			      VfatStatistics[l]->Fill(PadVFat.vfatiid);
			      EFFpadVsIdNormalized->Fill(l,1.);	
			      EFFpadVsIdNormalized_sectorX[Rsector]->Fill(l,1.);	
			      if(l==1)
				if(verbosity)
				  std::cout<<" ->->->->-> CLUSTER PAD EFF. OK"<<std::endl;

			    }
			  else
			    {
			     
			      T2Hit thehitPad=GiveExtrapolatedHit(l,refhitv);

			      if(HitInQuarterAnalyzed(thehitPad)==true)
				{

				  EFFpadVsId->Fill(l,0.);		
				  //std::cout<<" Work with Extrapol"<<std::endl;
				  T2GeometryUtil converter;
				  T2GeometryUtil::vfatid_channel PadVFat= converter.PadVfatsIdFromRecoHit(thehitPad); 
				  if(PadVFat.vfatiid==-1)
				    std::cout<<"Problem in Pad vfat association"<<thehitPad.GetHitNumPad()<<" "<<thehitPad.GetHitNumPad()<<std::endl;
 		     
				  VFATEFF[l]->Fill(PadVFat.vfatiid,0);
				  VfatStatistics[l]->Fill(PadVFat.vfatiid);
				  if(verbosity)
				    std::cout<<" @@ extrapol PadVFat.vfatiid="<<PadVFat.vfatiid<<std::endl;
				  int absVFsymb=(SelectedHalf*10+l)*100+PadVFat.vfatiid;
				  padvfatmapped=IsVfatMapped(absVFsymb);
				  
			      	

				  if(UseUncorrupetdEventMap==true)  //default option
				    {
				      T2VfatInformation::const_iterator itvf;
				      itvf= t2vfatinfoptr->find(absVFsymb);
				      if(itvf==t2vfatinfoptr->end())
					{
					  if(verbosity)
					    std::cout<<" Code Error itvf==t2vfatinfoptr->end()"<<std::endl;
				      
					}
				      else
					if(RsectorDead==false)
					if(itvf->second==0)
					  {
					    EFFpadVsIdNormalized->Fill(l,0.);  
					    EFFpadVsIdNormalized_sectorX[Rsector]->Fill(l,0.);
					    
					    VFATEFFNormalized[l]->Fill(PadVFat.vfatiid,0);
					    if(l==1)
					      if(verbosity)
						std::cout<<" ->->->->-> CLUSTER PAD INEFFICIENT. l="<<l<<std::endl;
					  }
					else
					  {
					    if(l==1)
					      if(verbosity)
						std::cout<<" ->->->->-> CLUSTER PAD EFF NOT CALCULATE BECAUSE OF CORRUPTION . "<<std::endl;
					  }
				      
				    }
				  else
				    if(RsectorDead==false)
				    if(padvfatmapped)//Analisi + grezza, normalizzi solo per i vfat morti.
				      {				  
					EFFpadVsIdNormalized->Fill(l,0.);
					EFFpadVsIdNormalized_sectorX[Rsector]->Fill(l,0.);
					VFATEFFNormalized[l]->Fill(PadVFat.vfatiid,0);
				      }
				}//if HitinquarterAnalyzed=true
			      
			    }
			  


			  if(stripongood==true)
			    {
			      touseVfat=StripHitgood;
			  
			      if(hitongood==true)
				{
				  touseVfat=StripAndPadHitgood;
				}
			    }

			  
			  
			  if((stripongood==true))
			    {
			      EFFstripVsId->Fill(l,1.);
			      T2GeometryUtil converter;			   
			      if(HitInQuarterAnalyzed(touseVfat)==true){
				T2GeometryUtil::vfatid_channel StripVFat=converter.StripVfatsIdFromRecoHit(touseVfat);
				if(StripVFat.vfatiid==-1)
				  std::cout<<"Problem HERE!"<<std::endl;
				
				VFATEFF[l]->Fill(StripVFat.vfatiid,1);
	
				 if(RsectorDead==false)
				   VFATEFFNormalized[l]->Fill(StripVFat.vfatiid,1);	
			
				 VfatStatistics[l]->Fill(StripVFat.vfatiid);
			
				if(StripVFat.vfatiid==15)
				  if(l==4)
				    if(verbosity)
				      std::cout<<" ~~~~~~~~~~~~~~~~~~~~ Vfat 15 on plane 4 ON"<<std::endl;
				
				if(verbosity)
				  if(l==4)
				    {
				      std::cout<<" ~~~~~~~~~~~~~~~~~~~~ Plane 4 STRIP ON"<<std::endl;
				      std::cout<<" ~~~~~~~~~~~~~ Strip-row: "<<(touseVfat.GetCluStripEntries()).at(0).rad_coord<<std::endl;
				      std::cout<<" ~~~~~~~~~~~~~ Strip-col: "<<(touseVfat.GetCluStripEntries()).at(0).ang_coord<<std::endl;
				      std::cout<<" ~~~~~~~~~~~~~VFAT STRIP id:"<<StripVFat.vfatiid<<std::endl;
				    }
				
				if(RsectorDead==false){
				  EFFstripVsIdNormalized->Fill(l,1.);
				  EFFstripVsIdNormalized_sectorX[Rsector]->Fill(l,1.);
				}
				
				if(l==1)
				  if(verbosity)
				    std::cout<<" ->->->->-> CLUSTER STRIP EFF. OK"<<std::endl;
			      }
			    }
			  else
			    {
				  EFFstripVsId->Fill(l,0.);
				  T2Hit thehitStrip=GiveExtrapolatedHit(l,refhitv);	

				  if(HitInQuarterAnalyzed(thehitStrip)==true){
				    
				    //Use only geometrical reco info
				    T2GeometryUtil converter;
				    T2GeometryUtil::vfatid_channel StripVFat= converter.StripVfatsIdFromRecoHit(thehitStrip);
				    if(StripVFat.vfatiid==-1)
				      std::cout<<"Problem HERE b!"<<std::endl;
				    
				    // vfatstripactive=GetVfatForHit(l,refhitv);	
				    if(StripVFat.vfatiid==15)
				      if(l==4)
					if(verbosity)
					  std::cout<<" ~~~~~~~~~~~~~~~~~~~~ Vfat 15 on plane 4 OFF. R:"<<thehitStrip.GetHitR()<<" Phi: "<<thehitStrip.GetHitPhi()<<std::endl;
				    
				    if(l==4)
				      if(verbosity)
					std::cout<<" ~~~~~~~~~~~~~~~~~~~~ Plane 4 STRIP OFF"<<std::endl;
				    
				    VFATEFF[l]->Fill(StripVFat.vfatiid,0);
				    VfatStatistics[l]->Fill(StripVFat.vfatiid);
				    if(verbosity)
				      std::cout<<" @@ extrapol StripVFat.vfatiid="<<StripVFat.vfatiid<<" plane "<<l<<std::endl;
				    int absVFsymb=(SelectedHalf*10+l)*100+StripVFat.vfatiid;
				    stripvfatmapped=IsVfatMapped(absVFsymb);//Found in the alive list
				    
				    
				    if(UseUncorrupetdEventMap==true)
				      {
					T2VfatInformation::const_iterator itvf;
					itvf= t2vfatinfoptr->find(absVFsymb);
					if(itvf==t2vfatinfoptr->end())
					  {
					    if(verbosity)
					      std::cout<<" Code Error itvf==t2vfatinfoptr->end()"<<std::endl;
					  }
					else
					   if(RsectorDead==false)
					  if(itvf->second==0)
					    {
					      EFFstripVsIdNormalized->Fill(l,0.);
					      EFFstripVsIdNormalized_sectorX[Rsector]->Fill(l,0.);
					      VFATEFFNormalized[l]->Fill(StripVFat.vfatiid,0);
					      if(l==1)
						if(verbosity)
						  std::cout<<" ->->->->-> CLUSTER STRIP INEFFICIENT. l="<<l<<std::endl;
					    }
					  else
					    {
					      if(l==1)
						if(verbosity)
						  std::cout<<" ->->->->-> CLUSTER STRIP EFF NOT CALCULATE BECAUSE OF CORRUPTION . "<<std::endl;
					    }
				      }
				    else
				       if(RsectorDead==false)
					 if(stripvfatmapped)
					   {
					     EFFstripVsIdNormalized->Fill(l,0.);
					     EFFstripVsIdNormalized_sectorX[Rsector]->Fill(l,0.);
					     VFATEFFNormalized[l]->Fill(StripVFat.vfatiid,0);
					   }
				    //}
				  }
				  
			    }
			  
			  
			  //std::cout<<" Debug-5"<<std::endl;

			  if((padongood==true)||(stripongood==true))
			    EFFClusterVsId->Fill(l,1.);
			  else //both are false
			    {			     
			      if((NoisyplaneForStrip==false)&&(NoisyplaneForPad==false))
				EFFClusterVsId->Fill(l,0.);			 
			    }


			  if(hitongood==true)
			    {
			      
			      EFFHitVsId->Fill(l,1.); 

			      if(RsectorDead==false)
				EFFHitVsIdNormalized->Fill(l,1.);

			      NoiseHitDet[l]->Fill(TotHITNoiseEnt.at(l)-1);
			      NoiseStripDet[l]->Fill(TotSTRIPNoiseEnt.at(l)-1);
			      NoisePadDet[l]->Fill(TotPADNoiseEnt.at(l)-1);


			      //std::cout<<TotHITNoiseEnt.at(l)<<std::endl;		
			      if (TotHITNoiseEnt.at(l)>1)
				Noise01HitDet[l]->Fill(1.);
			      else
				Noise01HitDet[l]->Fill(0.);
			      			    
			      if (TotPADNoiseEnt.at(l)>1)
				Noise01PadDet[l]->Fill(1.);
			      else
				Noise01PadDet[l]->Fill(0.);

			      if (TotSTRIPNoiseEnt.at(l)>1)
				Noise01StripDet[l]->Fill(1.);
			      else
				Noise01StripDet[l]->Fill(0.);
  
			    }
			  else
			    {
			      EFFHitVsId->Fill(l,0.);
			      T2Hit thehitExtrap=GiveExtrapolatedHit(l,refhitv);
			       if(HitInQuarterAnalyzed(thehitExtrap)==true){
			      
			      if(verbosity)
				std::cout<<" Plane "<<l<<" with "<<" NO class1 entries"<<std::endl;
			      T2GeometryUtil converter;
			      T2GeometryUtil::vfatid_channel StripVFatC1= converter.StripVfatsIdFromRecoHit(thehitExtrap);
			       if(StripVFatC1.vfatiid==-1)
				std::cout<<"Problem HERE!c"<<std::endl;
			      
			      T2GeometryUtil::vfatid_channel PadVFatC1= converter.PadVfatsIdFromRecoHit(thehitExtrap);
			      
			      if(verbosity)
				{
				  std::cout<<" @@ extrapol PadVFatC1.vfatiid="<<PadVFatC1.vfatiid<<std::endl;
				  std::cout<<" @@ extrapol StripVFatC1.vfatiid="<<StripVFatC1.vfatiid<<std::endl;

				  if(padongood)
				    std::cout<<" padongood: true ";
				  else
				    std::cout<<" padongood: false ";

				  if(stripongood)
				    std::cout<<" stripongood: true ";
				  else
				    std::cout<<" stripongood: false ";

				  if(hitongood)
				    std::cout<<" hitongood: true ";
				  else
				    std::cout<<" hitongood: false ";

				  std::cout<<"  "<<std::endl;

				}
			      //put EffHitNormalized to 0 only when striphit and padhit mapped 
			      int absVFsymbstr=(SelectedHalf*10+l)*100+StripVFatC1.vfatiid;
			      int absVFsymbpad=(SelectedHalf*10+l)*100+PadVFatC1.vfatiid;
			      
			      
			      if(UseUncorrupetdEventMap==true)
				{
				  T2VfatInformation::const_iterator itvfp;
				  T2VfatInformation::const_iterator itvfs;
				  itvfs= t2vfatinfoptr->find(absVFsymbstr);
				  itvfp= t2vfatinfoptr->find(absVFsymbpad);
				  if((itvfs==t2vfatinfoptr->end())||(itvfp==t2vfatinfoptr->end()))
				    {
				      if(verbosity)
					std::cout<<" Code Error itvf==t2vfatinfoptr->end()"<<std::endl;
				    }
				  else
				    if((itvfs->second==0)&&(itvfp->second==0))
				       if(RsectorDead==false)
					 EFFHitVsIdNormalized->Fill(l,0.);
				}
			      else // analisi + rozza (iniziale, tiene solo conto del mapping)
				{
				  if(RsectorDead==false){
				    if((IsVfatMapped(absVFsymbpad))&&(IsVfatMapped(absVFsymbstr)))				      
				      EFFHitVsIdNormalized->Fill(l,0.);				  
				  else
				    {
				      if((IsVfatMapped(absVFsymbstr))==false)
					{
					  if(verbosity)
					    std::cout<<" stripvfat not mapped"<<std::endl;
					}
				      if((IsVfatMapped(absVFsymbpad))==false)
					{
					  if(verbosity)
					    std::cout<<" padvfat not mapped"<<std::endl;				  
					}
				    }
				  }
				}
			      
			    }

			    }

			    

			  if(hiton==true)
			    EFFHitVsIdNoCut->Fill(l,1.);
			  else
			    EFFHitVsIdNoCut->Fill(l,0.);			

			  
		    }//if HitInDet and RSectorDead==false



		      //std::cout<<" Debug-6"<<std::endl;



		       //----------------------
		      // Histrograms Noise DR
		      //----------------------
		      //Average_Cl1HitNoiseNumb_vsDet
				/*
		      int NumberOfCl1HitNoiseWhenTracks=0; 
		      for(T2HitCollection::const_iterator ithit = t2hitcoll->begin(); ithit != t2hitcoll->end(); ithit++)
			{
			
			 	 
			  unsigned int  symb=(*ithit).GetHitPlane()*2+(*ithit).GetHitPlaneSide();
			  
			  bool IsATrackingHit=false;
			  
			  if(HitIsInTrackColl((*ithit),t2trackVectorALL))
			    IsATrackingHit=true;
			  			

			  if((symb==l)&&(IsATrackingHit==false))
			    {
			      
			      if(hitongood)
				{
				  NumberOfCl1HitNoiseWhenTracks++;
				  if(((*ithit).GetHitNumStrip()>0)&&((*ithit).GetHitNumPad()>0))
				    {
				      Class1NoiseDRDPHISep[l]->Fill(StripAndPadHitgood.GetHitR()-(*ithit).GetHitR(),StripAndPadHitgood.GetHitPhi()-(*ithit).GetHitPhi());
				      TrackingHitVsClass1NoiseDRSep[l]->Fill(StripAndPadHitgood.GetHitR(),StripAndPadHitgood.GetHitR()-(*ithit).GetHitR());
				    }
				  else
				    {
				      if((*ithit).GetHitNumStrip()>0)
					{
					  TrackingHitRVsStripCluNoiseDRSep[l]->Fill(StripAndPadHitgood.GetHitR(),StripAndPadHitgood.GetHitR()-(*ithit).GetHitR());
					}
				      else
					{
					  TrackingHitRVsPadCluNoiseDRSep[l]->Fill(StripAndPadHitgood.GetHitR(),StripAndPadHitgood.GetHitR()-(*ithit).GetHitR());
					}
				    }
				  
				}
			    
			      
			      
			    }
			
			}
		     
		      Average_Cl1HitNoiseNumb_vsDet->Fill(l,NumberOfCl1HitNoiseWhenTracks);
		     
		      for(unsigned int h=0;h<planeindexsecond.size();h++)
			{
			  double dr=0.;
			  for(unsigned int u=0;u<planeindexfirst.size();u++)
			    {
			      if(planeindexfirst.at(u)==planeindexsecond.at(h))
				dr=firstHit.at(u).GetHitR();
			    }
			  dr=dr-secondHit.at(h).GetHitR();
			  DRHit12->Fill(planeindexsecond.at(h),dr);
			  if(planeindexsecond.at(h)<10)
			    drHit12[planeindexsecond.at(h)]->Fill(dr);
			  else
			    if(verbosity)
			      std::cout<<" Warning, drhistogram"<<std::endl;
			}
		      //---------------------
		      */

		    
		    
		    }//if(countref>=Effgoodhitnumber)


			   
			  


		  //  alreadyprint0=false;
		}   // for-loop on l-detector

	      
			





	      // Other global histograms
	      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	      std::vector<unsigned int> numhitinplaneGood;
	      /*
	      for(unsigned int k=0;k<10;k++)
		{
		  numhitinplaneGood.push_back(0);
		}
		
	      for(T2HitCollection::const_iterator ithit = t2hitcoll->begin(); ithit != t2hitcoll->end(); ithit++){
		unsigned int symb=(*ithit).GetHitPlane()*2+(*ithit).GetHitPlaneSide();
		
		numhitinplaneGood.at(symb)++;
	      }
		
	      for(unsigned int k=0;k<10;k++)
		{
		  NumAllHitvsIdGoodevt->Fill(k,numhitinplaneGood.at(k));  
		}
	      */
	      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



	      }//if(TrackInMatrixRow.goodhitnumber>=Effgoodhitnumber)
	  
	    
	  }//if inserted
       } //if Track in right half	   	  
    }//end Track::event loop
      /*if(verbosity)*/ // std::cout<<" Event complete"<<std::endl;
      

     
      if(EventUtilizedForEffi)
	for(T2PadClusterCollection::const_iterator itpad = t2padclcoll->begin(); itpad != t2padclcoll->end(); itpad++){
	  
	  myT2Det = itpad->first;
	  myhalftele= myT2Det.halfTelescope();
	  myplaneside= myT2Det.planeSide();
	  myplane= myT2Det.plane();
	  unsigned int myarm=myT2Det.arm();
	  
	  double rr=0.; double ff=0.; double xx=0.; double yy=0.;
	  
	  unsigned int  symb= myplane*2 + myplaneside + myhalftele*10 + myarm*20;
	
	  if((symb/10)==SelectedHalf){

	    vector<T2Cluster> padClv = itpad->second;
	    for(unsigned int k=0;k<padClv.size();k++){
	      
	      rr=padClv[k].GetClusterR();
	      ff=padClv[k].GetClusterPhi()*3.14159/180.;
	      
	      xx=rr*cos(ff);
	      yy=rr*sin(ff);
	      //std::cout<<" Woila'"<<std::endl;
	      PadClusterXYForEffiEvents[(symb%10)]->Fill(xx,yy);
	      
	      if(symb<10)
		if(padClv[k].GetClusterPhi()>96)
		  if(padClv[k].GetClusterPhi()<264)
		    std::cout<<"Pippi-Bug is living in Pad Clu Phi="<<padClv[k].GetClusterPhi()<<" plane:"<<symb<<" entries:"<<padClv[k].GetNoOfEntries()<<std::endl;
	    	     
	    }
	  }
	  
	}			  

    }//if numevent<maxevents
  
  //std::cout<<" Debug-7"<<std::endl;





//std::cout<<" Event Track Size= "<<trackCollection->size()<<std::endl;

  //------------------------------
  //Begin Noise Det selected Histograms
  //------------------------------
  NoiseEventsStat->Fill(1);
 
  std::vector<unsigned int> Allcl1hitPlanecounter;
  std::vector<unsigned int> AllStripCluPlanecounter;
  std::vector<unsigned int> AllPadCluPlanecounter;
  for(unsigned int k=0;k<10;k++)
    {
      Allcl1hitPlanecounter.push_back(0);
      AllStripCluPlanecounter.push_back(0);
      AllPadCluPlanecounter.push_back(0);
    }

  for(T2HitCollection::const_iterator ithit = t2hitcoll->begin(); ithit != t2hitcoll->end(); ithit++)
    {
      
      uint32_t rawiddet=(*ithit).GetHitDetRawId();
      unsigned int  symb=RawtoSymb(rawiddet);  
      std::vector<cluster_entry> entriespadcl= (*ithit).ClusterPad_entries;
     
      unsigned int theplane=symb%10;
      if(entriespadcl.size()>0)
	PadCluSizeVsPlaneAll2->Fill(theplane,entriespadcl.size());
      if((symb/10)==SelectedHalf)
	{
	  if(((*ithit).GetHitNumStrip()>0)&&((*ithit).GetHitNumPad()>0))
	    {
	      Allcl1hitPlanecounter.at(theplane)++;
	      AllStripCluPlanecounter.at(theplane)++;
	      AllPadCluPlanecounter.at(theplane)++;
	      PadCluSizeVsPlaneAll->Fill(theplane,(*ithit).GetHitNumPad());
	      StipCluSizeVsPlaneAll->Fill(theplane,(*ithit).GetHitNumStrip());	
	    }
	  else 
	    {
	      if((*ithit).GetHitNumStrip()>0)
		{
		  AllStripCluPlanecounter.at(theplane)++;
		  StipCluSizeVsPlaneAll->Fill(theplane,(*ithit).GetHitNumStrip());	
		}
	      if((*ithit).GetHitNumPad()>0)
		{
		  AllPadCluPlanecounter.at(theplane)++;
		  PadCluSizeVsPlaneAll->Fill(theplane,(*ithit).GetHitNumPad());  
		}
	    }
	}



      if(((*ithit).GetHitNumStrip()>0)&&((*ithit).GetHitNumPad()>0))
	Class1HitPadStripCLSCorrel->Fill((*ithit).GetHitNumStrip(),(*ithit).GetHitNumPad());

      if((symb==DetForNoiseStudies))
	if(((*ithit).GetHitPhi()>PhiMinForNoiseStudies)&&((*ithit).GetHitPhi()<PhiMaxForNoiseStudies))
	  {
	    //if((*ithit).GetHitNumPad()>0)
	    if(trackCollection->size()==0)
	      if(((*ithit).GetHitNumStrip()>0)&&((*ithit).GetHitNumPad()>0))
		 if(verbosity) std::cout<<" R: "<<(*ithit).GetHitR()<<" Phi: "<<(*ithit).GetHitPhi()<<" Z: "<<(*ithit).GetHitZ()<<" Class1"<<std::endl;
	      else
		 if(verbosity) std::cout<<" R: "<<(*ithit).GetHitR()<<" Phi: "<<(*ithit).GetHitPhi()<<" Z: "<<(*ithit).GetHitZ()<<" Class0"<<std::endl;

	    if(((*ithit).GetHitNumStrip()>0)&&((*ithit).GetHitNumPad()==0))
	      {
		RNoiseStrip->Fill((*ithit).GetHitR()); 
		CLSRNoiseStrip->Fill((*ithit).GetHitNumStrip());   
	      }
	    
	    if(((*ithit).GetHitNumPad()>0)&&((*ithit).GetHitNumStrip()==0))
	      {
		PhiNoisePad->Fill((*ithit).GetHitPhi()); 
		CLSPhiNoisePad->Fill((*ithit).GetHitNumPad()); 
	      }
	    
	    if(((*ithit).GetHitNumPad()>0)&&((*ithit).GetHitNumStrip()>0))
	      {
		PhiNoiseCl1Hit->Fill((*ithit).GetHitPhi()); 
		RNoiseCl1Hit->Fill((*ithit).GetHitR()); 
		CLSPadNoiseCl1Hit->Fill((*ithit).GetHitNumPad());
		CLSStripNoiseCl1Hit->Fill((*ithit).GetHitNumStrip());
	      }
	    
	    if(((*ithit).GetHitNumPad()>0)&&((*ithit).GetHitNumStrip()>0))
	      {
		Cl1Hit_ClsStripVsR->Fill((*ithit).GetHitR(),(*ithit).GetHitNumStrip());//TProfile
		C11Hit_ClsPadVsR->Fill((*ithit).GetHitR(),(*ithit).GetHitNumPad());//TProfile
		C11Hit_ClsPadVsClsStrip->Fill((*ithit).GetHitNumPad(),(*ithit).GetHitNumStrip());//Th2F-3D		
		C1Hit_DphiVsPadCLS->Fill((*ithit).GetHitDPhi(),(*ithit).GetHitNumPad());//Th2F-3D
		

		
		std::vector<cluster_entry> entriespadcl= (*ithit).ClusterPad_entries;
		//	std::cout<<" ||||||||||----- entriespadcl.size()= "<<ithit->ClusterPad_entries.size()<<"-------||||||||"<<std::endl;
		//	std::cout<<" ||||||||||----- entriesstripcl.size()= "<<ithit->ClusterStrip_entries.size()<<"-------||||||||"<<std::endl;
		//	std::cout<<" |||||||||| entriespadcl.size()= "<<entriespadcl.size()<<"||||||||"<<std::endl;
		int rmin=4000;
		int rmax=-1;
		int cmin=4000;
		int cmax=-1;
		for(unsigned int y=0;y<entriespadcl.size();y++)
		  { 
		    //std::cout<<" r-c"<<entriespadcl.at(y).rad_coord<<"-"<<entriespadcl.at(y).ang_coord<<std::endl;
		    if(entriespadcl.at(y).rad_coord<=rmin)
		      {
			rmin=entriespadcl.at(y).rad_coord;
		      }
		    if(entriespadcl.at(y).rad_coord>=rmax)
		      {
			rmax=entriespadcl.at(y).rad_coord;
		      }
		    if(entriespadcl.at(y).ang_coord<=cmin)
		      {
			cmin=entriespadcl.at(y).ang_coord;
		      }

		    if(entriespadcl.at(y).ang_coord>=cmax)
		      {
			cmax=entriespadcl.at(y).ang_coord;
		      }
		    
		    

		  }
	
		if(entriespadcl.size()>0)
		  {
		    // std::cout<<" cmax end cmin: "<<cmax<<" "<<cmin<<std::endl;
		    // std::cout<<" (cmax-cmin)-(rmax-rmin)"<<(cmax-cmin)<<"-"<<(rmax-rmin)<<std::endl;
		    C1Hit_PadClsVsPadNumCol->Fill((*ithit).GetHitNumPad(),(cmax-cmin));
		    C1Hit_PadClsVsPadNumRow->Fill((*ithit).GetHitNumPad(),(rmax-rmin));
		    C1Hit_Pad_PadNumRowVsPadNumCol->Fill((rmax-rmin),(cmax-cmin));
		    //C11Hit_ClsPadVsClsPadCol->Fill((*ithit).GetHitNumPad(),(*ithit).GetAngularSpread());//Th2F-3D
		  }


		C1Hit_Dphi->Fill((*ithit).GetHitDPhi());
	      }
	  }
          
    }	//End T2Hit loop

  for(unsigned int k=0;k<10;k++)
    {
      NumPadCluVsPlane->Fill(k,AllPadCluPlanecounter.at(k));
      NumStripCluVsPlane->Fill(k,AllStripCluPlanecounter.at(k));
      NumCl1HitVsPlane->Fill(k,Allcl1hitPlanecounter.at(k));
    }
  //------------------------------
  //End Noise Det selected Histograms
  //------------------------------













 
  if(numevent==MaxEvents)  // fai l'analisi
    {
     
     
      for(unsigned int u=0;u<10;u++)
	{
	 if(verbosity) std::cout<<" Det_"<<u<<"||";
	}
      

      std::cout<<" MaxEvent Reached: Eventi della matrice "<< MatriceTracce.size()<<std::endl;
      // if(verbosity) std::cout << "Alignment will work with "<< MaxEvents <<" events and " <<roadXfit.size()<<"-"<<roadYfit.size()<<" X-Y Roads"<<" having hits in "<<alldetid.size() <<" different detectors."<<std::endl;

      // std::cout<<DoALign<<std::endl;

      /*

      //CONSTRUCTION OF HISTOGRAMS

      for(unsigned int k=0; k<roadXfit.size();k++)
	{
	  if(roadXfit.size()==origroadXfit.size())
	    {
	      for(unsigned int j=0; j<roadXfit.at(k).size();j++)
		{
		  
		  if(roadXfit.at(k).at(j).GetHitDetRawId()==origroadXfit.at(k).at(j).GetHitDetRawId())
		    {
		      unsigned int symbolic=RawtoSymb(roadXfit.at(k).at(j).GetHitDetRawId());		      
		      double diffX=roadXfit.at(k).at(j).GetHitX()-origroadXfit.at(k).at(j).GetHitX();		      
		      DetIdvsXShOnlyXfit->Fill(symbolic,diffX);
		      
		    }
		}
	    }
	}


      for(unsigned int k=0; k<roadYfit.size();k++)
	{
	  if(roadYfit.size()==origroadYfit.size())
	    {
	      for(unsigned int j=0; j<roadYfit.at(k).size();j++)
		{
		  
		  if(roadYfit.at(k).at(j).GetHitDetRawId()==origroadYfit.at(k).at(j).GetHitDetRawId())
		    {
		     
		      unsigned int symbolic=RawtoSymb(roadYfit.at(k).at(j).GetHitDetRawId());
		      double diffY=roadYfit.at(k).at(j).GetHitY()-origroadYfit.at(k).at(j).GetHitY();		      
		      DetIdvsYShOnlyYfit->Fill(symbolic,diffY);
		    }
		}
	    }
	}

      
      //START Quality Plot after alignment
      for(unsigned int k=0; k<roadYfit.size();k++)
	{
	  ProbYHistoAl->Fill(chi2Y(roadYfit.at(k),UseJointProb).at(0));
	}
       for(unsigned int k=0; k<roadXfit.size();k++)
	{
	  ProbXHistoAl->Fill(chi2X(roadXfit.at(k),UseJointProb).at(0));
	}

       //Fill R Resolution Plot (2)
       
       if(roadYfit.size()==roadXfit.size())
       for(unsigned int k=0; k<roadYfit.size();k++)
	{

	  std::vector<T2Hit> refhitv;
	  for(unsigned int j=0;j<roadYfit.at(k).size();j++)
	    {
	      T2Hit thehit = roadYfit.at(k).at(j);
	      thehit.SetHitX(roadXfit.at(k).at(j).GetHitX());
	      refhitv.push_back(thehit);
	    }

	   std::vector<double> drphi4Alignment;

	   TMatrixD covmatr(4,4);
	   covmatr.Zero();
	   double chi2corre;	      
	   std::vector<double> correFit=MyLinearfitCorr(refhitv,covmatr,chi2corre);
	   //std::vector<double> correFit=OLDMyLinearfit(std::vector<T2Hit> hitvec)
	   //
	   std::vector<double> support2;
	   
	   for(unsigned int c=0;c<2;c++)
	     {
	       support2.push_back(correFit[c]);
	     }
	   

	   for(unsigned int c=0;c<2;c++)
	     {
	       support2.push_back(correFit[2+c]);
	     }

	   for (unsigned int c=0;c<2;c++)
	     {
	       support2.push_back(covmatr(c,c));
	     }
	   

	   for (unsigned int c=0;c<2;c++)
	     {
	       support2.push_back(covmatr(2+c,2+c));
	     }
	   support2.push_back(covmatr(0,1));
	   support2.push_back(covmatr(2,3));
	      
	   std::vector<double> paramvectForAlign;
	   paramvectForAlign=support2;
	  
	   if(useRZforResol==2)
	     paramvectForAlign=MyLinearfit(refhitv,UseJointProb); 
	    //paramvectForAlign=MyLinearfit(refhitv,UseJointProb); 
	    //  if(chi2X(refhitv)>0.05)
	    //  if(chi2Y(refhitv)>0.05)
	   if(refhitv.size()>=(HitNumb4Align-1))
	     for (unsigned int u=0;u<refhitv.size();u++)
	       {
		 drphi4Alignment=ResiduiRPhi(paramvectForAlign,refhitv.at(u).GetHitX(), refhitv.at(u).GetHitY(),refhitv.at(u).GetHitZ());
		 unsigned int  symb=refhitv.at(u).GetHitPlane()*2+refhitv.at(u).GetHitPlaneSide();
		 if(refhitv.at(u).GetHitNumStrip()<=4)
		   if(refhitv.at(u).GetHitNumPad()<=3)
		     {
		       if(refhitv.at(u).GetHitR()<AlignmentHitRMax)
			 {
			   DXAlignDet[symb]->Fill(drphi4Alignment.at(0));
			   DYAlignDet[symb]->Fill(drphi4Alignment.at(1));
				   
			   //DRResol[symb]->Fill(sqrt(drphi4Alignment.at(1)*drphi4Alignment.at(1)+drphi4Alignment.at(0)*drphi4Alignment.at(0)));
			   double expX=(paramvectForAlign[0]*refhitv.at(u).GetHitZ()+paramvectForAlign[1]);
			   
			   //vedi pag 130  lib blu;

			  
			  // vect[0]=a_rz;
			  // vect[1]=b_rz;
			  // vect[2]=a_phiz;
			  // vect[3]=b_phiz;

			  // vect[4]=e_a_rz;
			  // vect[5]=e_b_rz;
			  // vect[6]=e_a_phiz;
			  // vect[7]=e_b_phiz;

			  // double correl=(-1.0)*Sz_r*(1.0/(S0_r*Szz_r-(Sz_r*Sz_r)));
			  // vect[8]=correl;
			  // correl=(-1.0)*Sz_phi*(1.0/(S0_phi*Szz_phi-(Sz_phi*Sz_phi)));
			  // vect[9]=correl;
			   

			   double phirad=refhitv.at(u).GetHitPhi()*3.14159/180.0;
			   double sigmax=cos(phirad)*cos(phirad)*0.12*0.12+sin(phirad)*sin(phirad)*0.015*0.015*refhitv.at(u).GetHitR()*refhitv.at(u).GetHitR();			     		
			   double sigmay=sin(phirad)*sin(phirad)*0.12*0.12+cos(phirad)*cos(phirad)*0.015*0.015*refhitv.at(u).GetHitR()*refhitv.at(u).GetHitR();
			 
			   HitXErr->Fill(sqrt(sigmax));
			   HitYErr->Fill(sqrt(sigmay)); 


			   double errorx=paramvectForAlign[5]*paramvectForAlign[5];
			   errorx=errorx+2*paramvectForAlign[8]*refhitv.at(u).GetHitZ()+paramvectForAlign[4]*paramvectForAlign[4]*refhitv.at(u).GetHitZ()*refhitv.at(u).GetHitZ();
			   errorx=sqrt(errorx);

			   //  std::cout<<" X corr="<<paramvectForAlign[8]<<" Z="<<refhitv.at(u).GetHitZ()<<" E_ax="<<paramvectForAlign[4]<<" E_bx="<<paramvectForAlign[5]<<std::endl;
			   double errory=paramvectForAlign[7]*paramvectForAlign[7];
			   errory=errory+2*paramvectForAlign[9]*refhitv.at(u).GetHitZ()+paramvectForAlign[6]*paramvectForAlign[6]*refhitv.at(u).GetHitZ()*refhitv.at(u).GetHitZ();
			   errory=sqrt(errory);
			   //std::cout<<" EexpX= "<<errorx<<" EexpY= "<<errory<<std::endl;
			   
			   /////////////////////////////////
			   double expY=(paramvectForAlign[2]*refhitv.at(u).GetHitZ()+paramvectForAlign[3]);  
			   double expPhi=atan(fabs(expY)/fabs(expX));
			   expPhi=expPhi*180.0/3.14159;

			   //if(expY<0)&&)
			   //  expPhi=360.0-expPhi;
				
			   //Added on 15-10-09
			   if((expY<0)&&(expX>0))
			     expPhi=360.0-expPhi;

			   if((expY>0)&&(expX<0))
			     expPhi=expPhi+90.;

			   if((expY<0)&&(expX<0))
			     expPhi=expPhi+180.;




			   double expR= sqrt(expX*expX+expY*expY);
				
			   double measuredphi=refhitv.at(u).GetHitPhi();
			   double mesuredR=refhitv.at(u).GetHitR();  

			   double measuredphi2=atan(fabs(refhitv.at(u).GetHitY())/fabs(refhitv.at(u).GetHitX()));
			   measuredphi2=measuredphi2*180.0/3.14159;
			   //if(refhitv.at(u).GetHitY()<0)
			     //  measuredphi2=360.0-measuredphi2;
			   //Added on 15-10-09
			   if((refhitv.at(u).GetHitY()<0)&&(refhitv.at(u).GetHitX()>0))
			     measuredphi2=360.0-measuredphi2;

			   if((refhitv.at(u).GetHitY()>0)&&(refhitv.at(u).GetHitX()<0))
			     measuredphi2=measuredphi2+90.;

			   if((refhitv.at(u).GetHitY()<0)&&(refhitv.at(u).GetHitX()<0))
			     measuredphi2=measuredphi2+180.;


			   double mesuredR2=sqrt(refhitv.at(u).GetHitX()*refhitv.at(u).GetHitX()+refhitv.at(u).GetHitY()*refhitv.at(u).GetHitY());
			   HitX->Fill(refhitv.at(u).GetHitX());
			   HitY->Fill(refhitv.at(u).GetHitY());  
			   
			  
			   CumulativeRResol->Fill(expR-mesuredR2);
			   CumulativePhiResol->Fill(expPhi-measuredphi2);
			   
			   CumulativeExpXUncert->Fill(errorx);
			   CumulativeExpYUncert->Fill(errory);

			   DRResol[symb]->Fill(expR-mesuredR2);  
			   DPhiResol[symb]->Fill(expPhi-measuredphi2);
				
			 }
		     }
	       }

	}
*/

       // for (unsigned int y=0;y<1000;y++)  
       //	 RandomGauss->Fill(0.2*static_cast<int>(ceil(RandGaussQ::shoot( 0,400*12 ))));

    }// if numevent==maxevent

  double chiRProb=0.;
  double chiPhiProb=0.;
  double chiPhiRed=0.;

  for(TrkCit=trackCollection->begin(); TrkCit!=trackCollection->end(); TrkCit++){
    //std::cout<<" Traccia "<<trkcounter<<"  |||| "<<std::endl;
    
/*    
     if(Testedcamera.size()<5)
      {
	Testedcamera.clear();
	Testedcamera.push_back(0);
	Testedcamera.push_back(1);
	Testedcamera.push_back(2);
	Testedcamera.push_back(3);
	Testedcamera.push_back(4);
	//std::cout<<" Warning: Testedcamera is expected to have 5 elements !!"<<std::endl;
      }
*/
	    
    chiPhiRed=(*TrkCit).ChiSquaredPhi()/((*TrkCit).GetHitEntries()-1);
   
    

    chiRProb=TMath::Prob((*TrkCit).ChiSquaredR(),((*TrkCit).GetHitEntries()-2));
    chiPhiProb=TMath::Prob((*TrkCit).ChiSquaredPhi(),((*TrkCit).GetHitEntries()-1));  
    Chi2phired->Fill(chiPhiRed);
    Chi2RProb->Fill(chiRProb);
    Chi2PhiProb->Fill(chiPhiProb);
    
  }


  // **************************************************** VFAT STUDIES ************************************************ //
  // **************************************************** VFAT STUDIES ************************************************ //
  if(VFATMonitoring){

  unsigned int pl=0;
  unsigned int vfpos=0;
//  unsigned int absid=0;
  unsigned int vfatMult=0;
  unsigned int NumTrk=0;

  unsigned int mult15_17=0;
  unsigned int mult15_19=0;
  unsigned int mult11_19=0;
  unsigned int numvfatComONH1=0;
  std::vector<unsigned int> VF15DispmultH1(5);
  std::vector<unsigned int> VF16DispmultH1(5);
  unsigned int VFMultPl10_VF15=0;
  unsigned int VFMultPl10_VF16=0;

  Double_t *H1NoiseNtupla = new double[NoisyStripVfatH1_PlIId.size()];
  Short_t *H1NoiseNtuplaS= new short[NoisyStripVfatH1_PlIId.size()];

  for(unsigned int kk=0; kk<NoisyStripVfatH1_PlIId.size();kk++){
    H1NoiseNtupla[kk]=0;
    H1NoiseNtuplaS[kk]=0;
  }

  for(itv=  DigiVFatptr->begin();itv!=DigiVFatptr->end(); ++itv)
    {
      T2DetId mydet=(*itv).first;  
      pl= RawtoSymb(mydet.calculateRawId(mydet.arm(),mydet.halfTelescope(),mydet.plane(),mydet.planeSide()));
      if(pl<40)
	{
	  
	  if(pl>=30)
	    NumTrk=count_t2trackVectorALL_H2;//H3
	  else{
	    if(pl>=20)
	      NumTrk=count_t2trackVectorALL_H3;//H2
	    else{
	      if(pl>=10)
		NumTrk=count_t2trackVectorALL_H0;//H1
	      else
		NumTrk=count_t2trackVectorALL_H1;//H0
	    }
	  }
	    

	  for(std::vector<T2DigiVfat>::const_iterator itvfat =(*itv).second.first; itvfat !=(*itv).second.second; itvfat++)
	    {
	      vfpos=itvfat->GetVfatPos();
	      //std::cout<<"Found event record with vfat position "<<vfpos<<" in plane "<<pl<<std::endl;
	      
	      //Id for  data VFAT

	      std::map<const unsigned int,unsigned int>::const_iterator itch;
	      vfatMult=0;
	      for(itch=(*itvfat).ChActMap.begin();itch != (*itvfat).ChActMap.end();itch++){
		
		if(((*itch).second)==1)//Active Channel
		  {
		    vFatCumulative[pl][vfpos]->Fill((*itch).first); //Fill The Bin corresponding to the channel
		    //std::cout<<"vFatPlots["<<absid<<"].ActiveChCumulative->Fill("<<(*itch).first<<");"<<std::endl;
		    vfatMult++;
		  }		
	      }

	      vFatMultiplicty[pl][vfpos]->Fill(vfatMult);
	      vFatMultiplictyVsTrk[pl][vfpos]->Fill(NumTrk,vfatMult);

	      if((pl==10)&&(vfpos==15))
		VFMultPl10_VF15=vfatMult;

	      if((pl==10)&&(vfpos==16))
		VFMultPl10_VF16=vfatMult;


	      

	      if(vfatMult>100){		
		TrkMultiplictyWhenEveryThingIsOn[pl][vfpos]->Fill(NumTrk);
		if(((pl/10)==1)&&((vfpos>14)||(vfpos<2))){
		  numvfatComONH1++;
		  std::pair<int,int> NoisyVF(pl,vfpos);
		  bool SaveForCorrel= find(NoisyStripVfatH1_PlIId.begin(),NoisyStripVfatH1_PlIId.end(),NoisyVF)!=NoisyStripVfatH1_PlIId.end();
		  if(SaveForCorrel)
		    {
		      int positionInNtupla=find(NoisyStripVfatH1_PlIId.begin(),NoisyStripVfatH1_PlIId.end(),NoisyVF)-NoisyStripVfatH1_PlIId.begin();
		     H1NoiseNtupla[positionInNtupla]=1;
		     H1NoiseNtuplaS[positionInNtupla]=1;
		    }
		}		  
	      }
	    
	      if((vfpos==11)&&(pl==19))
		mult11_19=vfatMult;

	      if((vfpos==15)&&(pl==17))
		mult15_17=vfatMult;

	      if((vfpos==15)&&(pl==19))
		mult15_19=vfatMult;

	      if((vfpos==15)&&(pl==11)) VF15DispmultH1[0]=vfatMult;
	      if((vfpos==15)&&(pl==13)) VF15DispmultH1[1]=vfatMult;
	      if((vfpos==15)&&(pl==15)) VF15DispmultH1[2]=vfatMult;
	      if((vfpos==15)&&(pl==17)) VF15DispmultH1[3]=vfatMult;
	      if((vfpos==15)&&(pl==19)) VF15DispmultH1[4]=vfatMult;
	      if((vfpos==16)&&(pl==11)) VF16DispmultH1[0]=vfatMult;
	      if((vfpos==16)&&(pl==13)) VF16DispmultH1[1]=vfatMult;
	      if((vfpos==16)&&(pl==15)) VF16DispmultH1[2]=vfatMult;
	      if((vfpos==16)&&(pl==17)) VF16DispmultH1[3]=vfatMult;
	      if((vfpos==16)&&(pl==19)) VF16DispmultH1[4]=vfatMult;

	      if((vfpos==1)&&(pl==10)) VFMultPl10_VF15=vfatMult;
	      if((vfpos==0)&&(pl==10)) VFMultPl10_VF16=vfatMult;
	    }

	  
	}
      else
	std::cout<<"Warning:  found a T2 plane with symb > 39"<<std::endl;      
    }

  
  for(unsigned int r=0;r<5;r++){
    Vfat15_Pl10andH1Vfat15_PlDispariMultCorrel[r]->Fill(VFMultPl10_VF15,VF15DispmultH1.at(r));
    Vfat16_Pl10andH1Vfat16_PlDispariMultCorrel[r]->Fill(VFMultPl10_VF16,VF16DispmultH1.at(r));
  }
 
  NumH1VfatCompletelyOnPerEvt->Fill(numvfatComONH1);
  
  
  Vfat15_Pl19and17MultCorrel->Fill(mult15_17,mult15_19);
  Vfat15_Pl19and11_Pl19MultCorrel->Fill(mult11_19,mult15_19);

  if(NumTrk<5)
    Vfat15_Pl19and17MultCorrel_LowMulti->Fill(mult15_17,mult15_19);
  
  // for(unsigned int h=0; h<NoisyStripVfatH1_PlIId.size();h++){
  // EvtVfat_Strip_WhenCompletelyOn->Fill();
  // }

  //Save Ntuple-Noise-Correlation for H1. Only Interesting noisy vfat specified in NoisyStripVfatH1_PlIId Will be taken into account.
  //The Ntuple order of the vfat is important and it is ordered according to the vfat numbering : first is vfat 0 in plane 11.

  EvtVfat_Strip_WhenCompletelyOn->Fill(H1NoiseNtupla);
  EvtVfat_Strip_WhenCompletelyOnS->Fill(H1NoiseNtupla);
  // std::cout<<"Filling With"<<std::endl;
  /*
  for(unsigned int uu=0;uu<NoisyStripVfatH1_PlIId.size();uu++)
    std::cout<<H1NoiseNtupla[uu]<<"-";
  */
  //std::cout<<std::endl;

  int sizzz=NoisyStripVfatH1_PlIId.size();

  Double_t *rand = new Double_t[sizzz];
  EvtVfat_Strip_WhenCompletelyOnS->GetRandom(rand,false);
  //std::cout<<"Extracting:"<<std::endl;
  /*
  for(unsigned int uu=0;uu<NoisyStripVfatH1_PlIId.size();uu++)
    std::cout<<rand[uu]<<"-";
  
  std::cout<<std::endl;
  */
  delete[] H1NoiseNtupla;
  delete[] H1NoiseNtuplaS;

  }//if(VFATMonitoring)
  
  


  /*::::::::::Local variable Declaration:::::::::*/
  //double lmyx=0.0;
  //double lmyy=0.0;
  //double lmyr=0.0;
  //double phipart=0.0;

  
  //double zdetshift;
  //double zglobhit;

  //int zmm[2][10]; //2 row for halftelescope. 10 coloumn for planes
  //double z1=13828.3;          //14035.605; //first Gem first drift gas zone (mm)
  //double planedist= 86.0;
  //double btbdist=24.6;   //25.0;
  //double ovdist=43.0;
  //double zinsidedet=1.5;           //4.5; From first Drift zone to RO board = 9mm


  // **************************************************** STRIP STUDIES ************************************************ //
  // **************************************************** STRIP STUDIES ************************************************ //
  // **************************************************** STRIP STUDIES ************************************************ //
  // **************************************************** STRIP STUDIES ************************************************ //
  // **************************************************** STRIP STUDIES ************************************************ //

  /*
  for(T2StripClusterCollection::const_iterator itstrip = t2strclcoll->begin(); itstrip != t2strclcoll->end(); itstrip++){
  
    vector<T2Cluster> stripClv = itstrip->second;
    for(unsigned int k=0;k<stripClv.size();k++){
      clusterstripentries->Fill(stripClv[k].GetNoOfEntries());
    
    }
  }
  */
  
  
 for(T2StripClusterCollection::const_iterator itstrip = t2strclcoll->begin(); itstrip != t2strclcoll->end(); itstrip++){
      myT2Det = itstrip->first;
      myhalftele= myT2Det.halfTelescope();
      myplaneside= myT2Det.planeSide();
      myplane= myT2Det.plane();
      unsigned int myarm=myT2Det.arm();
      vector<T2Cluster> stripClv = itstrip->second;
      unsigned int  symb= myplane*2 + myplaneside + myhalftele*10 + myarm*20;
      for(unsigned int k=0;k<stripClv.size();k++){
	clusterstripentries->Fill(stripClv[k].GetNoOfEntries());
	rr=stripClv[k].GetClusterR();
	ff=stripClv[k].GetClusterPhi()*3.14159/180.;
	xx=rr*cos(ff);
	yy=rr*sin(ff);
	StripCluRadiographyXY[symb]->Fill(xx,yy);
      }  
     
  }

 // ****************************************************  pad STUDIES ************************************************ //
  // **************************************************** pad STUDIES ************************************************ //
  // **************************************************** pad STUDIES ************************************************ //
  // **************************************************** pad STUDIES ************************************************ //
  // **************************************************** pad STUDIES ************************************************ //

 unsigned int num11=0;
 unsigned int num16=0;
 
for(T2PadClusterCollection::const_iterator itpad = t2padclcoll->begin(); itpad != t2padclcoll->end(); itpad++){

      myT2Det = itpad->first;
      myhalftele= myT2Det.halfTelescope();
      myplaneside= myT2Det.planeSide();
      myplane= myT2Det.plane();
      unsigned int myarm=myT2Det.arm();

      vector<T2Cluster> padClv = itpad->second;
      unsigned int  symb= myplane*2 + myplaneside + myhalftele*10 + myarm*20;
      if(symb==11)
	num11=padClv.size();
      if(symb==16)
	num16=padClv.size();   
 	  
      for(unsigned int k=0;k<padClv.size();k++){
	
	rr=padClv[k].GetClusterR();
	ff=padClv[k].GetClusterPhi()*3.14159/180.;
	
        xx=rr*cos(ff);
	yy=rr*sin(ff);
	//std::cout<<" Woila'"<<std::endl;
	PlaneHitRadiographyXY[symb]->Fill(xx,yy);
	if(yy>0)
	  PadClusterXYRadiographyYPlus[symb]->Fill(xx,yy);
	else
	  PadClusterXYRadiographyYMinus[symb]->Fill(xx,yy);
	
	PadClusterR[symb]->Fill(rr);
	

	if(symb<10)
	  if(padClv[k].GetClusterPhi()>96)
	    if(padClv[k].GetClusterPhi()<264)
	      std::cout<<"Pippi-Bug is living in Pad Clu Phi="<<padClv[k].GetClusterPhi()<<" plane:"<<symb<<" entries:"<<padClv[k].GetNoOfEntries()<<std::endl;
	  
        clusterpadentries->Fill(padClv[k].GetNoOfEntries());
      }
     
}

 PlaneNumPadClu11_vs_16->Fill(num11,num16);
 



  // **************************************************** HIT STUDIES *********************************************** //
  // **************************************************** HIT STUDIES ************************************************ //
  // **************************************************** HIT STUDIES ************************************************ //
  // **************************************************** HIT STUDIES ************************************************ //
  // **************************************************** HIT STUDIES ************************************************ //

  // COMPARISON RECOHIT LOCALHIT
  // now we run for each plane the digi-simulation

 std::vector<unsigned int> numhitinplane;

for(unsigned int k=0;k<10;k++)
   {
     numhitinplane.push_back(0);
   }
 

 
for(T2HitCollection::const_iterator ithit = t2hitcoll->begin(); ithit != t2hitcoll->end(); ithit++){
		       
  bool trovato=false;
  for(unsigned int u=0;u<alldetZ.size();u++)
    {
      if(fabs(alldetZ.at(u)-(*ithit).GetHitZ())<0.001)
	trovato=true;
    }
  if(trovato==false)
    alldetZ.push_back((*ithit).GetHitZ());

  C12AllDetIdPlane->Fill((*ithit).GetHitPlane());
  
  AllHitZ->Fill((*ithit).GetHitZ());
  AllHitR->Fill((*ithit).GetHitR());
  
  if((*ithit).GetHitClass()==1)
    {
      AllHitZCL1->Fill((*ithit).GetHitZ());
      AllHitRCL1->Fill((*ithit).GetHitR());
    }
    
  if(((*ithit).GetHitNumStrip()>0)&&((*ithit).GetHitNumPad()==0))
    {
      AllHitZStripOnly->Fill((*ithit).GetHitZ());
      AllHitRStripOnly->Fill((*ithit).GetHitR());
    }

  if(((*ithit).GetHitNumPad()>0)&&((*ithit).GetHitNumStrip()==0))
    {
      AllHitZPadOnly->Fill((*ithit).GetHitZ());
      AllHitRPadOnly->Fill((*ithit).GetHitR());
    }
  
 
    
  C12AllDetIdPlaneS->Fill((*ithit).GetHitPlaneSide());

  unsigned int symb=(*ithit).GetHitPlane()*2+(*ithit).GetHitPlaneSide();
  
  numhitinplane.at(symb)++;

  

}

 
 //std::cout<<" Detector Z:  ";
for(unsigned int k=0;k<alldetZ.size();k++)
  {
    //std::cout<<alldetZ.at(k) <<" - ";
    meanz=meanz+alldetZ.at(k);
  }

// std::cout<<" MediaZ= "<<meanz<<std::endl; 

for(unsigned int k=0;k<10;k++)
  {
    NumAllHitvsIdAllevt->Fill(k,numhitinplane.at(k));
   
  }


//-----------------
//Digi Studies.Commented because problem in data saving by root
//------------------

// if(simufile){
  for(itp= PadDigiptr->begin(); itp!=PadDigiptr->end(); ++itp)
    {
    
      T2DetId mydet=(*itp).first;  
      
      unsigned int symb= RawtoSymb(mydet.calculateRawId(mydet.arm(),mydet.halfTelescope(),mydet.plane(),mydet.planeSide()));
      // symb= mydet.plane()*2 + mydet.planeSide() + mydet.arm()*20 + mydet.halfTelescope()*10;
      //std::cout<<" Symb: "<<symb<<std::endl;
      // if((symb<10)&&((symb/10)==SelectedHalf))
      if(((symb/10)==SelectedHalf))
	{
	//  unsigned int matrixnumb;
	  unsigned int occup=0;
	 
	  for(std::vector<T2PadDigi>::const_iterator itpad =(*itp).second.first; itpad !=(*itp).second.second; itpad++)
	    {
	      occup++;
	   /*   if((*itpad).getPadNr()==0)
		{
		  matrixnumb=(((*itpad).getRow()+1)+24*((*itpad).getCol()));   //Prima del dataraw
		  //	  std::cout<<" Caso-1"<<std::endl;
		}
	      else
		{
		  matrixnumb=(*itpad).getPadNr();	     
	  //	  std::cout<<" Caso-2"<<std::endl;
		}
	      //PAdNr[symb]->Fill(matrixnumb);
	      */
	    }
	   DigiPadOccupancy->Fill((symb%10),(occup/1560.0)*100.);
	
	}
    }

 for(its= StripDigiptr->begin(); its!=StripDigiptr->end(); ++its)
    {
    
      T2DetId mydet=(*its).first;  
      
      unsigned int symb= RawtoSymb(mydet.calculateRawId(mydet.arm(),mydet.halfTelescope(),mydet.plane(),mydet.planeSide()));
      // symb= mydet.plane()*2 + mydet.planeSide() + mydet.arm()*20 + mydet.halfTelescope()*10;
       
      if(((symb/10)==SelectedHalf))
	{
	  unsigned int occup=0;
//	  unsigned int matrixnumb;
	  
	  for(std::vector<T2StripDigi>::const_iterator itstrip =(*its).second.first; itstrip !=(*its).second.second; itstrip++)
	    {
	      occup++;
	      // if(symb==3)
	      //	std::cout<<"  "<<(*itstrip).getStripNr()<<" "<<(*itstrip).getRow()<<" "<<(*itstrip).getCol()<<std::endl;
		//
	      //matrixnumb=(((*itstrip).getRow())+256*((*itstrip).getCol()));   //Prima del dataraw	    
//	      matrixnumb=(*itstrip).getStripNr();
	      //StripNr[symb]->Fill(matrixnumb);
	    }	 
	  DigiStripOccupancy->Fill((symb%10),(occup/512.0)*100);
	}
    }
 //}

numevent++;

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



// ------------ method called once each job just before starting event loop  ------------
void T2AnalyzerRaw::beginJob()
{

  std::cout<<" Begin Job"<<std::endl;
  char *cmsswPath = getenv("CMSSW_BASE");
 
  numRsectEffi=8;
  numPhisectEffi=24;

  for(unsigned int i=0; i<40;i++){
     std::vector<std::vector<double> >  Geometry_PadEfficiency_Num1P;
     std::vector<std::vector<double> >  Geometry_PadEfficiency_Den1P;
     std::vector<std::vector<double> >  Geometry_PadEfficiency1P;
  
     std::vector<std::vector<double> >  Geometry_StripEfficiency_Num1P;
     std::vector<std::vector<double> >  Geometry_StripEfficiency_Den1P;
     std::vector<std::vector<double> >  Geometry_StripEfficiency1P;

     for( int lr=0;lr<numRsectEffi; lr++){
       std::vector<double> dummiphi;

       for( int lp=0;lp<numPhisectEffi; lp++){
	 dummiphi.push_back(0.0);
       }
       
       Geometry_PadEfficiency_Num1P.push_back(dummiphi);
       Geometry_PadEfficiency_Den1P.push_back(dummiphi);
       Geometry_PadEfficiency1P.push_back(dummiphi);
  
       Geometry_StripEfficiency_Num1P.push_back(dummiphi);
       Geometry_StripEfficiency_Den1P.push_back(dummiphi);
       Geometry_StripEfficiency1P.push_back(dummiphi);
     }
     
     Geometry_PadEfficiency_Num.push_back(Geometry_PadEfficiency_Num1P);
     Geometry_PadEfficiency_Den.push_back(Geometry_PadEfficiency_Den1P);
     Geometry_PadEfficiency.push_back(Geometry_PadEfficiency1P);
     
     Geometry_StripEfficiency_Num.push_back(Geometry_StripEfficiency_Num1P);
     Geometry_StripEfficiency_Den.push_back(Geometry_StripEfficiency_Den1P);
     Geometry_StripEfficiency.push_back(Geometry_StripEfficiency1P);
     
  }
std::cout<<" Geo Effi Vet Initialized"<<std::endl;

  //LoadDeadSectorFile
  //if (cmsswPath && skipEventFileName[0] != '/')      
  //skipEventFileName= string(cmsswPath) + string("/src/") + skipEventFileName;
  
if(DeadSectFileName != "")   
  if(cmsswPath && DeadSectFileName[0] != '/'){
         
      DeadSectFileName= string(cmsswPath) + string("/src/") + DeadSectFileName;
    
    std::string line;
    ifstream myfile(DeadSectFileName.c_str());

    if(myfile.is_open())
      {
	int detid=-1;
	int sectordead=-1; 
	if(verbosity)
	  std::cout<<" Event skipped: "<<std::endl;
	
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
		if(verbosity)
		  std::cout<<oneword<<std::endl;

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
      std::cout<<" Dead sectors file not loaded or not found"<<std::endl;
    
  }
  

  //Load the vector of the event to be skipped.
  if(skipSelectedEvents)
    {
      AllCorruptedEvents.clear();
      if (cmsswPath && skipEventFileName[0] != '/')      
	skipEventFileName= string(cmsswPath) + string("/src/") + skipEventFileName;
      
      std::string line;
      ifstream myfile(skipEventFileName.c_str());
		      
      if(myfile.is_open())
	{

	  if(verbosity)
	    std::cout<<" Event skipped: "<<std::endl;

	  while (! myfile.eof() )
	    {
	      getline (myfile,line);
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
		  
		  bool verbosity=false;
		  if(verbosity)
		    std::cout<<oneword<<std::endl;
		  
		  
		  int evttoskip = atoi(tokens.at(m).c_str());
		  AllCorruptedEvents.push_back(evttoskip);
		  
		}

	    }
	  
	}                 
    }
  //End of skipselected event;

  if(cmsswPath && xmlfilenameUsed_NotDead[0] != '/')
    xmlfilenameUsed_NotDead = string(cmsswPath) + string("/src/") + xmlfilenameUsed_NotDead;

  onemap=new DAQInformationSourceXML_a(xmlfilenameUsed_NotDead);
  Map_NOTDEAD=onemap->produceMap();
  
  //MapProducer_NotDeadUsedVfat = new DAQInformationSourceXML_a(xmlfilenameUsed_NotDead);
  // Map_NOTDEAD=MapProducer_NotDeadUsedVfat->produceMap();

  
  if(cmsswPath && xmlfilenameFull[0] != '/')
    xmlfilenameFull=string(cmsswPath) + string("/src/") + xmlfilenameFull;

  MapProducer_FullVFATs = new DAQInformationSourceXML_a(xmlfilenameFull);
  Map_FullVfats=MapProducer_FullVFATs->produceMap();
  
  
 for(unsigned int iii=0;iii<40;iii++) 
  for(unsigned int VFPOS0_16=0;VFPOS0_16<17;VFPOS0_16++)
    {
      
      unsigned int vfatsymb= iii*100 + VFPOS0_16;
     
	
      unsigned short printId = 0;

     
      
      map<unsigned int, VFATRegisters>::const_iterator dit1 = Map_NOTDEAD->readoutIdToRegisters.find(vfatsymb);	     
      if (dit1 != Map_NOTDEAD->readoutIdToRegisters.end()) {
	
	printId = dit1->second.GetFullChipID();  
	printf("@@@@ Class Not-Dead Map: Vfat Position:  ChipID=0x%x  vfatsymb=%d   \n",printId,vfatsymb);
	 
      }

      map<unsigned int, VFATRegisters>::const_iterator dit3 = Map_FullVfats->readoutIdToRegisters.find(vfatsymb);	     
      if (dit3 != Map_FullVfats->readoutIdToRegisters.end()) {
	
	printId = dit3->second.GetFullChipID();  
	printf("???? Class test Full Map: Vfat Position:  ChipID=0x%x  vfatsymb=%d   \n",printId,vfatsymb);
	 
      }

    }


  /*
  MapProducer_NotDeadUsedVfat = new DAQInformationSourceXML_a(xmlfilenameUsed_NotDead);
  Map_NotDeadUsedVfat=MapProducer_NotDeadUsedVfat->produceMap();
  map<unsigned int, VFATRegisters>::const_iterator mapIterator;
  for( mapIterator = Map_NotDeadUsedVfat->readoutIdToRegisters.begin(); mapIterator != Map_NotDeadUsedVfat->readoutIdToRegisters.end(); mapIterator++ ) {
    int symbvfid=mapIterator->first;
    vfatsSymbId.push_back(symbvfid);
  }
  std::cout<<" Found "<<vfatsSymbId.size()<<" vfats in the map utilized for efficiency caluclation"<<std::endl;
  */




  /*
// vfatsSymbId;
  map<unsigned int, VFATRegisters>::const_iterator dit = Map_NOTDEAD->readoutIdToRegisters.find(symbvfat);
  
	
  unsigned short printId = 0;
  if(dit != Map_NOTDEAD->readoutIdToRegisters.end()) 
    {
      toret=true;
    }
  else 
    {
      toret=false;
      
      std::string xmlfull;
      xmlfull="/home/mirko/SL/WorkingArea/CMSSW_311SVN2/CMSSW_3_1_1/src/TotemRawData/RawToDigi/python/T2GeoMapIP5_4quarter_vmea_cmssw.xml";
      onemapFull=new DAQInformationSourceXML_a(xmlfull);
  
      testmapNOTDEADFull=onemapFull->produceMap();
           
      map<unsigned int, VFATRegisters>::const_iterator ditFull = testmapFull->readoutIdToRegisters.find(symbvfat);
	
      unsigned short printId = 0;
      if (ditFull != testmapFull->readoutIdToRegisters.end()) {
	
	printId = ditFull->second.GetFullChipID();  		  
	printf("Vfat Position:  ChipID=0x%x  vfatsymb=%d   NOT MAPPED (plane %d, iid %d)\n",printId,symbvfat,symbvfat/100,symbvfat%100);		  
      }
      
    }
*/



  




  TH1::AddDirectory(kFALSE);
  //std::vector<double> alldetZ;
  meanz=0;
  //NumEvent 8 is the monitor 9 if you start from 0 so start from 1.
  numevent=1;
  matrixentries=0;
  trkcounter=0;
  countegood=0;

  for(unsigned int y=0;y<10;y++)
    {
      HITNoiseEnt.push_back(0);
      PADNoiseEnt.push_back(0);
      STRIPNoiseEnt.push_back(0);

      TotHITNoiseEnt.push_back(0);
      TotPADNoiseEnt.push_back(0);
      TotSTRIPNoiseEnt.push_back(0);      

    }
  
 
    NumEventNonCorruptedUtilized= std::auto_ptr<TH1F>(new TH1F("NumEventNonCorruptedUtilized","NumEvent without any non-dead vfat corruption",2, -0.5, 1.5));

    CluStripentriesGoodHit = std::auto_ptr<TH1F>(new TH1F("CluStripentriesGoodHit","Strip cluster-size (Cl1 hit, noise excluded from DR-D#phi cut)",20, -0.5, 19.5));
    CluStripentriesGoodHit->SetDirectory(0); 
    CluPadentriesGoodHit = std::auto_ptr<TH1F>(new TH1F("CluPadentriesGoodHit","Pad cluster-size (Cl1 hit, noise excluded from DR-D#phi cut)",20, -0.5, 19.5));
    CluPadentriesGoodHit->SetDirectory(0); 


    AllHitRPadOnly= std::auto_ptr<TH1F> (new TH1F("AllHitZPadOnly","AllHitZPadOnly",120, 30., 150.));
    AllHitZPadOnly= std::auto_ptr<TH1F> (new TH1F("AllHitZPadOnly","AllHitZPadOnly",2000,-14500,14500));
    AllHitRStripOnly= std::auto_ptr<TH1F> (new TH1F("AllHitRStripOnly","AllHitRStripOnly",120, 30., 150.));
    AllHitZStripOnly= std::auto_ptr<TH1F> (new TH1F("AllHitZStripOnly","AllHitZStripOnly",2000,-14500,14500)); 
    AllHitRCL1= std::auto_ptr<TH1F> (new TH1F("AllHitRCL1","AllHitRCL1",120, 30., 150.));
    AllHitZCL1= std::auto_ptr<TH1F> (new TH1F("AllHitZCL1","AllHitZCL1",2000,-14500,14500)); 

    AllHitZ=std::auto_ptr<TH1F>(new TH1F("AllHitZ","AllHitZ",2000,-14500,14500));
    AllHitZ->SetDirectory(0);
    AllHitR= std::auto_ptr<TH1F>(new TH1F("AllHitR","AllHitR",120, 30., 150.));
    AllHitR->SetDirectory(0); 


  // Chi2PhiProbLogy= std::auto_ptr<TCanvas>(new TCanvas("Chi2PhiProblog","Azimuthal #chi^{2} probability",400,400));
  // Chi2PhiProbLogy->SetDirectory(0);
  // Chi2RProbLogy= std::auto_ptr<TCanvas>(new TCanvas("Chi2RProblog","Radial #chi^{2} probability",400,400));
  //Chi2RProbLogy->SetDirectory(0);
  chiRredREF = std::auto_ptr<TH1F>(new TH1F("chiRredREF","Chi2ridotto riferimento",100, 0., 50.));
  clusterstripentries = std::auto_ptr<TH1F>(new TH1F("clusterstripentries","Number of strips in strip-clusters",20, -0.5, 19.5));
  clusterstripentries->SetDirectory(0); 
  clusterpadentries = std::auto_ptr<TH1F>(new TH1F("clusterpadentries","Number of pads in pad-clusters",20, -0.5, 19.5));
  clusterpadentries->SetDirectory(0); 
  diffphiGUNHIT= std::auto_ptr<TH1F>(new TH1F("diffphiGUNHIT","#Phi Generated particle - #Phi Reconstructed Hit",120, -30, 30));
  diffphiGUNHIT->SetDirectory(0); 
  RLocHRecoH= std::auto_ptr<TH1F>(new  TH1F("RLocHRecoH"," Reconstructed Hit R - GEANT4 Hit R",100,-1.5, 1.5));
  RLocHRecoH->SetDirectory(0); 
  diffRCluHit= std::auto_ptr<TH1F>(new  TH1F("diffRCluHit"," Cluster R - GEANT4 Hit R",100,-1.5, 1.5));
  diffRCluHit->SetDirectory(0); 
  Trketa= std::auto_ptr<TH1F>(new TH1F("Trketa","Reconstructed Track #eta (no cut used)",300,1.0,8.0));
  Trketa->SetDirectory(0); 
  Trkphi= std::auto_ptr<TH1F>(new TH1F("Trkphi","Reconstructed Track #phi (no cut used)",361,0,360));
  Trkphi->SetDirectory(0); 

    
  TrkphiALL= std::auto_ptr<TH1F>(new TH1F("TrkphiALL","Reconstructed Track #phi, ALL quarters",361,0,360));
  TrketaALL= std::auto_ptr<TH1F>(new TH1F("TrketaALL","Reconstructed Track #eta ",200,-20,20));
  NumhitinTrackALL= std::auto_ptr<TH1F>(new TH1F("NumhitinTrackALL","Number of hits forming the track",21,-0.5,20.5));
  TrkQuarterIdALL= std::auto_ptr<TH1F>(new TH1F("TrkQuarterIdALL","Id of the first plane in Track",40,-0.5,39.5));

  

  NumTrackALL_ONEQuarter= std::auto_ptr<TH1F>(new TH1F("NumTrackALL_ONEQuarter","NumTrackALL_ONEQuarter",151,-0.5,150.5));
  NumTrackALL_ONEQuarter->SetDirectory(0); 

  NumTrackALL_EveryQuarter= std::auto_ptr<TH1F>(new TH1F("NumTrackALL_EveryQuarter","NumTrackALL_EveryQuarter",151,-0.5,150.5));
  NumTrackALL_EveryQuarter->SetDirectory(0); 

  

  DPhiGoodTrk= std::auto_ptr<TH1F>(new TH1F("DPhiGoodTrk","#phi Track - #phi GEANT4 Track (Only primary tracks are considered)",51,-25.5,25.5));
  DPhiGoodTrk->SetDirectory(0); 
  DEtaGoodTrk= std::auto_ptr<TH1F>(new TH1F("DEtaGoodTrk","#eta Track - #phi GEANT4 Track (Only primary tracks are considered)",200,-1.0,1.0));
  DEtaGoodTrk->SetDirectory(0); 
  diffphiCluGun = std::auto_ptr<TH1F>(new TH1F("diffphiCluGun","Cluster #Phi - #Phi Generated particle",81, -40.5, 40.5));
  diffphiCluGun->SetDirectory(0); 
  diffphiCluHit= std::auto_ptr<TH1F>(new TH1F("diffphiCluHit","Cluster #Phi - GEANT4 Hit #Phi",50, -5, 5));
  diffphiCluHit->SetDirectory(0); 
  Trketagood= std::auto_ptr<TH1F>(new TH1F("Trketagood","Reconstructed Track #eta (#chi^{2}_{p} cut)",200,1.0,8.0));
  Trketagood->SetDirectory(0); 
  //Trkphigood= std::auto_ptr<TH1F>(new TH1F("Trkphigood","Reconstructed Track #phi (#chi^{2}_{p} cut)",361,0,360));
  Trkphigood= std::auto_ptr<TH1F>(new TH1F("Trkphigood","Reconstructed Track #phi",361,0,360));
  Trkphigood->SetDirectory(0); 

  TrkphiRZPlus= std::auto_ptr<TH1F>(new TH1F("TrkphiRZPlus","Reconstructed Track #phi",361,0,360));
  TrkphiRZMinus= std::auto_ptr<TH1F>(new TH1F("TrkphiRZMinus","Reconstructed Track #phi",361,0,360));

  Chi2RProb = std::auto_ptr<TH1F>(new TH1F("Chi2RProb","Radial #chi^{2} probability",15000,0.0,1.05));
  Chi2RProb->SetDirectory(0); 
  Chi2PhiProb = std::auto_ptr<TH1F>(new TH1F("Chi2PhiProb","Azimuthal #chi^{2} probability",15000,0.0,1.05));
  Chi2PhiProb->SetDirectory(0);

  IstoDisp= std::auto_ptr<TH1F>(new TH1F("IstoDisp","Displ.",100,0.0,10.0));
  IstoDisp->SetDirectory(0);
   
  DigiPadOccupancy= std::auto_ptr<TProfile>(new TProfile("DigiPadOccupancy","Plane Pad occupancy",10,-0.5,9.5,""));
  DigiPadOccupancy->GetYaxis()->SetTitle("Occupancy (%)");
  DigiPadOccupancy->GetXaxis()->SetTitle("Detector plane");
  DigiPadOccupancy->SetDirectory(0);

  DigiStripOccupancy= std::auto_ptr<TProfile>(new TProfile("DigiStripOccupancy","Plane Strip occupancy",10,-0.5,9.5,""));
  DigiStripOccupancy->GetYaxis()->SetTitle("Occupancy (%)");
  DigiStripOccupancy->GetXaxis()->SetTitle("Detector plane");
  DigiStripOccupancy->SetDirectory(0);

  BunchStatusBit= std::auto_ptr<TH2F>(new TH2F("BunchStatusBit","BunchStatusBit",1025,-0.5,1024.5,1025,-0.5,1024.5));
  BunchStatusBit->SetDirectory(0);

  Class1HitPadStripCLSCorrel= std::auto_ptr<TH2D>(new TH2D("Class1HitPadStripCLSCorrel","Pad-Strip Cluster size correlation (Class1 Hit)",60,-0.5,59.5,60,-0.5,59.5));
  Class1HitPadStripCLSCorrel->GetYaxis()->SetTitle("Pad CL Size");
  Class1HitPadStripCLSCorrel->GetXaxis()->SetTitle("Strip CL Size"); 
  Class1HitPadStripCLSCorrel->SetDirectory(0);  


  SingleParticleEfficiencyCut = std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiencyCut","# Tracks reconstructed vs generated particle #eta (all cuts)",34,4.35,7.75,""));
  SingleParticleEfficiencyCut->SetDirectory(0); 
  SingleParticleEfficiency = std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiency","# Tracks reconstructed vs generated particle #eta (all track)",34,4.35,7.75,""));
  SingleParticleEfficiency->SetDirectory(0); 
  MultiParticleEfficiencyCut= std::auto_ptr<TProfile>(new TProfile("MultiParticleEfficiencyCut","<# Tracks> vs # Charged Particles (all track cuts)",12,-0.5,11.5,""));
  MultiParticleEfficiencyCut->SetDirectory(0); 
  MultiParticleEfficiencyCutNorm= std::auto_ptr<TProfile>(new TProfile("MultiParticleEfficiencyCutNorm","<# Tracks / # Charged Particles > vs # Charged Particles (all track cuts)",12,-0.5,11.5,"")); 
  MultiParticleEfficiencyCutNorm->SetDirectory(0); 
   


  //CosmTracks= std::auto_ptr<TCanvas>(new TCanvas("c"));
//CosmTracks->SetDirectory(0); 
/*
CosmTracks0exl= std::auto_ptr<TCanvas>(new TCanvas("c0"));
CosmTracks1exl= std::auto_ptr<TCanvas>(new TCanvas("c1"));
CosmTracks2exl= std::auto_ptr<TCanvas>(new TCanvas("c2"));
CosmTracks3exl= std::auto_ptr<TCanvas>(new TCanvas("c3"));
CosmTracks4exl= std::auto_ptr<TCanvas>(new TCanvas("c4"));
*/

Tr4EffDet0= std::auto_ptr<TProfile>(new TProfile("Tr4EffDet0","Efficiency respect to detector 0",12,-1.5,10.5,""));
Tr4EffDet0->SetDirectory(0); 
Tr4EffDet0->SetXTitle("Detector plane");
Tr4EffDet1= std::auto_ptr<TProfile>(new TProfile("Tr4EffDet1","Efficiency respect to detector 1",12,-1.5,10.5,""));
Tr4EffDet1->SetDirectory(0); 
Tr4EffDet1->SetXTitle("Detector plane");
Tr4EffDet2= std::auto_ptr<TProfile>(new TProfile("Tr4EffDet2","Efficiency respect to detector 2",12,-1.5,10.5,""));
Tr4EffDet2->SetDirectory(0); 
Tr4EffDet2->SetXTitle("Detector plane");
Tr4EffDet3= std::auto_ptr<TProfile>(new TProfile("Tr4EffDet3","Efficiency respect to detector 3",12,-1.5,10.5,""));
Tr4EffDet3->SetDirectory(0); 
Tr4EffDet3->SetXTitle("Detector plane");
Tr4EffDet4= std::auto_ptr<TProfile>(new TProfile("Tr4EffDet4","Efficiency respect to detector 4",12,-1.5,10.5,""));
Tr4EffDet4->SetDirectory(0); 
Tr4EffDet4->SetXTitle("Detector plane");

AllDetIdPlane= std::auto_ptr<TH1F>(new TH1F("AllDetIdPlane","AllDetIdPlane",10,-0.5,9.5));
AllDetIdPlane->SetDirectory(0); 
AllDetIdPlaneS= std::auto_ptr<TH1F>(new TH1F("AllDetIdPlaneS","AllDetIdPlaneS",10,-0.5,9.5));
AllDetIdPlaneS->SetDirectory(0); 

C12AllDetIdPlane= std::auto_ptr<TH1F>(new TH1F("C12AllDetIdPlane","C12AllDetIdPlane",10,-0.5,9.5));
C12AllDetIdPlane->SetDirectory(0); 
C12AllDetIdPlaneS= std::auto_ptr<TH1F>(new TH1F("C12AllDetIdPlaneS","C12AllDetIdPlaneS",10,-0.5,9.5));
C12AllDetIdPlaneS->SetDirectory(0); 


NumAllHitvsIdAllevt= std::auto_ptr<TProfile>(new TProfile("NumAllHitvsIdAllevt","NumAllHitvsIdAllevt",10,-0.5,9.5));
NumAllHitvsIdAllevt->SetDirectory(0); 
NumAllHitvsIdGoodevt= std::auto_ptr<TProfile>(new TProfile("NumAllHitvsIdGoodevt","NumAllHitvsIdGoodevt",10,-0.5,9.5));
 NumAllHitvsIdGoodevt->SetDirectory(0);


 AllDetId= std::auto_ptr<TH1F>(new TH1F("AllDetId","AllDetId",10,-0.5,9.5));
 AllDetId->SetDirectory(0); 

 drvsIddet= std::auto_ptr<TProfile>(new TProfile("drvsIddet","|dr| medio vs id",12,-1.5,10.5,""));
 drvsIddet->SetDirectory(0); 
 alldetdr= std::auto_ptr<TH1F>(new TH1F("alldetdr","dr hit",100,-5.0,5.0));
 alldetdr->SetDirectory(0); 
 dphivsIddet= std::auto_ptr<TProfile>(new TProfile("dphivsIddet","|dphi| medio vs id",12,-1.5,10.5,""));
 dphivsIddet->SetDirectory(0); 
 alldetdphi= std::auto_ptr<TH1F>(new TH1F("alldetdphi","dphi hit",100,-10.0,10.0));
 alldetdphi->SetDirectory(0); 

Tr4EffDetAll= std::auto_ptr<TProfile>(new TProfile("Tr4EffDetAll","Efficiency of detectors (Class 1 Hit required), 4 hit tracks",12,-1.5,10.5,""));
Tr4EffDetAll->SetDirectory(0); 
Tr4EffDetAll->SetXTitle("Detector plane");
//Tr4EffDetAll->SetDirectory(0); 

TrnEffDetAll=  std::auto_ptr<TProfile>(new TProfile(" TrnEffDetAll","Efficiency of detectors, n hit tracks",12,-1.5,10.5,""));
TrnEffDetAll->SetDirectory(0);
TrnEffDetAll->SetXTitle("Detector plane");
//Trnm1EffDetAll->SetDirectory(0);
Trnm1EffDetAll=  std::auto_ptr<TProfile>(new TProfile(" Trnm1EffDetAll","Efficiency of detectors, n-1 hit tracks",12,-1.5,10.5,""));
Trnm1EffDetAll->SetDirectory(0);
Trnm1EffDetAll->SetXTitle("Detector plane");
//Trnm1EffDetAll->SetDirectory(0);



TrnEffDetAllpad=std::auto_ptr<TProfile>(new TProfile("TrnEffDetAllpad","Efficiency (pad) of detectors, n hit tracks",12,-1.5,10.5,""));
TrnEffDetAllpad->SetDirectory(0);
Trnm1EffDetAllpad=  std::auto_ptr<TProfile>(new TProfile("Trnm1EffDetAllstrip","Efficiency (strip) of detectors, n-1 hit tracks",12,-1.5,10.5,""));
Trnm1EffDetAllpad->SetDirectory(0);
TrnEffDetAllstrip=std::auto_ptr<TProfile>(new TProfile("TrnEffDetAllpad","Efficiency (pad) of detectors, n hit tracks",12,-1.5,10.5,""));
TrnEffDetAllstrip->SetDirectory(0);
Trnm1EffDetAllstrip=  std::auto_ptr<TProfile>(new TProfile("Trnm1EffDetAllstrip","Efficiency (strip) of detectors, n-1 hit tracks",12,-1.5,10.5,""));
Trnm1EffDetAllstrip->SetDirectory(0);





Tr4EffDet0pad= std::auto_ptr<TProfile>(new TProfile("Tr4EffDet0pad","Efficiency respect to detector 0",12,-1.5,10.5,""));
Tr4EffDet0pad->SetDirectory(0); 
Tr4EffDet0pad->SetXTitle("Detector plane");
Tr4EffDet1pad= std::auto_ptr<TProfile>(new TProfile("Tr4EffDet1pad","Efficiency respect to detector 1",12,-1.5,10.5,""));
Tr4EffDet1pad->SetDirectory(0); 
Tr4EffDet1pad->SetXTitle("Detector plane");
Tr4EffDet2pad= std::auto_ptr<TProfile>(new TProfile("Tr4EffDet2pad","Efficiency respect to detector 2",12,-1.5,10.5,""));
Tr4EffDet2pad->SetDirectory(0); 
Tr4EffDet2pad->SetXTitle("Detector plane");
Tr4EffDet3pad= std::auto_ptr<TProfile>(new TProfile("Tr4EffDet3pad","Efficiency respect to detector 3",12,-1.5,10.5,""));
Tr4EffDet3pad->SetDirectory(0); 
Tr4EffDet3pad->SetXTitle("Detector plane");
Tr4EffDet4pad= std::auto_ptr<TProfile>(new TProfile("Tr4EffDet4pad","Efficiency respect to detector 4",12,-1.5,10.5,""));
Tr4EffDet4pad->SetDirectory(0); 
Tr4EffDet4pad->SetXTitle("Detector plane");
Tr4EffDetAllpad= std::auto_ptr<TProfile>(new TProfile("Tr4EffDetAllpad","Efficiency of detectors (Pad-Cluster required), 4 hit tracks",12,-1.5,10.5,""));
Tr4EffDetAllpad->SetDirectory(0); 
Tr4EffDetAllpad->SetXTitle("Detector plane");
//Tr4EffDetAllpad->SetDirectory(0); 


Tr4EffDet0strip= std::auto_ptr<TProfile>(new TProfile("Tr4EffDet0strip","Efficiency respect to detector 0",12,-1.5,10.5,""));
Tr4EffDet0strip->SetDirectory(0); 
Tr4EffDet0strip->SetXTitle("Detector plane");
Tr4EffDet1strip= std::auto_ptr<TProfile>(new TProfile("Tr4EffDet1strip","Efficiency respect to detector 1",12,-1.5,10.5,""));
Tr4EffDet1strip->SetDirectory(0); 
Tr4EffDet1strip->SetXTitle("Detector plane");
Tr4EffDet2strip= std::auto_ptr<TProfile>(new TProfile("Tr4EffDet2strip","Efficiency respect to detector 2",12,-1.5,10.5,""));
Tr4EffDet2strip->SetDirectory(0); 
Tr4EffDet2strip->SetXTitle("Detector plane");
Tr4EffDet3strip= std::auto_ptr<TProfile>(new TProfile("Tr4EffDet3strip","Efficiency respect to detector 3",12,-1.5,10.5,""));
Tr4EffDet3strip->SetDirectory(0); 
Tr4EffDet3strip->SetXTitle("Detector plane");
Tr4EffDet4strip= std::auto_ptr<TProfile>(new TProfile("Tr4EffDet4strip","Efficiency respect to detector 4",12,-1.5,10.5,""));
Tr4EffDet4strip->SetDirectory(0); 
Tr4EffDet4strip->SetXTitle("Detector plane");
Tr4EffDetAllstrip= std::auto_ptr<TProfile>(new TProfile("Tr4EffDetAllstrip","Efficiency of detectors (Strip-Cluster required), 4 hit tracks",12,-1.5,10.5,""));
Tr4EffDetAllstrip->SetDirectory(0); 
Tr4EffDetAllstrip->SetXTitle("Detector plane");
//Tr4EffDetAllstrip->SetDirectory(0);




Tr3EffDet0= std::auto_ptr<TProfile>(new TProfile("Tr3EffDet0","Efficiency respect to detector 0",12,-1.5,10.5,""));
Tr3EffDet0->SetDirectory(0); 
Tr3EffDet0->SetXTitle("Detector plane");
Tr3EffDet1= std::auto_ptr<TProfile>(new TProfile("Tr3EffDet1","Efficiency respect to detector 1",12,-1.5,10.5,""));
Tr3EffDet1->SetDirectory(0); 
Tr3EffDet1->SetXTitle("Detector plane");
Tr3EffDet2= std::auto_ptr<TProfile>(new TProfile("Tr3EffDet2","Efficiency respect to detector 2",12,-1.5,10.5,""));
Tr3EffDet2->SetDirectory(0); 
Tr3EffDet2->SetXTitle("Detector plane");
Tr3EffDet3= std::auto_ptr<TProfile>(new TProfile("Tr3EffDet3","Efficiency respect to detector 3",12,-1.5,10.5,""));
Tr3EffDet3->SetDirectory(0); 
Tr3EffDet3->SetXTitle("Detector plane");
Tr3EffDet4= std::auto_ptr<TProfile>(new TProfile("Tr3EffDet4","Efficiency respect to detector 4",12,-1.5,10.5,""));
Tr3EffDet4->SetDirectory(0); 
Tr3EffDet4->SetXTitle("Detector plane");
Tr3EffDetAll= std::auto_ptr<TProfile>(new TProfile("Tr3EffDetAll","Efficiency of detectors (Class 1 Hit required), 3 hit tracks",12,-1.5,10.5,""));
Tr3EffDetAll->SetDirectory(0); 
Tr3EffDetAll->SetXTitle("Detector plane");
//Tr3EffDetAll->SetDirectory(0); 

Tr3EffDet0pad= std::auto_ptr<TProfile>(new TProfile("Tr3EffDet0pad","Efficiency respect to detector 0",12,-1.5,10.5,""));
Tr3EffDet0pad->SetDirectory(0); 
Tr3EffDet0pad->SetXTitle("Detector plane");
Tr3EffDet1pad= std::auto_ptr<TProfile>(new TProfile("Tr3EffDet1pad","Efficiency respect to detector 1",12,-1.5,10.5,""));
Tr3EffDet1pad->SetDirectory(0); 
Tr3EffDet1pad->SetXTitle("Detector plane");
Tr3EffDet2pad= std::auto_ptr<TProfile>(new TProfile("Tr3EffDet2pad","Efficiency respect to detector 2",12,-1.5,10.5,""));
Tr3EffDet2pad->SetDirectory(0); 
Tr3EffDet2pad->SetXTitle("Detector plane");
Tr3EffDet3pad= std::auto_ptr<TProfile>(new TProfile("Tr3EffDet3pad","Efficiency respect to detector 3",12,-1.5,10.5,""));
Tr3EffDet3pad->SetDirectory(0); 
Tr3EffDet3pad->SetXTitle("Detector plane");
Tr3EffDet4pad= std::auto_ptr<TProfile>(new TProfile("Tr3EffDet4pad","Efficiency respect to detector 4",12,-1.5,10.5,""));
Tr3EffDet4pad->SetDirectory(0); 
Tr3EffDet4pad->SetXTitle("Detector plane");
Tr3EffDetAllpad= std::auto_ptr<TProfile>(new TProfile("Tr3EffDetAllpad","Efficiency of detectors (Pad-Cluster required), 3 hit tracks",12,-1.5,10.5,""));
Tr3EffDetAllpad->SetDirectory(0); 
Tr3EffDetAllpad->SetXTitle("Detector plane");
//Tr3EffDetAllpad->SetDirectory(0); 


Tr3EffDet0strip= std::auto_ptr<TProfile>(new TProfile("Tr3EffDet0strip","Efficiency respect to detector 0",12,-1.5,10.5,""));
Tr3EffDet0strip->SetDirectory(0); 
Tr3EffDet0strip->SetXTitle("Detector plane");
Tr3EffDet1strip= std::auto_ptr<TProfile>(new TProfile("Tr3EffDet1strip","Efficiency respect to detector 1",12,-1.5,10.5,""));
Tr3EffDet1strip->SetDirectory(0); 
Tr3EffDet1strip->SetXTitle("Detector plane");
Tr3EffDet2strip= std::auto_ptr<TProfile>(new TProfile("Tr3EffDet2strip","Efficiency respect to detector 2",12,-1.5,10.5,""));
Tr3EffDet2strip->SetDirectory(0); 
Tr3EffDet2strip->SetXTitle("Detector plane");
Tr3EffDet3strip= std::auto_ptr<TProfile>(new TProfile("Tr3EffDet3strip","Efficiency respect to detector 3",12,-1.5,10.5,""));
Tr3EffDet3strip->SetDirectory(0); 
Tr3EffDet3strip->SetXTitle("Detector plane");
Tr3EffDet4strip= std::auto_ptr<TProfile>(new TProfile("Tr3EffDet4strip","Efficiency respect to detector 4",12,-1.5,10.5,""));
Tr3EffDet4strip->SetDirectory(0); 
Tr3EffDet4strip->SetXTitle("Detector plane");
Tr3EffDetAllstrip= std::auto_ptr<TProfile>(new TProfile("Tr3EffDetAllstrip","Efficiency of detectors (Strip-Cluster required), 3 hit tracks",12,-1.5,10.5,""));
Tr3EffDetAllstrip->SetDirectory(0); 
Tr3EffDetAllstrip->SetXTitle("Detector plane");



Tr2EffDet0= std::auto_ptr<TProfile>(new TProfile("Tr2EffDet0","Efficiency respect to detector 0",12,-1.5,10.5,""));
Tr2EffDet0->SetDirectory(0); 
Tr2EffDet0->SetXTitle("Detector plane");
Tr2EffDet1= std::auto_ptr<TProfile>(new TProfile("Tr2EffDet1","Efficiency respect to detector 1",12,-1.5,10.5,""));
Tr2EffDet1->SetDirectory(0); 
Tr2EffDet1->SetXTitle("Detector plane");
Tr2EffDet2= std::auto_ptr<TProfile>(new TProfile("Tr2EffDet2","Efficiency respect to detector 2",12,-1.5,10.5,""));
Tr2EffDet2->SetDirectory(0); 
Tr2EffDet2->SetXTitle("Detector plane");
Tr2EffDet3= std::auto_ptr<TProfile>(new TProfile("Tr2EffDet3","Efficiency respect to detector 3",12,-1.5,10.5,""));
Tr2EffDet3->SetDirectory(0); 
Tr2EffDet3->SetXTitle("Detector plane");
Tr2EffDet4= std::auto_ptr<TProfile>(new TProfile("Tr2EffDet4","Efficiency respect to detector 4",12,-1.5,10.5,""));
Tr2EffDet4->SetDirectory(0); 
Tr2EffDet4->SetXTitle("Detector plane");
Tr2EffDetAll= std::auto_ptr<TProfile>(new TProfile("Tr2EffDetAll","Efficiency of detectors",12,-1.5,10.5,""));
Tr2EffDetAll->SetDirectory(0); 
Tr2EffDetAll->SetXTitle("Detector plane");
//Tr2EffDetAll->SetDirectory(0); 

Tr2EffDet0pad= std::auto_ptr<TProfile>(new TProfile("Tr2EffDet0pad","Efficiency respect to detector 0",12,-1.5,10.5,""));
Tr2EffDet0pad->SetDirectory(0); 
Tr2EffDet0pad->SetXTitle("Detector plane");
Tr2EffDet1pad= std::auto_ptr<TProfile>(new TProfile("Tr2EffDet1pad","Efficiency respect to detector 1",12,-1.5,10.5,""));
Tr2EffDet1pad->SetDirectory(0); 
Tr2EffDet1pad->SetXTitle("Detector plane");
Tr2EffDet2pad= std::auto_ptr<TProfile>(new TProfile("Tr2EffDet2pad","Efficiency respect to detector 2",12,-1.5,10.5,""));
Tr2EffDet2pad->SetDirectory(0); 
Tr2EffDet2pad->SetXTitle("Detector plane");
Tr2EffDet3pad= std::auto_ptr<TProfile>(new TProfile("Tr2EffDet3pad","Efficiency respect to detector 3",12,-1.5,10.5,""));
Tr2EffDet3pad->SetDirectory(0); 
Tr2EffDet3pad->SetXTitle("Detector plane");
Tr2EffDet4pad= std::auto_ptr<TProfile>(new TProfile("Tr2EffDet4pad","Efficiency respect to detector 4",12,-1.5,10.5,""));
Tr2EffDet4pad->SetDirectory(0); 
Tr2EffDet4pad->SetXTitle("Detector plane");
Tr2EffDetAllpad= std::auto_ptr<TProfile>(new TProfile("Tr2EffDetAllpad","Efficiency of detectors",12,-1.5,10.5,""));
Tr2EffDetAllpad->SetDirectory(0); 
Tr2EffDetAllpad->SetXTitle("Detector plane");
//Tr2EffDetAllpad->SetDirectory(0); 


Tr2EffDet0strip= std::auto_ptr<TProfile>(new TProfile("Tr2EffDet0strip","Efficiency respect to detector 0",12,-1.5,10.5,""));
Tr2EffDet0strip->SetDirectory(0); 
Tr2EffDet0strip->SetXTitle("Detector plane");
Tr2EffDet1strip= std::auto_ptr<TProfile>(new TProfile("Tr2EffDet1strip","Efficiency respect to detector 1",12,-1.5,10.5,""));
Tr2EffDet1strip->SetDirectory(0); 
Tr2EffDet1strip->SetXTitle("Detector plane");
Tr2EffDet2strip= std::auto_ptr<TProfile>(new TProfile("Tr2EffDet2strip","Efficiency respect to detector 2",12,-1.5,10.5,""));
Tr2EffDet2strip->SetDirectory(0); 
Tr2EffDet2strip->SetXTitle("Detector plane");
Tr2EffDet3strip= std::auto_ptr<TProfile>(new TProfile("Tr2EffDet3strip","Efficiency respect to detector 3",12,-1.5,10.5,""));
Tr2EffDet3strip->SetDirectory(0); 
Tr2EffDet3strip->SetXTitle("Detector plane");
Tr2EffDet4strip= std::auto_ptr<TProfile>(new TProfile("Tr2EffDet4strip","Efficiency respect to detector 4",12,-1.5,10.5,""));
Tr2EffDet4strip->SetDirectory(0); 
Tr2EffDet4strip->SetXTitle("Detector plane");
Tr2EffDetAllstrip= std::auto_ptr<TProfile>(new TProfile("Tr2EffDetAllstrip","Efficiency of detectors",12,-1.5,10.5,""));
Tr2EffDetAllstrip->SetDirectory(0); 
Tr2EffDetAllstrip->SetXTitle("Detector plane");
//Tr2EffDetAllstrip->SetDirectory(0);

RandomGauss=std::auto_ptr<TH1F>(new TH1F("RandomGauss","RandomGauss",1000,-5000,5000));
RandomGauss->SetDirectory(0);




Chi2phired=std::auto_ptr<TH1F>(new TH1F("Chi2phired","Reduced #chi^{2}-#Phi",200,0.,50.));
Chi2phired->SetDirectory(0);

Chi2red=std::auto_ptr<TH1F>(new TH1F("Chi2red","Reduced #chi^{2}-R",200,0.,50.));
Chi2red->SetDirectory(0);

ActiveSymbDet=std::auto_ptr<TH1F>(new TH1F("ActiveSymbDet","Symb det inside tracks",20,-0.5,19.5));
ActiveSymbDet->SetDirectory(0);
PolarAngles=std::auto_ptr<TH1F>(new TH1F("PolarAngles","Track #theta_{x} angle (rad)",20,-0.25,0.25));
PolarAngles->SetDirectory(0);

NumhitinTrack=std::auto_ptr<TH1F>(new TH1F("NumhitinTrack","# hit in Track, in the analyzed quarter",18,-0.5,17.5));
NumhitinTrack->SetDirectory(0);
NumhitinTrackAligned=std::auto_ptr<TH1F>(new TH1F("NumhitinTrackAligned","# hit in Matrix-Track Aligned",18,-0.5,17.5));
NumhitinTrackAligned->SetDirectory(0);

NumhitinTrackGood=std::auto_ptr<TH1F>(new TH1F("NumhitinTrackGood","# hit in Track",12,-0.5,11.5));
NumhitinTrackGood->SetDirectory(0);

TracksplanestatN= std::auto_ptr<TH1F>(new TH1F("TracksplanestatN","Number of reference tracks used for the efficiency calculation",12,-1.5,10.5));
TracksplanestatNm1= std::auto_ptr<TH1F>(new TH1F("TracksplanestatNm1","Number of tracks candidate for efficiency of x axis plane",12,-1.5,10.5));
Tracksplanestat3= std::auto_ptr<TH1F>(new TH1F("Tracksplanestat3","Number of tracks candidate for efficiency of x axis plane",12,-1.5,10.5));
TracksplanestatN->SetXTitle("Detector Id");
Tracksplanestat3->SetDirectory(0);
Tracksplanestat3->SetXTitle("plane");

Tracksplanestat4= std::auto_ptr<TH1F>(new TH1F("Tracksplanestat4","Number of tracks candidate for efficiency of x axis plane",12,-1.5,10.5));

Tracksplanestat4->SetDirectory(0);
Tracksplanestat4->SetXTitle("plane");

 EventHemisphere= std::auto_ptr<TH1D>(new TH1D("EventHemisphere","EventHemisphere",4,-1.5,2.5)); 
  Count_t2trackVectorALL_H0= std::auto_ptr<TH1D>(new TH1D("Count_t2trackVectorALL_H0","Count_t2trackVectorALL_H0",200,-0.5,199.5)); 
  Count_t2trackVectorALL_H1= std::auto_ptr<TH1D>(new TH1D("Count_t2trackVectorALL_H1","Count_t2trackVectorALL_H1",200,-0.5,199.5)); 
  Count_t2trackVectorALL_H2= std::auto_ptr<TH1D>(new TH1D("Count_t2trackVectorALL_H2","Count_t2trackVectorALL_H2",200,-0.5,199.5)); 
  Count_t2trackVectorALL_H3= std::auto_ptr<TH1D>(new TH1D("Count_t2trackVectorALL_H3","Count_t2trackVectorALL_H3",200,-0.5,199.5)); 
    


 NumPadCluVsPlane=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlane","<# Pad Clu> vs Plane ",12,-1.5,10.5)); 
 NumStripCluVsPlane= std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlane"," <# Strip Clu> vs Plane ",12,-1.5,10.5)); 
 NumCl1HitVsPlane=std::auto_ptr<TProfile> (new TProfile("NumCl1HitVsPlane","<# Cl1Hit> vs Plane",12,-1.5,10.5)); 
 PadCluSizeVsPlaneAll=std::auto_ptr<TProfile> (new TProfile("PadCluSizeVsPlaneAll","Pad Cluster Size vs Plane",12,-1.5,10.5)); 
 StipCluSizeVsPlaneAll=std::auto_ptr<TProfile> (new TProfile("StipCluSizeVsPlaneAll","Strip Cluster Size vs Plane",12,-1.5,10.5));

 PadCluSizeVsPlaneAll2=std::auto_ptr<TProfile> (new TProfile("PadCluSizeVsPlaneAll2","Pad Cluster Size vs Plane from vect",12,-1.5,10.5));
  PadCluMultiplVsPlane= std::auto_ptr<TProfile>(new TProfile("PadCluMultiplVsPlane","Pad Cluster multiplicity vs Plane leftright req",40,-0.5,39.5,""));
  PadCluMultiplVsPlane->SetDirectory(0);



NumPadCluVsPlaneAll3H0=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H0","Num Pad Cluster  vs Plane -cluinfo- H0",12,-1.5,10.5));
 PadCluSizeVsPlaneAll3H0=std::auto_ptr<TProfile> (new TProfile("PadCluSizeVsPlaneAll3H0","Pad Cluster Size vs Plane -cluinfo- H0",12,-1.5,10.5));

NumPadCluVsPlaneAll3H1=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H1","Num Pad Cluster  vs Plane -cluinfo- H1",12,-1.5,10.5));
 PadCluSizeVsPlaneAll3H1=std::auto_ptr<TProfile> (new TProfile("PadCluSizeVsPlaneAll3H1","Pad Cluster Size vs Plane -cluinfo- H1",12,-1.5,10.5));

NumPadCluVsPlaneAll3H2=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H2","Num Pad Cluster  vs Plane -cluinfo- H2",12,-1.5,10.5));
 PadCluSizeVsPlaneAll3H2=std::auto_ptr<TProfile> (new TProfile("PadCluSizeVsPlaneAll3H2","Pad Cluster Size vs Plane -cluinfo- H2",12,-1.5,10.5));

 NumPadCluVsPlaneAll3H3=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H3","Num Pad Cluster  vs Plane -cluinfo- H3",12,-1.5,10.5));
 PadCluSizeVsPlaneAll3H3=std::auto_ptr<TProfile> (new TProfile("PadCluSizeVsPlaneAll3H3","Pad Cluster Size vs Plane -cluinfo- H3",12,-1.5,10.5));


NumPadCluVsPlaneAll3H0_Cutted_LOWMultipl=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H0_Cutted_LOWMultipl","Pad Cluster Size vs Plane -cluinfo- H0",12,-1.5,10.5));
NumPadCluVsPlaneAll3H1_Cutted_LOWMultipl=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H1_Cutted_LOWMultipl","Pad Cluster Size vs Plane -cluinfo- H1",12,-1.5,10.5));
NumPadCluVsPlaneAll3H2_Cutted_LOWMultipl=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H2_Cutted_LOWMultipl","Pad Cluster Size vs Plane -cluinfo- H2",12,-1.5,10.5));
NumPadCluVsPlaneAll3H3_Cutted_LOWMultipl=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H3_Cutted_LOWMultipl","Pad Cluster Size vs Plane -cluinfo- H3",12,-1.5,10.5));




NumStripCluVsPlaneAll3H0=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H0","Num Strip Cluster  vs Plane -cluinfo- H0",12,-1.5,10.5));
 StripCluSizeVsPlaneAll3H0=std::auto_ptr<TProfile> (new TProfile("StripCluSizeVsPlaneAll3H0","Strip Cluster Size vs Plane -cluinfo- H0",12,-1.5,10.5));

NumStripCluVsPlaneAll3H1=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H1","Num Strip Cluster  vs Plane -cluinfo- H1",12,-1.5,10.5));
 StripCluSizeVsPlaneAll3H1=std::auto_ptr<TProfile> (new TProfile("StripCluSizeVsPlaneAll3H1","Strip Cluster  vs Plane -cluinfo- H1",12,-1.5,10.5));

NumStripCluVsPlaneAll3H2=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H2","Num Strip Cluster Size vs Plane -cluinfo- H2",12,-1.5,10.5));
 StripCluSizeVsPlaneAll3H2=std::auto_ptr<TProfile> (new TProfile("StripCluSizeVsPlaneAll3H2","Strip Cluster  vs Plane -cluinfo- H2",12,-1.5,10.5));

NumStripCluVsPlaneAll3H3=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H3","Num Strip Cluster Size vs Plane -cluinfo- H3",12,-1.5,10.5));
 StripCluSizeVsPlaneAll3H3=std::auto_ptr<TProfile> (new TProfile("StripCluSizeVsPlaneAll3H3","Strip Cluster  vs Plane -cluinfo- H3",12,-1.5,10.5));


NumPadCluVsPlaneAll3H0_Cutted=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H0_Cutted","Num Pad Cluster  vs Plane -cluinfo- H0 cls<5",12,-1.5,10.5));
NumPadCluVsPlaneAll3H1_Cutted=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H1_Cutted","Num Pad Cluster  vs Plane -cluinfo- H1 cls<5",12,-1.5,10.5));
NumPadCluVsPlaneAll3H2_Cutted=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H2_Cutted","Num Pad Cluster  vs Plane -cluinfo- H2 cls<5",12,-1.5,10.5));
NumPadCluVsPlaneAll3H3_Cutted=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H3_Cutted","Num Pad Cluster  vs Plane -cluinfo- H3 cls<5",12,-1.5,10.5));


NumStripCluVsPlaneAll3H0_Cutted=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H0_Cutted","Num Strip Cluster  vs Plane -cluinfo- H0 cls<5",12,-1.5,10.5));
NumStripCluVsPlaneAll3H1_Cutted=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H1_Cutted","Num Strip Cluster  vs Plane -cluinfo- H1 cls<5",12,-1.5,10.5));
NumStripCluVsPlaneAll3H2_Cutted=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H2_Cutted","Num Strip Cluster  vs Plane -cluinfo- H2 cls<5",12,-1.5,10.5));
NumStripCluVsPlaneAll3H3_Cutted=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H3_Cutted","Num Strip Cluster  vs Plane -cluinfo- H3 cls<5",12,-1.5,10.5));




NumPadCluVsPlaneAll3H0_CuttedPlus=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H0_CuttedPlus","Num Pad Cluster  vs Plane -cluinfo- H0 cls<5",12,-1.5,10.5));
NumPadCluVsPlaneAll3H1_CuttedPlus=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H1_CuttedPlus","Num Pad Cluster  vs Plane -cluinfo- H1 cls<5",12,-1.5,10.5));
NumPadCluVsPlaneAll3H2_CuttedPlus=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H2_CuttedPlus","Num Pad Cluster  vs Plane -cluinfo- H2 cls<5",12,-1.5,10.5));
NumPadCluVsPlaneAll3H3_CuttedPlus=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H3_CuttedPlus","Num Pad Cluster  vs Plane -cluinfo- H3 cls<5",12,-1.5,10.5));
NumPadCluVsPlaneAll3H0_CuttedMinus=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H0_CuttedMinus","Num Pad Cluster  vs Plane -cluinfo- H0 cls<5",12,-1.5,10.5));
NumPadCluVsPlaneAll3H1_CuttedMinus=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H1_CuttedMinus","Num Pad Cluster  vs Plane -cluinfo- H1 cls<5",12,-1.5,10.5));
NumPadCluVsPlaneAll3H2_CuttedMinus=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H2_CuttedMinus","Num Pad Cluster  vs Plane -cluinfo- H2 cls<5",12,-1.5,10.5));
NumPadCluVsPlaneAll3H3_CuttedMinus=std::auto_ptr<TProfile> (new TProfile("NumPadCluVsPlaneAll3H3_CuttedMinus","Num Pad Cluster  vs Plane -cluinfo- H3 cls<5",12,-1.5,10.5));



NumStripCluVsPlaneAll3H0_CuttedPlus=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H0_CuttedPlus","Num Strip Cluster  vs Plane -cluinfo- H0 cls<5",12,-1.5,10.5));
NumStripCluVsPlaneAll3H1_CuttedPlus=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H1_CuttedPlus","Num Strip Cluster  vs Plane -cluinfo- H1 cls<5",12,-1.5,10.5));
NumStripCluVsPlaneAll3H2_CuttedPlus=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H2_CuttedPlus","Num Strip Cluster  vs Plane -cluinfo- H2 cls<5",12,-1.5,10.5));
NumStripCluVsPlaneAll3H3_CuttedPlus=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H3_CuttedPlus","Num Strip Cluster  vs Plane -cluinfo- H3 cls<5",12,-1.5,10.5));
NumStripCluVsPlaneAll3H0_CuttedMinus=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H0_CuttedMinus","Num Strip Cluster  vs Plane -cluinfo- H0 cls<5",12,-1.5,10.5));
NumStripCluVsPlaneAll3H1_CuttedMinus=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H1_CuttedMinus","Num Strip Cluster  vs Plane -cluinfo- H1 cls<5",12,-1.5,10.5));
NumStripCluVsPlaneAll3H2_CuttedMinus=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H2_CuttedMinus","Num Strip Cluster  vs Plane -cluinfo- H2 cls<5",12,-1.5,10.5));
NumStripCluVsPlaneAll3H3_CuttedMinus=std::auto_ptr<TProfile> (new TProfile("NumStripCluVsPlaneAll3H3_CuttedMinus","Num Strip Cluster  vs Plane -cluinfo- H3 cls<5",12,-1.5,10.5));







 CumulativeNumPadCluAll3H0= std::auto_ptr<TH1F>(new TH1F("CumulativeNumPadCluAll3H0","H0 Cumulative Pad-Cluster Multiplicity per Plane",1001 ,-0.5 ,1000.5 ));
 CumulativePadCluSizeAll3H0= std::auto_ptr<TH1F>(new TH1F("CumulativePadCluSizeAll3H0","H0 Cumulative Pad-Cluster Size ",501 ,-0.5 ,500.5 ));
 CumulativeNumStripCluAll3H0= std::auto_ptr<TH1F>(new TH1F("CumulativeNumStripCluAll3H0","H0 Cumulative Strip-Cluster Multiplicity per Plane ",1001 ,-0.5 ,1000.5 ));
 CumulativeStripCluSizeAll3H0= std::auto_ptr<TH1F>(new TH1F("CumulativeStripCluSizeAll3H0","H0 Cumulative Strip-Cluster Size ",501 ,-0.5 ,500.5 ));


CumulativePadCluSize_UsedInEffi_HX= std::auto_ptr<TH1F>(new TH1F("CumulativePadCluSize_UsedInEffi_HX","Hx Cumulative Pad-Cluster Multiplicity per Plane Used In Effi",1001 ,-0.5 ,1000.5 ));
CumulativeStripCluSize_UsedInEffi_HX= std::auto_ptr<TH1F>(new TH1F("CumulativeStripCluSize_UsedInEffi_HX","Hx Cumulative Strip-Cluster Size Used in Effi",501 ,-0.5 ,500.5 ));


																	     

PhiNoiseCl1Hit= std::auto_ptr<TH1F>(new TH1F("PhiNoiseCl1Hit","PhiNoiseCl1Hit ",360 ,0. ,360. ));
PhiNoiseCl1Hit->SetDirectory(0);
CLSPadNoiseCl1Hit= std::auto_ptr<TH1F>(new TH1F("CLSPadNoiseCl1Hit","Class 1 Hit Pad Cluster Size",21 ,-0.5 ,20.5 ));
CLSPadNoiseCl1Hit->SetXTitle("Pad cl-size");
CLSPadNoiseCl1Hit->SetDirectory(0);

CLSStripNoiseCl1Hit= std::auto_ptr<TH1F>(new TH1F("CLSStripNoiseCl1Hit","Class 1 Hit Strip Cluster Size",21 ,-0.5 ,20.5 ));
CLSStripNoiseCl1Hit->SetXTitle("Strip cl-size");
CLSStripNoiseCl1Hit->SetDirectory(0);
CLSPhiNoisePad= std::auto_ptr<TH1F>(new TH1F("CLSPhiNoisePad","CLSPhiNoisePad ",21 , -0.5,20.5 ));
CLSPhiNoisePad->SetDirectory(0);
PhiNoisePad= std::auto_ptr<TH1F>(new TH1F("PhiNoisePad","PhiNoisePad ",360 ,0. ,360. ));
PhiNoisePad->SetDirectory(0);
RNoiseStrip= std::auto_ptr<TH1F>(new TH1F("RNoiseStrip","RNoiseStrip ",220,40,150));
RNoiseStrip->SetDirectory(0);
RNoiseCl1Hit= std::auto_ptr<TH1F>(new TH1F("RNoiseCl1Hit"," Noise Hit-R (CL-1 Hit)",25,40,150));
RNoiseCl1Hit->SetXTitle("Hit R (mm)");
RNoiseCl1Hit->SetDirectory(0);
CLSRNoiseStrip= std::auto_ptr<TH1F>(new TH1F("CLSRNoiseStrip"," CLSRNoiseStrip",21 ,-0.5 ,20.5  ));
CLSRNoiseStrip->SetDirectory(0);

 //->Fill((*ithit).GetHitR(),(*ithit).GetHitNumStrip());//TProfile
Cl1Hit_ClsStripVsR= std::auto_ptr<TProfile>(new TProfile("Cl1Hit_ClsStripVsR","Cl1Hit_ClsStripVsR ",22,40.,150.));                     
Cl1Hit_ClsStripVsR->SetDirectory(0);

C11Hit_ClsPadVsR= std::auto_ptr<TProfile>(new TProfile("Cl1Hit_ClsPadVsR","Cl1Hit_ClsPadVsR ",22,40.,150.));//TProfile
C11Hit_ClsPadVsR->SetDirectory(0);

C1Hit_DphiVsPadCLS= std::auto_ptr<TH2F>(new TH2F("C1Hit_DphiVsPadCLS","C1Hit_DphiVsPadCLS",40,0.,40.,30,-0.5,29.5)); 
C1Hit_DphiVsPadCLS->SetDirectory(0);
//C1Hit_DphiVsPadCLS->SetOption("lego");


 C11Hit_ClsPadVsClsStrip= std::auto_ptr<TH2F>(new TH2F("ClsPadVsClsStrip","Cluster Pad Vs Cluster Strip (C1-Hit)",20,-0.5,19.5,30,-0.5,29.5)); 
C11Hit_ClsPadVsClsStrip->SetDirectory(0);
C11Hit_ClsPadVsClsStrip->SetXTitle("Pad CLS");
C11Hit_ClsPadVsClsStrip->SetYTitle("Strip CLS");


C11Hit_ClsPadVsClsPadCol= std::auto_ptr<TH2F>(new TH2F("C11Hit_ClsPadVsClsPadCol","Cluster Pad Vs Cluster Pad #Col (C1-Hit)",20,-0.5,19.5,20,-0.5,19.5)); 
C11Hit_ClsPadVsClsPadCol->SetDirectory(0);
C11Hit_ClsPadVsClsPadCol->SetXTitle("Pad CLS");
C11Hit_ClsPadVsClsPadCol->SetYTitle("Pad CLS Col");









C1Hit_PadClsVsPadNumCol= std::auto_ptr<TH2F>(new TH2F("C1Hit_PadClsVsPadNumCol","C1Hit_PadClsVsPadNumCol",20,-0.5,19.5,15,-0.5,14.5));
C1Hit_PadClsVsPadNumCol->SetDirectory(0);
C1Hit_PadClsVsPadNumRow= std::auto_ptr<TH2F>(new TH2F("C1Hit_PadClsVsPadNumRow","C1Hit_PadClsVsPadNumRow",20,-0.5,19.5,15,-0.5,14.5));
C1Hit_PadClsVsPadNumRow->SetDirectory(0);
C1Hit_Pad_PadNumRowVsPadNumCol= std::auto_ptr<TH2F>(new TH2F("C1Hit_Pad_PadNumRowVsPadNumCol","C1Hit_Pad_PadNumRowVsPadNumCol",20,-0.5,19.5,15,-0.5,14.5));
C1Hit_Pad_PadNumRowVsPadNumCol->SetDirectory(0);

//C11Hit_ClsPadVsClsStrip->SetOption("lego");
C1Hit_Dphi= std::auto_ptr<TH1F>(new TH1F("C1Hit_Dphi","C1Hit_Dphi",15,0,30.));
C1Hit_Dphi->SetDirectory(0);

EFFClusterVsId= std::auto_ptr<TProfile>(new TProfile("EFFClusterVsId","Cluster Efficiency vs Detector-Id",12,-1.5,10.5));
EFFClusterVsId->SetDirectory(0);

 EFFstripVsIdNoCut= std::auto_ptr<TProfile>(new TProfile("EFFstripVsIdNoCut","Strip Efficiency vs Detector-Id (no cut)",12,-1.5,10.5));
EFFstripVsIdNoCut->SetDirectory(0);
 EFFstripVsId= std::auto_ptr<TProfile>(new TProfile("EFFstripVsId","Strip Efficiency vs Detector-Id",12,-1.5,10.5));
EFFstripVsId->SetDirectory(0);
 EFFpadVsIdNoCut= std::auto_ptr<TProfile>(new TProfile("EFFpadVsIdNoCut","Pad Efficiency vs Detector-Id (no cut)",12,-1.5,10.5));
EFFpadVsIdNoCut->SetDirectory(0);
 EFFpadVsId= std::auto_ptr<TProfile>(new TProfile("EFFpadVsId","Pad Efficiency vs Detector-Id",12,-1.5,10.5));
EFFpadVsId->SetDirectory(0);
 EFFHitVsIdNoCut= std::auto_ptr<TProfile>(new TProfile("EFFHitVsIdNoCut","Hit Efficiency vs Detector-Id (no cut)",12,-1.5,10.5));
EFFHitVsIdNoCut->SetDirectory(0);
 EFFHitVsId= std::auto_ptr<TProfile>(new TProfile("EFFHitVsId","Hit Efficiency vs Detector-Id ",12,-1.5,10.5));
EFFHitVsId->SetDirectory(0);


 EFFHitVsIdNormalized= std::auto_ptr<TProfile>(new TProfile("EFFHitVsIdNormalized","Hit Efficiency vs Detector-Id Normalized",12,-1.5,10.5));
 EFFHitVsIdNormalized->SetDirectory(0);

 EFFstripVsIdNormalized= std::auto_ptr<TProfile>(new TProfile("EFFstripVsIdNormalized","Strip Efficiency vs Detector-Id Normalized",12,-1.5,10.5));
 EFFstripVsIdNormalized->SetDirectory(0);  

 EFFpadVsIdNormalized= std::auto_ptr<TProfile>(new TProfile("EFFpadVsIdNormalized","Pad Efficiency vs Detector-Id Normalized",12,-1.5,10.5));
 EFFpadVsIdNormalized->SetDirectory(0);



 


 DRHit12= std::auto_ptr<TProfile>(new TProfile("DRHit12","First two hits radial separation vs DetId",12,-1.5,10.5));
 DRHit12->SetDirectory(0);


 EFF_DRstripNoCutdet0= std::auto_ptr<TH1F>(new TH1F("EFF_DRstripNoCutdet0","Strip Dr-shift (no cut) respect to the expected position (Det 0)",100,-10,10));  
EFF_DRstripNoCutdet0->SetDirectory(0);
 EFF_DRstripdet0= std::auto_ptr<TH1F>(new TH1F(" EFF_DRstripdet0"," ",100,-10,10)); 
EFF_DRstripdet0->SetDirectory(0);
 EFF_DRpadNoCutdet0= std::auto_ptr<TH1F>(new TH1F(" EFF_DRpadNoCutdet0"," ",100,-10,10)); 
EFF_DRpadNoCutdet0->SetDirectory(0);
 EFF_DRpaddet0= std::auto_ptr<TH1F>(new TH1F("EFF_DRpaddet0"," ",100,-10,10)); 
EFF_DRpaddet0->SetDirectory(0);
 EFF_DRHitNoCutdet0 = std::auto_ptr<TH1F>(new TH1F("EFF_DRHitNoCutdet0","Hit Dr-shift (no cut) respect to the expected position (Det 0) ",100,-10,10)); 
EFF_DRHitNoCutdet0->SetDirectory(0);
 EFF_DRHitdet0= std::auto_ptr<TH1F>(new TH1F("EFF_DRHitdet0"," ",100,-10,10)); 
EFF_DRHitdet0->SetDirectory(0);

 EFF_DRstripNoCutdet2= std::auto_ptr<TH1F>(new TH1F("EFF_DRstripNoCutdet2","Strip Dr-shift (no cut) respect to the expected position (Det 2) ",100,-10,10)); 
EFF_DRstripNoCutdet2->SetDirectory(0);
 EFF_DRstripdet2 = std::auto_ptr<TH1F>(new TH1F("EFF_DRstripdet2","Strip Dr-shift respect to the expected position (Det 2) ",100,-10,10)); 
EFF_DRstripdet2->SetDirectory(0);
 EFF_DRpadNoCutdet2= std::auto_ptr<TH1F>(new TH1F("EFF_DRpadNoCutdet2"," ",100,-10,10)); 
EFF_DRpadNoCutdet2->SetDirectory(0);
 EFF_DRpaddet2= std::auto_ptr<TH1F>(new TH1F("EFF_DRpaddet2"," ",100,-10,10)); 
EFF_DRpaddet2->SetDirectory(0);
 EFF_DRHitNoCutdet2= std::auto_ptr<TH1F>(new TH1F("EFF_DRHitNoCutdet2","Hit Dr-shift (no cut) respect to the expected position (Det 2) ",100,-10,10)); 
EFF_DRHitNoCutdet2->SetDirectory(0);
 EFF_DRHitdet2 = std::auto_ptr<TH1F>(new TH1F("EFF_DRHitdet2","Hit Dr-shift respect to the expected position (Det 2) ",80,-4,4)); 
EFF_DRHitdet2->SetDirectory(0);

 EFF_DRstripNoCutdet4= std::auto_ptr<TH1F>(new TH1F(" EFF_DRstripNoCutdet4"," ",100,-10,10)); 
EFF_DRstripNoCutdet4->SetDirectory(0);
 EFF_DRstripdet4  = std::auto_ptr<TH1F>(new TH1F("EFF_DRstripdet4"," ",100,-10,10)); 
EFF_DRstripdet4->SetDirectory(0);
 EFF_DRpadNoCutdet4 = std::auto_ptr<TH1F>(new TH1F("EFF_DRpadNoCutdet4"," ",100,-10,10));
EFF_DRpadNoCutdet4->SetDirectory(0);
 EFF_DRpaddet4= std::auto_ptr<TH1F>(new TH1F("EFF_DRpaddet4"," ",100,-10,10)); 
EFF_DRpaddet4->SetDirectory(0);
 EFF_DRHitNoCutdet4= std::auto_ptr<TH1F>(new TH1F("EFF_DRHitNoCutdet4"," ",100,-10,10)); 
EFF_DRHitNoCutdet4->SetDirectory(0);
 EFF_DRHitdet4= std::auto_ptr<TH1F>(new TH1F("EFF_DRHitdet4"," ",100,-10,10)); 
EFF_DRHitdet4->SetDirectory(0);


 EFF_DRstripNoCutdet6= std::auto_ptr<TH1F>(new TH1F("EFF_DRstripNoCutdet6"," ",100,-10,10)); 
EFF_DRstripNoCutdet6->SetDirectory(0);
 EFF_DRstripdet6= std::auto_ptr<TH1F>(new TH1F("EFF_DRstripdet6"," ",100,-10,10)); 
EFF_DRstripdet6->SetDirectory(0);
 EFF_DRpadNoCutdet6= std::auto_ptr<TH1F>(new TH1F("EFF_DRpadNoCutdet6"," ",100,-10,10)); 
EFF_DRpadNoCutdet6->SetDirectory(0);
 EFF_DRpaddet6= std::auto_ptr<TH1F>(new TH1F("EFF_DRpaddet6"," ",100,-10,10)); 
EFF_DRpaddet6->SetDirectory(0);
 EFF_DRHitNoCutdet6= std::auto_ptr<TH1F>(new TH1F("EFF_DRHitNoCutdet6"," ",100,-10,10)); 
EFF_DRHitNoCutdet6->SetDirectory(0);
 EFF_DRHitdet6= std::auto_ptr<TH1F>(new TH1F("EFF_DRHitdet6"," ",100,-10,10)); 
EFF_DRHitdet6->SetDirectory(0);


 EFF_DRstripNoCutdet8= std::auto_ptr<TH1F>(new TH1F("EFF_DRstripNoCutdet8"," ",100,-10,10)); 
EFF_DRstripNoCutdet8->SetDirectory(0);
 EFF_DRstripdet8= std::auto_ptr<TH1F>(new TH1F("EFF_DRstripdet8"," ",100,-10,10)); 
EFF_DRstripdet8->SetDirectory(0);
 EFF_DRpadNoCutdet8= std::auto_ptr<TH1F>(new TH1F("EFF_DRpadNoCutdet8"," ",100,-10,10)); 
EFF_DRpadNoCutdet8->SetDirectory(0);
 EFF_DRpaddet8= std::auto_ptr<TH1F>(new TH1F("EFF_DRpaddet8"," ",100,-10,10)); 
EFF_DRpaddet8->SetDirectory(0);
 EFF_DRHitNoCutdet8= std::auto_ptr<TH1F>(new TH1F("EFF_DRHitNoCutdet8"," ",100,-10,10)); 
EFF_DRHitNoCutdet8->SetDirectory(0);
 EFF_DRHitdet8= std::auto_ptr<TH1F>(new TH1F("EFF_DRHitdet8"," ",100,-10,10)); 
EFF_DRHitdet8->SetDirectory(0);




 NoiseEventsStat= std::auto_ptr<TH1F>(new TH1F("NoiseEventsStat"," ",5,0,4));
NoiseEventsStat->SetDirectory(0);




 EFF_DRstripNoCutdet1= std::auto_ptr<TH1F>(new TH1F("EFF_DRstripNoCutdet1"," ",100,-10,10));  
EFF_DRstripNoCutdet1->SetDirectory(0);
 EFF_DRstripdet1= std::auto_ptr<TH1F>(new TH1F(" EFF_DRstripdet1"," ",100,-10,10)); 
EFF_DRstripdet1->SetDirectory(0);
 EFF_DRpadNoCutdet1= std::auto_ptr<TH1F>(new TH1F(" EFF_DRpadNoCutdet1"," ",100,-10,10)); 
EFF_DRpadNoCutdet1->SetDirectory(0);
 EFF_DRpaddet1= std::auto_ptr<TH1F>(new TH1F("EFF_DRpaddet1"," ",100,-10,10)); 
EFF_DRpaddet1->SetDirectory(0);
 EFF_DRHitNoCutdet1 = std::auto_ptr<TH1F>(new TH1F("EFF_DRHitNoCutdet1"," ",100,-10,10)); 
EFF_DRHitNoCutdet1->SetDirectory(0);
 EFF_DRHitdet1= std::auto_ptr<TH1F>(new TH1F("EFF_DRHitdet1"," ",100,-10,10)); 
EFF_DRHitdet1->SetDirectory(0);

 EFF_DRstripNoCutdet3= std::auto_ptr<TH1F>(new TH1F("EFF_DRstripNoCutdet3"," ",100,-10,10)); 
EFF_DRstripNoCutdet3->SetDirectory(0);
 EFF_DRstripdet3 = std::auto_ptr<TH1F>(new TH1F("EFF_DRstripdet3"," ",100,-10,10)); 
EFF_DRstripdet3->SetDirectory(0);
 EFF_DRpadNoCutdet3= std::auto_ptr<TH1F>(new TH1F("EFF_DRpadNoCutdet3"," ",100,-10,10)); 
EFF_DRpadNoCutdet3->SetDirectory(0);
 EFF_DRpaddet3= std::auto_ptr<TH1F>(new TH1F("EFF_DRpaddet3"," ",100,-10,10)); 
EFF_DRpaddet3->SetDirectory(0);
 EFF_DRHitNoCutdet3= std::auto_ptr<TH1F>(new TH1F("EFF_DRHitNoCutdet3"," ",100,-10,10)); 
EFF_DRHitNoCutdet3->SetDirectory(0);
 EFF_DRHitdet3 = std::auto_ptr<TH1F>(new TH1F("EFF_DRHitdet3"," ",100,-10,10)); 
EFF_DRHitdet3->SetDirectory(0);

 EFF_DRstripNoCutdet5= std::auto_ptr<TH1F>(new TH1F(" EFF_DRstripNoCutdet5"," ",100,-10,10)); 
EFF_DRstripNoCutdet5->SetDirectory(0);
 EFF_DRstripdet5  = std::auto_ptr<TH1F>(new TH1F("EFF_DRstripdet5"," ",100,-10,10)); 
EFF_DRstripdet5->SetDirectory(0);
 EFF_DRpadNoCutdet5 = std::auto_ptr<TH1F>(new TH1F("EFF_DRpadNoCutdet5"," ",100,-10,10));
EFF_DRpadNoCutdet5->SetDirectory(0);
 EFF_DRpaddet5= std::auto_ptr<TH1F>(new TH1F("EFF_DRpaddet5"," ",100,-10,10)); 
EFF_DRpaddet5->SetDirectory(0);
 EFF_DRHitNoCutdet5= std::auto_ptr<TH1F>(new TH1F("EFF_DRHitNoCutdet5"," ",100,-10,10)); 
EFF_DRHitNoCutdet5->SetDirectory(0);
 EFF_DRHitdet5= std::auto_ptr<TH1F>(new TH1F("EFF_DRHitdet5"," ",100,-10,10)); 
EFF_DRHitdet5->SetDirectory(0);


 EFF_DRstripNoCutdet7= std::auto_ptr<TH1F>(new TH1F("EFF_DRstripNoCutdet7"," ",100,-10,10)); 
EFF_DRstripNoCutdet7->SetDirectory(0);
 EFF_DRstripdet7= std::auto_ptr<TH1F>(new TH1F("EFF_DRstripdet7"," ",100,-10,10)); 
EFF_DRstripdet7->SetDirectory(0);
 EFF_DRpadNoCutdet7= std::auto_ptr<TH1F>(new TH1F("EFF_DRpadNoCutdet7"," ",100,-10,10)); 
EFF_DRpadNoCutdet7->SetDirectory(0);
 EFF_DRpaddet7= std::auto_ptr<TH1F>(new TH1F("EFF_DRpaddet7"," ",100,-10,10)); 
EFF_DRpaddet7->SetDirectory(0);
 EFF_DRHitNoCutdet7= std::auto_ptr<TH1F>(new TH1F("EFF_DRHitNoCutdet7"," ",100,-10,10)); 
EFF_DRHitNoCutdet7->SetDirectory(0);
 EFF_DRHitdet7= std::auto_ptr<TH1F>(new TH1F("EFF_DRHitdet7"," ",100,-10,10)); 
EFF_DRHitdet7->SetDirectory(0);


 EFF_DRstripNoCutdet9= std::auto_ptr<TH1F>(new TH1F("EFF_DRstripNoCutdet9"," ",100,-10,10)); 
EFF_DRstripNoCutdet9->SetDirectory(0);
 EFF_DRstripdet9= std::auto_ptr<TH1F>(new TH1F("EFF_DRstripdet9"," ",100,-10,10)); 
EFF_DRstripdet9->SetDirectory(0);
 EFF_DRpadNoCutdet9= std::auto_ptr<TH1F>(new TH1F("EFF_DRpadNoCutdet9"," ",100,-10,10)); 
EFF_DRpadNoCutdet9->SetDirectory(0);
 EFF_DRpaddet9= std::auto_ptr<TH1F>(new TH1F("EFF_DRpaddet9"," ",100,-10,10)); 
EFF_DRpaddet9->SetDirectory(0);
 EFF_DRHitNoCutdet9= std::auto_ptr<TH1F>(new TH1F("EFF_DRHitNoCutdet9"," ",100,-10,10)); 
EFF_DRHitNoCutdet9->SetDirectory(0);
 EFF_DRHitdet9= std::auto_ptr<TH1F>(new TH1F("EFF_DRHitdet9"," ",100,-10,10)); 
EFF_DRHitdet9->SetDirectory(0);






 EFF_DPhistripNoCutdet0= std::auto_ptr<TH1F>(new TH1F("EFF_DPhistripNoCutdet0"," ",100,-10,10));  
EFF_DPhistripNoCutdet0->SetDirectory(0);
 EFF_DPhistripdet0= std::auto_ptr<TH1F>(new TH1F(" EFF_DPhistripdet0"," ",100,-10,10)); 
EFF_DPhistripdet0->SetDirectory(0);
 EFF_DPhipadNoCutdet0= std::auto_ptr<TH1F>(new TH1F(" EFF_DPhipadNoCutdet0"," ",100,-10,10)); 
EFF_DPhipadNoCutdet0->SetDirectory(0);
 EFF_DPhipaddet0= std::auto_ptr<TH1F>(new TH1F("EFF_DPhipaddet0"," ",100,-10,10)); 
EFF_DPhipaddet0->SetDirectory(0);
 EFF_DPhiHitNoCutdet0 = std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitNoCutdet0"," ",100,-10,10)); 
EFF_DPhiHitNoCutdet0->SetDirectory(0);
 EFF_DPhiHitdet0= std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitdet0"," ",100,-10,10)); 
EFF_DPhiHitdet0->SetDirectory(0);

 EFF_DPhistripNoCutdet2= std::auto_ptr<TH1F>(new TH1F("EFF_DPhistripNoCutdet2"," ",100,-10,10)); 
EFF_DPhistripNoCutdet2->SetDirectory(0);
 EFF_DPhistripdet2 = std::auto_ptr<TH1F>(new TH1F("EFF_DPhistripdet2"," ",100,-10,10)); 
EFF_DPhistripdet2->SetDirectory(0);
 EFF_DPhipadNoCutdet2= std::auto_ptr<TH1F>(new TH1F("EFF_DPhipadNoCutdet2"," ",100,-10,10)); 
EFF_DPhipadNoCutdet2->SetDirectory(0);
 EFF_DPhipaddet2= std::auto_ptr<TH1F>(new TH1F("EFF_DPhipaddet2"," ",100,-10,10)); 
EFF_DPhipaddet2->SetDirectory(0);
 EFF_DPhiHitNoCutdet2 = std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitNoCutdet2","Hit D#phi-shift (no cut) respect to the expected position (Det 2) ",100,-10,10)); 
EFF_DPhiHitNoCutdet2->SetDirectory(0);
 EFF_DPhiHitdet2 = std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitdet2","Hit D#phi-shift respect to the expected position (Det 2) ",20,-4,4)); 
EFF_DPhiHitdet2->SetDirectory(0);


 EFF_DPhistripNoCutdet4= std::auto_ptr<TH1F>(new TH1F(" EFF_DPhistripNoCutdet4"," ",100,-10,10));
EFF_DPhistripNoCutdet4->SetDirectory(0);
 EFF_DPhistripdet4  = std::auto_ptr<TH1F>(new TH1F("EFF_DPhistripdet4"," ",100,-10,10)); 
EFF_DPhistripdet4->SetDirectory(0);
 EFF_DPhipadNoCutdet4 = std::auto_ptr<TH1F>(new TH1F("EFF_DPhipadNoCutdet4"," ",100,-10,10)); 
EFF_DPhipadNoCutdet4->SetDirectory(0);
 EFF_DPhipaddet4= std::auto_ptr<TH1F>(new TH1F("EFF_DPhipaddet4"," ",100,-10,10)); 
EFF_DPhipaddet4->SetDirectory(0);
 EFF_DPhiHitNoCutdet4= std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitNoCutdet4"," ",100,-10,10)); 
EFF_DPhiHitNoCutdet4->SetDirectory(0);
 EFF_DPhiHitdet4= std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitdet4"," ",100,-10,10)); 
EFF_DPhiHitdet4->SetDirectory(0);

 EFF_DPhistripNoCutdet6= std::auto_ptr<TH1F>(new TH1F("EFF_DPhistripNoCutdet6"," ",100,-10,10)); 
EFF_DPhistripNoCutdet6->SetDirectory(0);
 EFF_DPhistripdet6= std::auto_ptr<TH1F>(new TH1F("EFF_DPhistripdet6"," ",100,-10,10)); 
EFF_DPhistripdet6->SetDirectory(0);
 EFF_DPhipadNoCutdet6= std::auto_ptr<TH1F>(new TH1F("EFF_DPhipadNoCutdet6"," ",100,-10,10)); 
EFF_DPhipadNoCutdet6->SetDirectory(0);
 EFF_DPhipaddet6= std::auto_ptr<TH1F>(new TH1F("EFF_DPhipaddet6"," ",100,-10,10)); 
EFF_DPhipaddet6->SetDirectory(0);
 EFF_DPhiHitNoCutdet6= std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitNoCutdet6"," ",100,-10,10)); 
EFF_DPhiHitNoCutdet6->SetDirectory(0);
 EFF_DPhiHitdet6= std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitdet6"," ",100,-10,10)); 
EFF_DPhiHitdet6->SetDirectory(0);

 EFF_DPhistripNoCutdet8= std::auto_ptr<TH1F>(new TH1F("EFF_DPhistripNoCutdet8"," ",100,-10,10));
EFF_DPhistripNoCutdet8->SetDirectory(0);
 EFF_DPhistripdet8= std::auto_ptr<TH1F>(new TH1F("EFF_DPhistripdet8"," ",100,-10,10)); 
EFF_DPhistripdet8->SetDirectory(0);
 EFF_DPhipadNoCutdet8= std::auto_ptr<TH1F>(new TH1F("EFF_DPhipadNoCutdet8"," ",100,-10,10)); 
EFF_DPhipadNoCutdet8->SetDirectory(0);
 EFF_DPhipaddet8= std::auto_ptr<TH1F>(new TH1F("EFF_DPhipaddet8"," ",100,-10,10)); 
EFF_DPhipaddet8->SetDirectory(0);
 EFF_DPhiHitNoCutdet8= std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitNoCutdet8"," ",100,-10,10)); 
EFF_DPhiHitNoCutdet8->SetDirectory(0);
 EFF_DPhiHitdet8= std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitdet8"," ",100,-10,10));
EFF_DPhiHitdet8->SetDirectory(0);




EFF_DPhistripNoCutdet1= std::auto_ptr<TH1F>(new TH1F("EFF_DPhistripNoCutdet1"," ",100,-10,10));  
EFF_DPhistripNoCutdet1->SetDirectory(0);
 EFF_DPhistripdet1= std::auto_ptr<TH1F>(new TH1F(" EFF_DPhistripdet1"," ",100,-10,10)); 
EFF_DPhistripdet1->SetDirectory(0);
 EFF_DPhipadNoCutdet1= std::auto_ptr<TH1F>(new TH1F(" EFF_DPhipadNoCutdet1"," ",100,-10,10)); 
EFF_DPhipadNoCutdet1->SetDirectory(0);
 EFF_DPhipaddet1= std::auto_ptr<TH1F>(new TH1F("EFF_DPhipaddet1"," ",100,-10,10)); 
EFF_DPhipaddet1->SetDirectory(0);
 EFF_DPhiHitNoCutdet1 = std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitNoCutdet1"," ",100,-10,10)); 
EFF_DPhiHitNoCutdet1->SetDirectory(0);
 EFF_DPhiHitdet1= std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitdet1"," ",100,-10,10)); 
EFF_DPhiHitdet1->SetDirectory(0);

 EFF_DPhistripNoCutdet3= std::auto_ptr<TH1F>(new TH1F("EFF_DPhistripNoCutdet3"," ",100,-10,10)); 
EFF_DPhistripNoCutdet3->SetDirectory(0);
 EFF_DPhistripdet3 = std::auto_ptr<TH1F>(new TH1F("EFF_DPhistripdet3"," ",100,-10,10)); 
EFF_DPhistripdet3->SetDirectory(0);
 EFF_DPhipadNoCutdet3= std::auto_ptr<TH1F>(new TH1F("EFF_DPhipadNoCutdet3"," ",100,-10,10)); 
EFF_DPhipadNoCutdet3->SetDirectory(0);
 EFF_DPhipaddet3= std::auto_ptr<TH1F>(new TH1F("EFF_DPhipaddet3"," ",100,-10,10)); 
EFF_DPhipaddet3->SetDirectory(0);
 EFF_DPhiHitNoCutdet3 = std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitNoCutdet3"," ",100,-10,10)); 
EFF_DPhiHitNoCutdet3->SetDirectory(0);
 EFF_DPhiHitdet3 = std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitdet3"," ",100,-10,10)); 
EFF_DPhiHitdet3->SetDirectory(0);


 EFF_DPhistripNoCutdet5= std::auto_ptr<TH1F>(new TH1F(" EFF_DPhistripNoCutdet5"," ",100,-10,10));
EFF_DPhistripNoCutdet5->SetDirectory(0);
 EFF_DPhistripdet5  = std::auto_ptr<TH1F>(new TH1F("EFF_DPhistripdet5"," ",100,-10,10)); 
EFF_DPhistripdet5->SetDirectory(0);
 EFF_DPhipadNoCutdet5 = std::auto_ptr<TH1F>(new TH1F("EFF_DPhipadNoCutdet5"," ",100,-10,10)); 
EFF_DPhipadNoCutdet5->SetDirectory(0);
 EFF_DPhipaddet5= std::auto_ptr<TH1F>(new TH1F("EFF_DPhipaddet5"," ",100,-10,10)); 
EFF_DPhipaddet5->SetDirectory(0);
 EFF_DPhiHitNoCutdet5= std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitNoCutdet5"," ",100,-10,10)); 
EFF_DPhiHitNoCutdet5->SetDirectory(0);
 EFF_DPhiHitdet5= std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitdet5"," ",100,-10,10)); 
EFF_DPhiHitdet5->SetDirectory(0);

 EFF_DPhistripNoCutdet7= std::auto_ptr<TH1F>(new TH1F("EFF_DPhistripNoCutdet7"," ",100,-10,10)); 
EFF_DPhistripNoCutdet7->SetDirectory(0);
 EFF_DPhistripdet7= std::auto_ptr<TH1F>(new TH1F("EFF_DPhistripdet7"," ",100,-10,10)); 
EFF_DPhistripdet7->SetDirectory(0);
 EFF_DPhipadNoCutdet7= std::auto_ptr<TH1F>(new TH1F("EFF_DPhipadNoCutdet7"," ",100,-10,10)); 
EFF_DPhipadNoCutdet7->SetDirectory(0);
 EFF_DPhipaddet7= std::auto_ptr<TH1F>(new TH1F("EFF_DPhipaddet7"," ",100,-10,10)); 
EFF_DPhipaddet7->SetDirectory(0);
 EFF_DPhiHitNoCutdet7= std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitNoCutdet7"," ",100,-10,10)); 
EFF_DPhiHitNoCutdet7->SetDirectory(0);
 EFF_DPhiHitdet7= std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitdet7"," ",100,-10,10)); 
EFF_DPhiHitdet7->SetDirectory(0);

 EFF_DPhistripNoCutdet9= std::auto_ptr<TH1F>(new TH1F("EFF_DPhistripNoCutdet9"," ",100,-10,10));
EFF_DPhistripNoCutdet9->SetDirectory(0);
 EFF_DPhistripdet9= std::auto_ptr<TH1F>(new TH1F("EFF_DPhistripdet9"," ",100,-10,10)); 
EFF_DPhistripdet9->SetDirectory(0);
 EFF_DPhipadNoCutdet9= std::auto_ptr<TH1F>(new TH1F("EFF_DPhipadNoCutdet9"," ",100,-10,10)); 
EFF_DPhipadNoCutdet9->SetDirectory(0);
 EFF_DPhipaddet9= std::auto_ptr<TH1F>(new TH1F("EFF_DPhipaddet9"," ",100,-10,10)); 
EFF_DPhipaddet9->SetDirectory(0);
 EFF_DPhiHitNoCutdet9= std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitNoCutdet9"," ",100,-10,10)); 
EFF_DPhiHitNoCutdet9->SetDirectory(0);
 EFF_DPhiHitdet9= std::auto_ptr<TH1F>(new TH1F("EFF_DPhiHitdet9"," ",100,-10,10));
EFF_DPhiHitdet9->SetDirectory(0);





AL_DRHitNoCutdet0= std::auto_ptr<TH1F>(new TH1F("AL_DRHitNoCutdet0"," ",40,-4,4));
AL_DRHitNoCutdet0->SetDirectory(0);

AL_DPhiHitNoCutdet0= std::auto_ptr<TH1F>(new TH1F("AL_DPhiHitNoCutdet0"," ",40,-4,4));
AL_DPhiHitNoCutdet0->SetDirectory(0);


AL_DRHitNoCutdet1= std::auto_ptr<TH1F>(new TH1F("AL_DRHitNoCutdet1"," ",40,-4,4));
AL_DRHitNoCutdet1->SetDirectory(0);



AL_DRHitNoCutdet2= std::auto_ptr<TH1F>(new TH1F("AL_DRHitNoCutdet2"," ",40,-4,4));
AL_DRHitNoCutdet2->SetDirectory(0);



AL_DRHitNoCutdet3= std::auto_ptr<TH1F>(new TH1F("AL_DRHitNoCutdet3"," ",40,-4,4));
AL_DRHitNoCutdet3->SetDirectory(0);

AL_DRHitNoCutdet4= std::auto_ptr<TH1F>(new TH1F("AL_DRHitNoCutdet4"," ",40,-4,4));
AL_DRHitNoCutdet4->SetDirectory(0);

AL_DRHitNoCutdet5= std::auto_ptr<TH1F>(new TH1F("AL_DRHitNoCutdet5"," ",40,-4,4));
AL_DRHitNoCutdet5->SetDirectory(0);

AL_DRHitNoCutdet6= std::auto_ptr<TH1F>(new TH1F("AL_DRHitNoCutdet6"," ",40,-4,4));
AL_DRHitNoCutdet6->SetDirectory(0);

AL_DRHitNoCutdet7= std::auto_ptr<TH1F>(new TH1F("AL_DRHitNoCutdet7"," ",40,-4,4));
AL_DRHitNoCutdet7->SetDirectory(0);

AL_DRHitNoCutdet8= std::auto_ptr<TH1F>(new TH1F("AL_DRHitNoCutdet8"," ",40,-4,4));
AL_DRHitNoCutdet8->SetDirectory(0);

AL_DRHitNoCutdet9= std::auto_ptr<TH1F>(new TH1F("AL_DRHitNoCutdet9"," ",40,-4,4));
AL_DRHitNoCutdet9->SetDirectory(0);



AL_DPhiHitNoCutdet1= std::auto_ptr<TH1F>(new TH1F("AL_DPhiHitNoCutdet1"," ",40,-4,4));
AL_DPhiHitNoCutdet1->SetDirectory(0);

AL_DPhiHitNoCutdet2= std::auto_ptr<TH1F>(new TH1F("AL_DPhiHitNoCutdet2"," ",40,-4,4));
AL_DPhiHitNoCutdet2->SetDirectory(0);

AL_DPhiHitNoCutdet3= std::auto_ptr<TH1F>(new TH1F("AL_DPhiHitNoCutdet3"," ",40,-4,4));
AL_DPhiHitNoCutdet3->SetDirectory(0);

AL_DPhiHitNoCutdet4= std::auto_ptr<TH1F>(new TH1F("AL_DPhiHitNoCutdet4"," ",40,-4,4));
AL_DPhiHitNoCutdet4->SetDirectory(0);

AL_DPhiHitNoCutdet5= std::auto_ptr<TH1F>(new TH1F("AL_DPhiHitNoCutdet5"," ",40,-4,4));
AL_DPhiHitNoCutdet5->SetDirectory(0);

AL_DPhiHitNoCutdet6= std::auto_ptr<TH1F>(new TH1F("AL_DPhiHitNoCutdet6"," ",40,-4,4));
AL_DPhiHitNoCutdet6->SetDirectory(0);

AL_DPhiHitNoCutdet7= std::auto_ptr<TH1F>(new TH1F("AL_DPhiHitNoCutdet7"," ",40,-4,4));
AL_DPhiHitNoCutdet7->SetDirectory(0);

AL_DPhiHitNoCutdet8= std::auto_ptr<TH1F>(new TH1F("AL_DPhiHitNoCutdet8"," ",40,-4,4));
AL_DPhiHitNoCutdet8->SetDirectory(0);

AL_DPhiHitNoCutdet9= std::auto_ptr<TH1F>(new TH1F("AL_DPhiHitNoCutdet9"," ",40,-4,4));
AL_DPhiHitNoCutdet9->SetDirectory(0);




 EFF_CLSHitdet0= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet0","CLSHitdet0 ",129,-0.5,129.5));
EFF_CLSHitdet0->SetDirectory(0);



EFF_CLSStripdet0= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet0","Strip Cluster size, det 0 ",129,-0.5,128.5)); 
EFF_CLSStripdet0->SetDirectory(0);

EFF_CLSPaddet0= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet0","CLSPaddet0 ",129,-0.5,128.5)); 
EFF_CLSPaddet0->SetDirectory(0);

 EFF_CLSHitdet2= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet2","CLSHitdet2 ",129,-0.5,128.5));
 EFF_CLSHitdet2->SetDirectory(0);
 EFF_CLSStripdet2= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet2","Strip Cluster size, det 2 ",129,-0.5,128.5));
 EFF_CLSStripdet2->SetDirectory(0);
 EFF_CLSPaddet2= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet2"," CLSPaddet2",129,-0.5,128.5)); 
EFF_CLSPaddet2->SetDirectory(0);

 EFF_CLSHitdet4= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet4","CLSHitdet4 ",129,-0.5,128.5)); 
EFF_CLSHitdet4->SetDirectory(0);
 EFF_CLSStripdet4= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet4","Strip Cluster size, det 4 ",129,-0.5,128.5));
 EFF_CLSStripdet4->SetDirectory(0);
 EFF_CLSPaddet4= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet4","CLSPaddet4 ",129,-0.5,128.5)); 
EFF_CLSPaddet4->SetDirectory(0);

 EFF_CLSHitdet6= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet6"," CLSHitdet6",129,-0.5,128.5)); 
EFF_CLSHitdet6->SetDirectory(0);
 EFF_CLSStripdet6= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet6","Strip Cluster size, det 6 ",129,-0.5,128.5));
 EFF_CLSStripdet6->SetDirectory(0);
 EFF_CLSPaddet6= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet6","CLSPaddet6 ",129,-0.5,128.5)); 
EFF_CLSPaddet6->SetDirectory(0);

 EFF_CLSHitdet8= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet8","CLSHitdet8 ",129,-0.5,128.5)); 
EFF_CLSHitdet8->SetDirectory(0);
 EFF_CLSStripdet8= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet8","Strip Cluster size, det 8 ",129,-0.5,128.5)); 
EFF_CLSStripdet8->SetDirectory(0);
 EFF_CLSPaddet8= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet8","CLSPaddet8 ",129,-0.5,128.5)); 
EFF_CLSPaddet8->SetDirectory(0);




EFF_CLSHitdet1= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet1","CLSHitdet1 ",129,-0.5,128.5));
 EFF_CLSHitdet1->SetDirectory(0);
 EFF_CLSStripdet1= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet1","Strip Cluster size, det 1 ",129,-0.5,128.5)); 
EFF_CLSStripdet1->SetDirectory(0);
 EFF_CLSPaddet1= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet1","CLSPaddet1 ",129,-0.5,128.5)); 
EFF_CLSPaddet1->SetDirectory(0);

 EFF_CLSHitdet3= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet3","CLSHitdet3 ",129,-0.5,128.5));
 EFF_CLSHitdet3->SetDirectory(0);
 EFF_CLSStripdet3= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet3"," CLSStripdet3",129,-0.5,128.5));
 EFF_CLSStripdet3->SetDirectory(0);
 EFF_CLSPaddet3= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet3","CLSPaddet3 ",129,-0.5,128.5)); 
EFF_CLSPaddet3->SetDirectory(0);

 EFF_CLSHitdet5= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet5","CLSHitdet5 ",129,-0.5,128.5)); 
EFF_CLSHitdet5->SetDirectory(0);
 EFF_CLSStripdet5= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet5","Strip Cluster size, det 5 ",129,-0.5,128.5));
 EFF_CLSStripdet5->SetDirectory(0);
 EFF_CLSPaddet5= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet5","CLSPaddet5 ",129,-0.5,128.5)); 
EFF_CLSPaddet5->SetDirectory(0);

 EFF_CLSHitdet7= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet7","CLSHitdet7 ",129,-0.5,128.5)); 
EFF_CLSHitdet7->SetDirectory(0);
 EFF_CLSStripdet7= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet7","Strip Cluster size, det 7 ",129,-0.5,128.5));
 EFF_CLSStripdet7->SetDirectory(0);
 EFF_CLSPaddet7= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet7","CLSPaddet7 ",129,-0.5,128.5)); 
EFF_CLSPaddet7->SetDirectory(0);

 EFF_CLSHitdet9= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet9","CLSHitdet9 ",129,-0.5,128.5)); 
EFF_CLSHitdet9->SetDirectory(0);
 EFF_CLSStripdet9= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet9","Strip Cluster size, det 9 ",129,-0.5,128.5)); 
EFF_CLSStripdet9->SetDirectory(0);
 EFF_CLSPaddet9= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet9"," CLSPaddet9",129,-0.5,128.5)); 
EFF_CLSPaddet9->SetDirectory(0);








 EFF_CLSHitdet0NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet0NoCut","CLSHitdet0NoCut ",129,-0.5,128.5));
 EFF_CLSHitdet0NoCut->SetDirectory(0);
 EFF_CLSStripdet0NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet0NoCut","Strip Cluster size, det 0 NoCut ",129,-0.5,128.5)); 
EFF_CLSStripdet0NoCut->SetDirectory(0);
 EFF_CLSPaddet0NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet0NoCut","Pad Cluster size, det 0 NoCut ",129,-0.5,128.5)); 
EFF_CLSPaddet0NoCut->SetDirectory(0);

 EFF_CLSHitdet2NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet2NoCut","CLSHitdet2NoCut ",129,-0.5,128.5)); 
EFF_CLSHitdet2NoCut->SetDirectory(0);
 EFF_CLSStripdet2NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet2NoCut","Strip Cluster size, det 2 NoCut ",129,-0.5,128.5));
 EFF_CLSStripdet2NoCut->SetDirectory(0);
 EFF_CLSPaddet2NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet2NoCut","Pad Cluster size, det 2 NoCut ",129,-0.5,128.5)); 
EFF_CLSPaddet2NoCut->SetDirectory(0);

 EFF_CLSHitdet4NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet4NoCut","CLSHitdet4NoCut ",129,-0.5,128.5)); 
EFF_CLSHitdet4NoCut->SetDirectory(0);
 EFF_CLSStripdet4NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet4NoCut","Strip Cluster size, det 4 NoCut ",129,-0.5,128.5));
 EFF_CLSStripdet4NoCut->SetDirectory(0);
 EFF_CLSPaddet4NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet4NoCut","Pad Cluster size, det 4 NoCut ",129,-0.5,128.5)); 
EFF_CLSPaddet4NoCut->SetDirectory(0);

 EFF_CLSHitdet6NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet6NoCut","CLSHitdet6NoCut ",129,-0.5,128.5)); 
EFF_CLSHitdet6NoCut->SetDirectory(0);
 EFF_CLSStripdet6NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet6NoCut","Strip Cluster size, det 6 NoCut ",129,-0.5,128.5)); 
EFF_CLSStripdet6NoCut->SetDirectory(0);
 EFF_CLSPaddet6NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet6NoCut","Pad Cluster size, det 6 NoCut ",129,-0.5,128.5)); 
EFF_CLSPaddet6NoCut->SetDirectory(0);

 EFF_CLSHitdet8NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet8NoCut","CLSHitdet8NoCut ",129,-0.5,128.5)); 
EFF_CLSHitdet8NoCut->SetDirectory(0);
 EFF_CLSStripdet8NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet8NoCut","Strip Cluster size, det 8 NoCut ",129,-0.5,128.5));
 EFF_CLSStripdet8NoCut->SetDirectory(0);
 EFF_CLSPaddet8NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet8NoCut","Pad Cluster size, det 8 NoCut ",129,-0.5,128.5)); 
EFF_CLSPaddet8NoCut->SetDirectory(0);


 EFF_CLSHitdet1NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet1NoCut"," ",129,-0.5,128.5));
 EFF_CLSHitdet1NoCut->SetDirectory(0);
 EFF_CLSStripdet1NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet1NoCut"," ",129,-0.5,128.5)); 
EFF_CLSStripdet1NoCut->SetDirectory(0);
 EFF_CLSPaddet1NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet1NoCut"," ",129,-0.5,128.5)); 
EFF_CLSPaddet1NoCut->SetDirectory(0);

 EFF_CLSHitdet3NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet3NoCut"," ",129,-0.5,128.5)); 
EFF_CLSHitdet3NoCut->SetDirectory(0);
 EFF_CLSStripdet3NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet3NoCut","Strip Cluster size, det 3 NoCut ",129,-0.5,128.5));
 EFF_CLSStripdet3NoCut->SetDirectory(0);
 EFF_CLSPaddet3NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet3NoCut","Pad Cluster size, det 3 NoCut ",129,-0.5,128.5)); 
EFF_CLSPaddet3NoCut->SetDirectory(0);

 EFF_CLSHitdet5NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet5NoCut","CLSHitdet5NoCut ",129,-0.5,128.5)); 
EFF_CLSHitdet5NoCut->SetDirectory(0);
 EFF_CLSStripdet5NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet5NoCut","Strip Cluster size, det 5 NoCut ",129,-0.5,128.5));
 EFF_CLSStripdet5NoCut->SetDirectory(0);
 EFF_CLSPaddet5NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet5NoCut","Pad Cluster size, det 5 NoCut ",129,-0.5,128.5)); 
EFF_CLSPaddet5NoCut->SetDirectory(0);

 EFF_CLSHitdet7NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet7NoCut"," CLSHitdet7NoCut",129,-0.5,128.5)); 
EFF_CLSHitdet7NoCut->SetDirectory(0);
 EFF_CLSStripdet7NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet7NoCut","Strip Cluster size, det 7 NoCut ",129,-0.5,128.5)); 
EFF_CLSStripdet7NoCut->SetDirectory(0);
 EFF_CLSPaddet7NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet7NoCut","Pad Cluster size, det 7 NoCut ",129,-0.5,128.5)); 
EFF_CLSPaddet7NoCut->SetDirectory(0);

 EFF_CLSHitdet9NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSHitdet9NoCut","CLSHitdet9NoCut ",129,-0.5,128.5)); 
EFF_CLSHitdet9NoCut->SetDirectory(0);
 EFF_CLSStripdet9NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSStripdet9NoCut","Strip Cluster size, det 9 NoCut ",129,-0.5,128.5));
 EFF_CLSStripdet9NoCut->SetDirectory(0);
 EFF_CLSPaddet9NoCut= std::auto_ptr<TH1F>(new TH1F("EFF_CLSPaddet9NoCut","Pad Cluster size, det 9 NoCut ",129,-0.5,128.5)); 
EFF_CLSPaddet9NoCut->SetDirectory(0);






EFF_C1CLSStripdet0= std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSStripdet0","Class 1 Hit Strip Cluster Size, Det 0",21,-0.5,20.5)) ;
EFF_C1CLSStripdet0->SetXTitle("Strip Cl-size");
EFF_C1CLSStripdet0->SetDirectory(0);

EFF_C1CLSPaddet0=  std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSPaddet0","Class 1 Hit Pad Cluster Size, Det 0",21,-0.5,20.5)) ;
EFF_C1CLSPaddet0->SetXTitle("Pad Cl-size");
EFF_C1CLSPaddet0->SetDirectory(0);

EFF_C1CLSStripdet1= std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSStripdet1"," ",129,-0.5,128.5)) ;
EFF_C1CLSStripdet1->SetDirectory(0);
EFF_C1CLSPaddet1=  std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSPaddet1"," ",129,-0.5,128.5)) ;
EFF_C1CLSPaddet1->SetDirectory(0);

EFF_C1CLSStripdet2= std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSStripdet2"," ",129,-0.5,128.5)) ;
EFF_C1CLSStripdet2->SetDirectory(0);
EFF_C1CLSPaddet2=  std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSPaddet2"," ",129,-0.5,128.5)) ;
EFF_C1CLSPaddet2->SetDirectory(0);

EFF_C1CLSStripdet3= std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSStripdet3"," ",129,-0.5,128.5)) ;
EFF_C1CLSStripdet3->SetDirectory(0);
EFF_C1CLSPaddet3=  std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSPaddet3"," ",129,-0.5,128.5)) ;
EFF_C1CLSPaddet3->SetDirectory(0);

EFF_C1CLSStripdet4= std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSStripdet4"," ",129,-0.5,128.5)) ;
EFF_C1CLSStripdet4->SetDirectory(0);
EFF_C1CLSPaddet4=  std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSPaddet4"," ",129,-0.5,128.5)) ;
EFF_C1CLSPaddet4->SetDirectory(0);

EFF_C1CLSStripdet5= std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSStripdet5"," ",129,-0.5,128.5)) ;
EFF_C1CLSStripdet5->SetDirectory(0);
EFF_C1CLSPaddet5=  std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSPaddet5"," ",129,-0.5,128.5)) ;
EFF_C1CLSPaddet5->SetDirectory(0);

EFF_C1CLSStripdet6= std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSStripdet6"," ",129,-0.5,128.5)) ;
EFF_C1CLSStripdet6->SetDirectory(0);
EFF_C1CLSPaddet6=  std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSPaddet6"," ",129,-0.5,128.5)) ;
EFF_C1CLSPaddet6->SetDirectory(0);

EFF_C1CLSStripdet7= std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSStripdet7"," ",129,-0.5,128.5)) ;
EFF_C1CLSStripdet7->SetDirectory(0);
EFF_C1CLSPaddet7=  std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSPaddet7"," ",129,-0.5,128.5)) ;
EFF_C1CLSPaddet7->SetDirectory(0);

EFF_C1CLSStripdet8= std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSStripdet8"," ",129,-0.5,128.5)) ;
EFF_C1CLSStripdet8->SetDirectory(0);
EFF_C1CLSPaddet8=  std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSPaddet8"," ",129,-0.5,128.5)) ;
EFF_C1CLSPaddet8->SetDirectory(0);

EFF_C1CLSStripdet9= std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSStripdet9"," ",129,-0.5,128.5)) ;
EFF_C1CLSStripdet9->SetDirectory(0);
EFF_C1CLSPaddet9=  std::auto_ptr<TH1F>(new TH1F("EFF_C1CLSPaddet9"," ",129,-0.5,128.5)) ;
EFF_C1CLSPaddet9->SetDirectory(0);

RelEffi_ForTracking=  std::auto_ptr<TH1F>(new TH1F("RelEffi_ForTracking"," Plane Contribution to the Efficiency Tracking",11,-0.5,10.5)) ;
RelEffi_ForTracking->SetDirectory(0);


PSimHitYD1= std::auto_ptr<TH1F>(new TH1F("PSimHitXD1","PSimHitYD1",300,-150.,150.));
PSimHitYD1->SetDirectory(0);
PSimHitXD1= std::auto_ptr<TH1F>(new TH1F("PSimHitXD1","PSimHitXD1",300,-150.,150.));
PSimHitXD1->SetDirectory(0);
 
 Plside0PSimHitXD0= std::auto_ptr<TH1F>(new TH1F(" Plside0PSimHitXD0"," Plside0PSimHitXD0",300,-150.,150.));
 Plside0PSimHitXD0->SetDirectory(0);
 Plside1PSimHitXD0= std::auto_ptr<TH1F>(new TH1F(" Plside1PSimHitXD0"," Plside1PSimHitXD0",300,-150.,150.));
 Plside1PSimHitXD0->SetDirectory(0);
 Plside0PSimHitYD0= std::auto_ptr<TH1F>(new TH1F(" Plside0PSimHitYD0"," Plsid00PSimHitYD0",300,-150.,150.));
 Plside0PSimHitYD0->SetDirectory(0);
 Plside1PSimHitYD0= std::auto_ptr<TH1F>(new TH1F(" Plside1PSimHitYD0"," Plside1PSimHitYD0",300,-150.,150.));
 Plside1PSimHitYD0->SetDirectory(0);





  DetIdvsXShOnlyXfit= std::auto_ptr<TProfile>(new TProfile("DetIdvsXShOnlyXfit","Det Id vs X Shift (OnlyXfit)",42,-1.5,40.5,""));
  DetIdvsXShOnlyXfit->SetDirectory(0);

  DetIdvsYShOnlyYfit= std::auto_ptr<TProfile>(new TProfile("DetIdvsYShOnlyYfit","Det Id vs Y Shift (OnlyYfit)",42,-1.5,40.5,""));
  DetIdvsYShOnlyYfit->SetDirectory(0);

ProbYHisto= std::auto_ptr<TH1F>(new TH1F("ProbYHisto","ProbYHisto",200,0.,1.));
ProbYHisto->SetDirectory(0);
ProbXHisto= std::auto_ptr<TH1F>(new TH1F("ProbXHisto","ProbXHisto",200,0.,1.));
ProbXHisto->SetDirectory(0);
EventCandidate= std::auto_ptr<TH1F>(new TH1F("EventCandidatesl","EventCandidatesl",3,-0.5,2.5));
ProbXHistoAl= std::auto_ptr<TH1F>(new TH1F("ProbXHistoAl","ProbXHistoAl",200,0.,1.));
ProbXHistoAl->SetDirectory(0);
ProbYHistoAl= std::auto_ptr<TH1F>(new TH1F("ProbYHistoAl","ProbYHistoAl",200,0.,1.));
ProbYHistoAl->SetDirectory(0);

AXHisto= std::auto_ptr<TH1F>(new TH1F("AXHisto","AXHisto",100,-0.1,0.1));
AXHisto->SetDirectory(0);
AYHisto= std::auto_ptr<TH1F>(new TH1F("AYHisto","AYHisto",100,-0.1,0.1));
AYHisto->SetDirectory(0);

AXError= std::auto_ptr<TH1F>(new TH1F("AXError","AXError",100,-0.05,0.05));
AXError->SetDirectory(0);
AYError= std::auto_ptr<TH1F>(new TH1F("AYError","AYError",100,-0.05,0.05));
AYError->SetDirectory(0);
BXError= std::auto_ptr<TH1F>(new TH1F("BXError","BXError",100,-100.,100.));
BXError->SetDirectory(0);
BYError= std::auto_ptr<TH1F>(new TH1F("BYError","BYError",100,-1000.,1000.));
BYError->SetDirectory(0);
 HitYErr= std::auto_ptr<TH1F>(new TH1F("HitYErr","HitYErr",100,-1,1));
HitYErr->SetDirectory(0);
HitXErr= std::auto_ptr<TH1F>(new TH1F("HitXErr","HitXErr",100,-1,1));
HitXErr->SetDirectory(0);

HitX= std::auto_ptr<TH1F>(new TH1F("HitX","HitX",300,-150,150));
HitX->SetDirectory(0);
HitY= std::auto_ptr<TH1F>(new TH1F("HitY","HitY",300,-150,150));
HitY->SetDirectory(0);

TrkChi2X= std::auto_ptr<TH1F>(new TH1F("TrkChi2X","TrkChi2X",800,0.,800.));
TrkChi2X->SetDirectory(0);
TrkChi2Y= std::auto_ptr<TH1F>(new TH1F("TrkChi2Y","TrkChi2Y",800,0.,800.));
TrkChi2Y->SetDirectory(0);


HitMatrixR= std::auto_ptr<TH1F>(new TH1F("HitMatrixR","Class 1 Hit Matrix R",240,30.,150.)); 
HitMatrixR->SetDirectory(0);
HitMatrixR->SetXTitle("R (mm)");

HitMatrixPhi= std::auto_ptr<TH1F>(new TH1F("HitMatrixPhi","Class 1 Hit Phi",360,0.,360.)); 
HitMatrixPhi->SetDirectory(0);
 HitMatrixPhi->SetXTitle("Deg");

 TrackPhi0R0= std::auto_ptr<TH2F>(new TH2F("TrackPhi0R0","Track #phi-R Entry Point utilized in efficiency",360,0.,360.,120,30.,150.)); 
 TrackPhi0R0->SetDirectory(0);
 //TrackPhi0R0->SetOption("lego");
 char sZname2[1024];
 char sZnamehist[1024];

  sprintf(sZnamehist,"Hit Noise level vs Det Id (DR= %f D#phi= %f)", Effmaxdrhit,Effmaxdphihit);
  HITNoiseLevelvsDet=  std::auto_ptr<TH1F>(new TH1F("HITNoiseLevelvsDet",sZnamehist,12,-0.5,11.5)); 
  HITNoiseLevelvsDet->SetDirectory(0);
  sprintf(sZnamehist,"Pad-Cluster Noise level vs Det Id (D#phi= %f)", Effmaxdphihit);
  PADNoiseLevelvsDet=  std::auto_ptr<TH1F>(new TH1F("PADNoiseLevelvsDet",sZnamehist,12,-0.5,11.5));  
  PADNoiseLevelvsDet->SetDirectory(0);
  sprintf(sZnamehist,"Strip-Cluster Noise level vs Det Id (DR= %f)", Effmaxdrhit);
  STRIPNoiseLevelvsDet=  std::auto_ptr<TH1F> (new TH1F("STRIPNoiseLevelvsDet",sZnamehist,12,-0.5,11.5)); 
  STRIPNoiseLevelvsDet->SetDirectory(0);

 
  Average_Cl1HitNoiseNumb_vsDet =  std::auto_ptr<TProfile>(new TProfile("Mean Number of Class1 Noise per event vs Det when Tracking Plane","Average_Cl1HitNoiseNumb_vsDet",12,-0.5,11.5)); 
  Average_Cl1HitNoiseNumb_vsDet->SetDirectory(0);

//std::auto_ptr<TH2F>  RPhiEvtHit[10];
//muEtaEnergy= std::auto_ptr<TH2F>(new TH2F("muEtaEnergy","#eta vs Energy muons",30,4.5,8.5,40,1.,400.));
//muEtaEnergy->SetDirectory(0);



  VfatsCorruptionFrequency2D=std::auto_ptr<TH2D>(new TH2D("VfatsCorruptionFrequency2D","VfatsCorruptionFrequency2D",40,-0.5,39.5,17,-0.5,16.5));
  VfatsCorruptionFrequency2D->SetDirectory(0);

  T2vfatinfoErrorMapNOTDEAD=std::auto_ptr<TH2D>(new TH2D("T2vfatinfoErrorMapNOTDEAD","Map of corrupted frame, only alive vfat",40,-0.5,39.5,17,-0.5,16.5));
  T2vfatinfoErrorMapNOTDEAD->GetXaxis()->SetTitle("Plane") ;
  T2vfatinfoErrorMapNOTDEAD->GetYaxis()->SetTitle("VFAT position") ;
  T2vfatinfoErrorMapNOTDEAD->SetDrawOption("text");
  T2vfatinfoErrorMapNOTDEAD->SetDirectory(0);

  T2vfatinfoErrorMap_DEADWithCorrectFrame=std::auto_ptr<TH2D>(new TH2D("T2vfatinfoErrorMap_DEADWithCorrectFrame","Map of dead vfat with correct frames",40,-0.5,39.5,17,-0.5,16.5));
  T2vfatinfoErrorMap_DEADWithCorrectFrame->GetXaxis()->SetTitle("Plane") ;
  T2vfatinfoErrorMap_DEADWithCorrectFrame->GetYaxis()->SetTitle("VFAT position") ;

  T2vfatinfoErrorMap=std::auto_ptr<TH2D>(new TH2D("T2vfatinfoErrorMap","Map of missing vfat frames",40,-0.5,39.5,17,-0.5,16.5));
  T2vfatinfoErrorMap->GetXaxis()->SetTitle("Plane") ;
  T2vfatinfoErrorMap->GetYaxis()->SetTitle("VFAT position") ;
  T2vfatinfoErrorMap->SetDrawOption("text");
  T2vfatinfoErrorMap->SetDirectory(0);

  T2vfatinfoErrorMapGeoQ1=std::auto_ptr<TH2D>(new TH2D("T2vfatinfoErrorMapGeoQ1","Geometrical view of corruption",13,-0.5,12.5,2,-0.5,1.5));
  T2vfatinfoErrorMapGeoQ1->GetXaxis()->SetTitle("Azimuthal Position");
  T2vfatinfoErrorMapGeoQ1->GetYaxis()->SetTitle("Radial Position");
  

  T2vfatinfoErrorMapFootPr=std::auto_ptr<TH2D>(new TH2D("T2vfatinfoErrorMapFootPr","T2vfatinfoErrorMapFootPr",40,-0.5,39.5,17,-0.5,16.5));
  T2vfatinfoErrorMapCRC=std::auto_ptr<TH2D>(new TH2D("T2vfatinfoErrorMapCRC","T2vfatinfoErrorMapCRC",40,-0.5,39.5,17,-0.5,16.5));
  T2vfatinfoErrorNotInData=std::auto_ptr<TH2D>(new TH2D("T2vfatinfoErrorNotInData","T2vfatinfoErrorNotInData",40,-0.5,39.5,17,-0.5,16.5));
  T2vfatinfoErrorMapNotIncluded=std::auto_ptr<TH2D>(new TH2D("T2vfatinfoErrorMapNotIncluded","T2vfatinfoErrorMapNotIncluded",40,-0.5,39.5,17,-0.5,16.5));
  T2vfatinfoErrorErrorMapping=std::auto_ptr<TH2D>(new TH2D("T2vfatinfoErrorErrorMapping","T2vfatinfoErrorErrorMapping",40,-0.5,39.5,17,-0.5,16.5));




 for(unsigned int m=0;m<4; m++)
   {

     sprintf(sZname2, "EFFpadVsIdNormalized_sector%d", m); 
     sprintf(sZnamehist, "Pad Efficiency vs Detector-Id Normalized R-Sector-%d", m); 
     EFFpadVsIdNormalized_sectorX[m]=std::auto_ptr<TProfile>(new TProfile(sZname2,sZnamehist,12,-1.5,10.));
     EFFpadVsIdNormalized_sectorX[m]->SetDirectory(0);

     sprintf(sZname2, "EFFstripVsIdNormalized_sector%d", m); 
     sprintf(sZnamehist, "Strip Efficiency vs Detector-Id Normalized R-Sector-%d", m); 
     EFFstripVsIdNormalized_sectorX[m]=std::auto_ptr<TProfile>(new TProfile(sZname2,sZnamehist,12,-1.5,10.));
     EFFstripVsIdNormalized_sectorX[m]->SetDirectory(0);
     
   }





 PlaneNumPadClu11_vs_16=std::auto_ptr<TH2D>(new TH2D("PlaneNumPadClu11_vs_16","PlaneNumPadClu11_vs_16",200,-0.5,199.5,200,-0.5,199.5));
 PlaneNumPadClu11_vs_16->GetXaxis()->SetTitle("H1 plane 1") ;
 PlaneNumPadClu11_vs_16->GetYaxis()-> SetTitle("H1 plane 6") ;

  for(unsigned int m=0;m<40; m++)
   {


     sprintf(sZname2, "HGeometry_StripEfficiency_Num %d", m); 
     sprintf(sZnamehist, "HGeometry_StripEfficiency_Num %d", m);
     
     HGeometry_StripEfficiency_Num[m]=std::auto_ptr<TH2D>(new TH2D(sZname2,sZnamehist,numRsectEffi,-0.5,numRsectEffi-0.5,numPhisectEffi,-0.5,numPhisectEffi-0.5));

     sprintf(sZname2, "HGeometry_StripEfficiency_Den %d", m); 
     sprintf(sZnamehist, "HGeometry_StripEfficiency_Den %d", m);
     HGeometry_StripEfficiency_Den[m]=std::auto_ptr<TH2D>(new TH2D(sZname2,sZnamehist,numRsectEffi,-0.5,numRsectEffi-0.5,numPhisectEffi,-0.5,numPhisectEffi-0.5)); 

     sprintf(sZname2, "HGeometry_PadEfficiency_Num %d", m); 
     sprintf(sZnamehist, "HGeometry_PadEfficiency_Num %d", m);
     HGeometry_PadEfficiency_Num[m]=std::auto_ptr<TH2D>(new TH2D(sZname2,sZnamehist,numRsectEffi,-0.5,numRsectEffi-0.5,numPhisectEffi,-0.5,numPhisectEffi-0.5));

     sprintf(sZname2, "HGeometry_PadEfficiency_Den %d", m); 
     sprintf(sZnamehist, "HGeometry_PadEfficiency_Den %d", m);
     HGeometry_PadEfficiency_Den[m]=std::auto_ptr<TH2D>(new TH2D(sZname2,sZnamehist,numRsectEffi,-0.5,numRsectEffi-0.5,numPhisectEffi,-0.5,numPhisectEffi-0.5)); 



     sprintf(sZname2, "PlaneHitRadiographyXY %d", m); 
     sprintf(sZnamehist, "Cluster Pad  X-Y in plane %d", m); 
     PlaneHitRadiographyXY[m]=std::auto_ptr<TH2D>(new TH2D(sZname2,sZnamehist,70,-144.5,144.5,70,-144.5,144.5));
     PlaneHitRadiographyXY[m]->GetXaxis()->SetTitle("X") ;
     PlaneHitRadiographyXY[m]->GetYaxis()-> SetTitle("Y") ;
     PlaneHitRadiographyXY[m]->SetDirectory(0);
     
     sprintf(sZname2, "PadClusterR %d", m); 
     sprintf(sZnamehist, "Cluster Pad  R in plane %d", m); 
     PadClusterR[m]=std::auto_ptr<TH1D>(new TH1D(sZname2,sZnamehist,110,40.,150.)); 
  
     sprintf(sZname2, "PadClusterXYRadiographyYMinus %d", m); 
     sprintf(sZnamehist, "Cluster Pad  XY Y<0 in plane %d", m); 
     PadClusterXYRadiographyYMinus[m]=std::auto_ptr<TH2D>(new TH2D(sZname2,sZnamehist,70,-144.5,144.5,70,-144.5,144.5));

     sprintf(sZname2, "PadClusterXYRadiographyYPlus %d", m); 
     sprintf(sZnamehist, "Cluster Pad  XY Y>0 in plane %d", m); 
     PadClusterXYRadiographyYPlus[m]=std::auto_ptr<TH2D>(new TH2D(sZname2,sZnamehist,70,-144.5,144.5,70,-144.5,144.5));



     sprintf(sZname2, "StripCluRadiographyXY %d", m); 
     sprintf(sZnamehist, "Cluster Strip  X-Y in plane %d", m); 
     StripCluRadiographyXY[m]=std::auto_ptr<TH2D>(new TH2D(sZname2,sZnamehist,70,-144.5,144.5,70,-144.5,144.5));
     StripCluRadiographyXY[m]->GetXaxis()->SetTitle("X") ;
     StripCluRadiographyXY[m]->GetYaxis()-> SetTitle("Y") ;
     StripCluRadiographyXY[m]->SetDirectory(0);


     sprintf(sZname2, "VfatsCorruptionFrequency %d", m); 
     sprintf(sZnamehist, "Vfats Corruption Frequency in Plane %d", m);
     VfatsCorruptionFrequency[m]=std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,18,-0.5,17.5));
     VfatsCorruptionFrequency[m]->SetDirectory(0);

     sprintf(sZname2, "PadClusterR_AllvsPlane %d", m); 
     sprintf(sZnamehist, "PadClusterR_AllvsPlane %d", m);
     PadClusterR_AllvsPlane[m]=std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,120,30.,150.)); 

     sprintf(sZname2, "PadClusterSize_AllvsPlane %d", m); 
     sprintf(sZnamehist, "PadClusterSize_AllvsPlane %d", m);
     PadClusterSize_AllvsPlane[m]=std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,28,-0.5,27.5));
     

    
     if(m<10)
       {
	 sprintf(sZname2, "PadClusterXYForEffiEvents %d", m); 
	 sprintf(sZnamehist, "Cluster Pad  X-Y in plane %d utilized for EFFI", m); 
	 PadClusterXYForEffiEvents[m]=std::auto_ptr<TH2D>(new TH2D(sZname2,sZnamehist,70,-144.5,144.5,70,-144.5,144.5));
	 PadClusterXYForEffiEvents[m]->GetXaxis()->SetTitle("X") ;
	 PadClusterXYForEffiEvents[m]->GetYaxis()-> SetTitle("Y") ;
	 PadClusterXYForEffiEvents[m]->SetDirectory(0);
       }
     
     
     if((m%10)==0)
       {
	 unsigned int hindex=m/10;
	 sprintf(sZname2, "HalfTeleTrkRadiographyXY %d", hindex); 
	 sprintf(sZnamehist, "Track X-Y entry point in quarter H %d", hindex); 
	 HalfTeleTrkRadiographyXY[hindex]=std::auto_ptr<TH2D>(new TH2D(sZname2,sZnamehist,70,-144.5,144.5,70,-144.5,144.5));
	 HalfTeleTrkRadiographyXY[hindex]->GetXaxis()->SetTitle("X") ;
	 HalfTeleTrkRadiographyXY[hindex]->GetYaxis()-> SetTitle("Y") ;
	 HalfTeleTrkRadiographyXY[hindex]->SetDirectory(0);
       }
   }


  


  for(unsigned int m=0;m<40; m++)
    for(unsigned int k=0;k<17; k++)
      {

	sprintf(sZname2, "vFatCumulative %d-%d", m,k); 
	sprintf(sZnamehist, "vFatCumulative %d-%d", m,k);	
	vFatCumulative[m][k]=std::auto_ptr<TH1D>(new TH1D(sZname2,sZnamehist,131,-0.5,130.5));

	sprintf(sZname2, "vFatMultiplicty %d-%d", m,k); 
	sprintf(sZnamehist, "vFatMultiplicty %d-%d", m,k);
	vFatMultiplicty[m][k]=std::auto_ptr<TH1D>(new TH1D(sZname2,sZnamehist,131,-0.5,130.5));
	
	sprintf(sZname2, "vFatMultiplictyVsTrk %d-%d", m,k); 
	sprintf(sZnamehist, "vFatMultiplictyVsTrk %d-%d", m,k);
	vFatMultiplictyVsTrk[m][k]=std::auto_ptr<TProfile>(new TProfile(sZname2,sZnamehist,151,-0.5,150.5));

	sprintf(sZname2, "TrkMultiplictyWhenEveryThingIsOn %d-%d", m,k); 
	sprintf(sZnamehist, "TrkMultiplictyWhenEveryThingIsOn %d-%d", m,k);
	TrkMultiplictyWhenEveryThingIsOn[m][k]=std::auto_ptr<TH1D>(new TH1D(sZname2,sZnamehist,131,-0.5,130.5));	
      }

  //Strip only, max is 4*10=40
  NumH1VfatCompletelyOnPerEvt=std::auto_ptr<TH1D> (new TH1D("NumH1VfatCompletelyOnPerEvt","NumH1VfatCompletelyOnPerEvt",41,-0.5,40.5));

  for(unsigned int k=0;k<5; k++)
    {
      sprintf(sZname2, "Vfat15_Pl10andH1Vfat15_PlDispariMultCorrel %d", k); 
      sprintf(sZnamehist, "Vfat15_Pl10andH1Vfat15_PlDispariMultCorrel %d",k);      
      Vfat15_Pl10andH1Vfat15_PlDispariMultCorrel[k]= std::auto_ptr<TH2D>(new TH2D(sZname2,sZnamehist,130,-0.5,129.5,130,-0.5,129.5));   


      sprintf(sZname2, "Vfat16_Pl10andH1Vfat16_PlDispariMultCorrel %d", k); 
      sprintf(sZnamehist, "Vfat16_Pl10andH1Vfat16_PlDispariMultCorrel %d",k);      
      Vfat16_Pl10andH1Vfat16_PlDispariMultCorrel[k]= std::auto_ptr<TH2D>(new TH2D(sZname2,sZnamehist,130,-0.5,129.5,130,-0.5,129.5));   
  	
    } 
 
  Vfat15_Pl19and11_Pl19MultCorrel= std::auto_ptr<TH2D>(new TH2D("Vfat15_Pl19and11_Pl19MultCorrel","Vfat15_Pl19and11_Pl19MultCorrel",130,-0.5,129.5,130,-0.5,129.5));   
  
  Vfat15_Pl19and17MultCorrel= std::auto_ptr<TH2D>(new TH2D("Vfat15_Pl19and17MultCorrel","Vfat15_Pl19and17MultCorrel",130,-0.5,129.5,130,-0.5,129.5));   

  Vfat15_Pl19and17MultCorrel_LowMulti= std::auto_ptr<TH2D>(new TH2D("Vfat15_Pl19and17MultCorrel_LowMulti","Vfat15_Pl19and17MultCorrel_LowMulti",130,-0.5,129.5,130,-0.5,129.5));   

  //-----------------------------------------------------------------------------------------------------
  // Noise correlation Histograms Begin
  
  string EvtVfat_Strip_WhenCompletelyOn_Title="HistoCorrelationNoise";
  


  NoisyStripVfatH1_PlIId.clear();
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(11,0));
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(11,1));  
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(11,15));
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(11,16));
  
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(13,0));
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(13,1));  
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(13,15));
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(13,16));
  
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(15,0));
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(15,1));  
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(15,15));
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(15,16));
  
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(17,0));
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(17,1));  
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(17,15));
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(17,16));
  
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(19,0));
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(19,1));  
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(19,15));
  NoisyStripVfatH1_PlIId.push_back(std::pair<int,int>(19,16));

  EvtVfat_Strip_IdWhenCompletelyOn_ID=std::auto_ptr<TH1D>(new TH1D("EvtVfat_Strip_IdWhenCompletelyOn_ID","EvtVfat_Strip_IdWhenCompletelyOn_ID",5000,-0.5,4999.5));//This is the (100*plane + iid) VFAT id for the noisy vfat considered in the noise-correlation
  unsigned int uniqueid=0;
  for(unsigned int g=0;g<NoisyStripVfatH1_PlIId.size();g++)
    {
      uniqueid=NoisyStripVfatH1_PlIId.at(g).first*100+NoisyStripVfatH1_PlIId.at(g).second;
      string axlabel="";
      std::stringstream out;
      out << uniqueid;
      axlabel= out.str();
      axlabel=axlabel+";";
      EvtVfat_Strip_WhenCompletelyOn_Title=EvtVfat_Strip_WhenCompletelyOn_Title+axlabel;
      EvtVfat_Strip_IdWhenCompletelyOn_ID->Fill(uniqueid);
   }

  std::cout<<"X axis label:"<<EvtVfat_Strip_WhenCompletelyOn_Title.c_str()<<std::endl;

  Int_t *bins = new Int_t[NoisyStripVfatH1_PlIId.size()];
  Double_t *xmin = new Double_t[NoisyStripVfatH1_PlIId.size()];
  Double_t *xmax = new Double_t[NoisyStripVfatH1_PlIId.size()];

  for(unsigned int kk=0; kk<NoisyStripVfatH1_PlIId.size();kk++){
    bins[kk]=2;
    xmin[kk]=-0.5;
    xmax[kk]=1.5;
  }
 
  
  Int_t dimension=NoisyStripVfatH1_PlIId.size();
  
  EvtVfat_Strip_WhenCompletelyOn = std::auto_ptr<THnSparseD>(new THnSparseD("EvtVfat_Strip_WhenCompletelyOn",EvtVfat_Strip_WhenCompletelyOn_Title.c_str(),dimension,bins,xmin,xmax));

  EvtVfat_Strip_WhenCompletelyOnS = std::auto_ptr<THnSparseS>(new THnSparseS("EvtVfat_Strip_WhenCompletelyOnS",EvtVfat_Strip_WhenCompletelyOn_Title.c_str(),dimension,bins,xmin,xmax));


for(unsigned int g=0;g<NoisyStripVfatH1_PlIId.size();g++)
    {
      uniqueid=NoisyStripVfatH1_PlIId.at(g).first*100+NoisyStripVfatH1_PlIId.at(g).second;
      TString axlabel;//="";
      std::stringstream out;
      out << uniqueid;
      axlabel= (TString) out.str();
      std::cout<<"Imosing axis name "<<g<<" "<<axlabel<<std::endl;
      EvtVfat_Strip_WhenCompletelyOn->GetAxis(g)->SetTitle(axlabel);//Name
      TAxis* currentAx=  EvtVfat_Strip_WhenCompletelyOn->GetAxis(g);
      const char* thetile = currentAx->GetTitle();
      std::cout<<"Axes Title: "<<g<<"-"<<thetile<<std::endl;
    }


  // Noise correlation Histograms End
  //-----------------------------------------------------------------------------------------------------
 

  /*
  They take name and title, the number of dimensions, and for each dimension
  the number of bins, the minimal, and the maximal value on the dimension's
  axis. A TH2 h("h","h",10, 0., 10., 20, -5., 5.) would correspond to
  Int_t bins[2] = {10, 20};
  Double_t xmin[2] = {0., -5.};
  Double_t xmax[2] = {10., 5.};
  THnSparse hs("hs", "hs", 2, bins, min, max);

  * Filling
  A THnSparse is filled just like a regular histogram, using
  THnSparse::Fill(x, weight), where x is a n-dimensional Double_t value.
  */

  
  for(unsigned int m=0;m<10; m++)
   {

     sprintf(sZname2, "TrkHitR_vsplane %d", m); 
     sprintf(sZnamehist, "TrkHitR_vsplane %d", m);
     TrkHitR_vsplane[m]=std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,120,30.,150.));

     sprintf(sZname2, "VfatEfficiency_Plane %d", m); 
     sprintf(sZnamehist, "Vfat Efficiency in Plane %d", m);
     VFATEFF[m]= std::auto_ptr<TProfile>(new TProfile(sZname2,sZnamehist,18,-0.5,17.5));
     VFATEFF[m]->SetDirectory(0); 


     sprintf(sZname2, "VfatEfficiency_Normalized_Plane %d", m); 
     sprintf(sZnamehist, "Vfat Efficiency in Plane %d (Normalized)", m);
     VFATEFFNormalized[m]= std::auto_ptr<TProfile>(new TProfile(sZname2,sZnamehist,18,-0.5,17.5));
     VFATEFFNormalized[m]->SetDirectory(0); 

     sprintf(sZname2, "VfatEffSTAT_Plane %d", m); 
     sprintf(sZnamehist, "Vfat statistics for Eff. calc, plane %d", m);
     VfatStatistics[m]= std::auto_ptr<TH1F> (new TH1F(sZname2,sZnamehist,18,-0.5,17.5)); 
     VfatStatistics[m]->SetDirectory(0); 

     sprintf(sZname2, "Class1NoiseDRDPHISep %d", m); 
     sprintf(sZnamehist, "Class1Noise DR-D#PHI Separation from tracking Hit, plane %d", m);
     Class1NoiseDRDPHISep[m] = std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,450,30,150,360,0,360));   

     sprintf(sZname2, "TrackingHitVsClass1NoiseDRSep %d", m); 
     sprintf(sZnamehist, "Tracking Hit R Vs Distance from Class1-Hit noise, plane %d", m);
     TrackingHitVsClass1NoiseDRSep[m]= std::auto_ptr<TProfile>(new TProfile(sZname2,sZnamehist,81,-140.5,140.5));

     sprintf(sZname2, "TrackingHitRVsStripCluNoiseDRSep %d", m); 
     sprintf(sZnamehist, "Tracking Hit R Vs Distance from Class 0 Strip Cluster noise, plane  %d", m);
     TrackingHitRVsStripCluNoiseDRSep[m]= std::auto_ptr<TProfile>(new TProfile(sZname2,sZnamehist,81,-140.5,140.5));

     sprintf(sZname2, "TrackingHitRVsPadCluNoiseDRSep %d", m); 
     sprintf(sZnamehist, "Tracking Hit R Vs Distance from Class 0 Pad Cluster noise, plane %d", m);
     TrackingHitRVsPadCluNoiseDRSep[m]= std::auto_ptr<TProfile>(new TProfile(sZname2,sZnamehist,81,-140.5,140.5));

 


       sprintf(sZname2, "RPhiHit Plane %d", m); 
       sprintf(sZnamehist, "RPhiHit Plane %d", m);
       RPhiEvtHit[m] = std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,360,0,360,480,30,150));
       RPhiEvtHit[m]->SetDirectory(0);

       sprintf(sZname2, "NoiseHitDet %d", m); 
       sprintf(sZnamehist, "Number of Class 1 Hit in det %d when the particle is detected", m);
       NoiseHitDet[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,200,-0.5,199.5));
       NoiseHitDet[m]->SetDirectory(0);
       sprintf(sZname2, "Noise01HitDet %d", m); 
       sprintf(sZnamehist, "Presence of Noisy Class 1 Hit in det %d when the particle is detected", m);
       Noise01HitDet[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,200,-0.5,199.5));
       Noise01HitDet[m]->SetDirectory(0);

 
       sprintf(sZname2, "NoiseStripDet %d", m); 
       sprintf(sZnamehist, "Number of Class 1 Strip in det %d when the particle is detected", m);
       NoiseStripDet[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,100,-0.5,99.5));
       NoiseStripDet[m]->SetDirectory(0);
       sprintf(sZname2, "Noise01StripDet %d", m); 
       sprintf(sZnamehist, "Presence of Noisy Class 1 Strip in det %d when the particle is detected", m);
       Noise01StripDet[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,100,-0.5,99.5));
       Noise01StripDet[m]->SetDirectory(0);

       
       sprintf(sZname2, "NoisePadDet %d", m); 
       sprintf(sZnamehist, "Number of Class 1 Pad in det %d when the particle is detected", m);
       NoisePadDet[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,500,-0.5,499.5));
       NoisePadDet[m]->SetDirectory(0);
       sprintf(sZname2, "Noise01PadDet %d", m); 
       sprintf(sZnamehist, "Presence of Noisy Class 1 Pad in det %d when the particle is detected", m);
       Noise01PadDet[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,500,-0.5,499.5));
       Noise01PadDet[m]->SetDirectory(0);


       sprintf(sZname2, "DXAlignDet %d", m); 
       sprintf(sZnamehist, "DX shifts det %d", m);
       DXAlignDet[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,20,-3,3));
       DXAlignDet[m]->SetDirectory(0);
       
       sprintf(sZname2, "DYAlignDet %d", m); 
       sprintf(sZnamehist, "DY shifts det %d", m);
       DYAlignDet[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,20,-3,3));
       DYAlignDet[m]->SetDirectory(0);

       
       sprintf(sZname2, "DRResol %d", m); 
       sprintf(sZnamehist, "R resolution det %d", m);
       DRResol[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,60,-2,2));
       DRResol[m]->SetDirectory(0);
       DRResol[m]->SetXTitle("#Delta R (mm)");

       sprintf(sZname2, "DPhiResol %d", m); 
       sprintf(sZnamehist, "Phi resolution det %d", m);
       // DPhiResol[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,30,-5,5));
       DPhiResol[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,60,-5,5));
       DPhiResol[m]->SetDirectory(0);
       DPhiResol[m]->SetXTitle("#Delta #Phi (deg)");

       sprintf(sZname2, "drHit12 %d", m); 
       sprintf(sZnamehist, "Particle Hit R - Noise Hit R (first two CL-1 hits), Det %d", m);
       drHit12[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,150,-150,150));
       drHit12[m]->SetDirectory(0);
       drHit12[m]->SetXTitle("DR (mm)");


       sprintf(sZname2, "DXResp0 %d", m); 
       sprintf(sZnamehist,"X schift det %d respect det 0", m);
       DXResp0[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,200,-10,10));
       DXResp0[m]->SetDirectory(0);
       sprintf(sZname2, "DYResp0 %d", m); 
       sprintf(sZnamehist,"Y schift det %d respect det 0", m);
       DYResp0[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,200,-10,10));
       DYResp0[m]->SetDirectory(0);

       sprintf(sZname2, "DXResp9 %d", m); 
       sprintf(sZnamehist,"X schift det %d respect det 9", m);
       DXResp9[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,200,-10,10));
       DXResp9[m]->SetDirectory(0);
       sprintf(sZname2, "DYResp9 %d", m); 
       sprintf(sZnamehist,"Y schift det %d respect det 9", m);
       DYResp9[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,200,-10,10));
       DYResp9[m]->SetDirectory(0);

       sprintf(sZname2, "PAdNr %d", m); 
       sprintf(sZnamehist,"PAdNr %d ", m);
       PAdNr[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,1600,-0.5,1599.5)) ;
       PAdNr[m]->SetDirectory(0);

       sprintf(sZname2, "StripNr %d", m); 
       sprintf(sZnamehist,"StripNr %d ", m);
       StripNr[m]= std::auto_ptr<TH1F>(new TH1F(sZname2,sZnamehist,560,-0.5,559.5)) ;
       StripNr[m]->SetDirectory(0);
   } 


   CumulativePhiResol= std::auto_ptr<TH1F>(new TH1F("CumulativePhiResol","Azimuthal Hit Resolution",20,-5,5));
   CumulativePhiResol->SetDirectory(0);
   CumulativeRResol= std::auto_ptr<TH1F>(new TH1F("CumulativeRResol","Radial Hit Resolution",40,-2,2));
   CumulativeRResol->SetDirectory(0);
   CumulativeExpXUncert= std::auto_ptr<TH1F>(new TH1F("CumulativeExpXUncert","CumulativeExpXUncert",100,-10,10));
   CumulativeExpYUncert= std::auto_ptr<TH1F>(new TH1F("CumulativeExpYUncert","CumulativeExpYUncert",100,-10,10));
   
   CumulativeRResol->SetXTitle("#Delta R (mm)");
   CumulativePhiResol->SetXTitle("#Delta #Phi (deg)");

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

  EFF_DPhiHitNoCutdet2->SetXTitle("D#phi (deg)");
  EFF_DRHitNoCutdet2->SetXTitle("Dr (mm)"); 
  PolarAngles->SetXTitle("(rad)"); 




  //std::cout<<"End Beg Job"<<std::endl;

}






































// ------------ method called once each job just after ending the event loop  ------------

void T2AnalyzerRaw::endJob()
{
 
   std::cout<<" End Job. Number of non-corrupted event: "<<countegood<<std::endl;
   
  std::string namefileEffSec=outputFileName;
  std::string sought1 = "root";
  int pos1=namefileEffSec.find(sought1); 
  namefileEffSec.replace(pos1, sought1.size(), "dat"); 
  std::cout<<" Opening file for store Sector Efficiency Plot "<<namefileEffSec.c_str()<<std::endl;

  //  fstream file_opSec(namefileEffSec.c_str(), ios::out); 
  //Warning!! Found very strange problem opening file with "in" flag,file is written only if it already exists in the path.
  //fstream file_op("test.dat",ios::in | ios::out);
  //  if(((file_opSec).is_open())==false)
  //  {
  //    std::cout<<" ERROR: Problem in opening file file_opSec"<<std::endl;
  //    (file_opSec).clear();
  //  }

  for(unsigned int j=0;j<40;j++)
    for( int r=0;r<numRsectEffi;r++)
      for( int p=0;p<numPhisectEffi;p++){

	if(Geometry_PadEfficiency_Den.at(j).at(r).at(p)>0)
	  Geometry_PadEfficiency.at(j).at(r).at(p)=Geometry_PadEfficiency_Num.at(j).at(r).at(p)/Geometry_PadEfficiency_Den.at(j).at(r).at(p);
	if(Geometry_StripEfficiency_Den.at(j).at(r).at(p)>0)
	  Geometry_StripEfficiency.at(j).at(r).at(p)=Geometry_StripEfficiency_Num.at(j).at(r).at(p)/Geometry_StripEfficiency_Den.at(j).at(r).at(p);
	

	//file_opSec<<"DetId: "<<j<<"  CluType: Pad   RSec: "<<r<<"   PhiSec: "<<p<<" Effi: "<< Geometry_PadEfficiency.at(j).at(r).at(p) <<" ) \n";
	//file_opSec<<"DetId: "<<j<<"  CluType: Strip   RSec: "<<r<<"   PhiSec: "<<p<<" Effi: "<< Geometry_StripEfficiency.at(j).at(r).at(p) <<" ) \n";
      }

  //  file_opSec.close(); 


  //std::vector<std::vector<std::vector<double> > > Geometry_PadEfficiency_Num;
  //std::vector<std::vector<std::vector<double> > > Geometry_PadEfficiency_Den;
  //std::vector<std::vector<std::vector<double> > > Geometry_PadEfficiency;
  
  //std::vector<std::vector<std::vector<double> > > Geometry_StripEfficiency_Num;
  //std::vector<std::vector<std::vector<double> > > Geometry_StripEfficiency_Den;
  //std::vector<std::vector<std::vector<double> > > Geometry_StripEfficiency;
  

 
 
  NumEventNonCorruptedUtilized->Fill(1,countegood);


  if(produceVfatEffiFile)
    {
      
      std::string namefiledis=outputFileName;
      std::string sought = "root";
      int pos=namefiledis.find(sought); 
      namefiledis.replace(pos, sought.size(), "dat"); 
     	
      std::cout<<" Opening file for store Vfat Efficiency Plot "<<namefiledis.c_str()<<std::endl;

      fstream file_op(namefiledis.c_str(), ios::out); 
      //Warning!! Found very strange problem opening file with "in" flag,file is written only if it already exists in the path.
      //fstream file_op("test.dat",ios::in | ios::out);
      if(((file_op).is_open())==false)
	{
	  std::cout<<" ERROR: Problem in opening file IN T2TrkBased().cc"<<std::endl;
	  (file_op).clear();
	}

      int planeID0=SelectedHalf*10; 
      for(unsigned int m=0;m<10; m++)
	{
	  int planeID=m+planeID0;
	  for(unsigned int iid=0;iid<=16;iid++)
	    {//VFATEFFNormalized
	      file_op<<"DetId: "<<planeID<<"  VFat_iid: "<<iid<<"   Eff: "<< VFATEFF[m]->GetBinContent(iid+1)<<"        pm: "<< VFATEFF[m]->GetBinError(iid+1)<<" ( stat: "<<VFATEFF[m]->GetBinEntries(iid+1)<<" ) \n";
	    }
	  //h.GetBinEntries(i); //return number of entries in bin i
	  //> h.GetBinContent(i); returns the mean value of bin i
	  //> h.GetBinError(i); returns the standard deeviation in bin i

	  //VFATEFF[m]->Write("");
	  //VfatStatistics[m]->Write("");
	}

       file_op.close(); 
    }

 






 
for(unsigned int y=0;y<10;y++)
    {

      if(TotHITNoiseEnt.at(y)>0)
	HITNoiseLevelvsDet->Fill(y,(double)(HITNoiseEnt.at(y))/(double)(TotHITNoiseEnt.at(y)));
      if(TotPADNoiseEnt.at(y)>0)
	PADNoiseLevelvsDet->Fill(y,(double)(PADNoiseEnt.at(y))/(double)(TotPADNoiseEnt.at(y)));
      if(TotSTRIPNoiseEnt.at(y)>0)
	STRIPNoiseLevelvsDet->Fill(y,(double)(STRIPNoiseEnt.at(y))/(double)(TotSTRIPNoiseEnt.at(y)));

    }


  
/*
 Tr4EffDetAll->Fill(Testedcamera.at(0),	Tr4EffDet0->GetBinContent(Testedcamera.at(0)+2));
 //Tr4EffDetAll->SetBinError(2,Tr4EffDet0->GetBinError(2));

 Tr4EffDetAll->Fill(Testedcamera.at(1),	Tr4EffDet1->GetBinContent(Testedcamera.at(1)+2));
 //Tr4EffDetAll->SetBinError(3,Tr4EffDet1->GetBinError(3));
 
 Tr4EffDetAll->Fill(Testedcamera.at(2),	Tr4EffDet2->GetBinContent(Testedcamera.at(2)+2));
 // Tr4EffDetAll->SetBinError(4,Tr4EffDet2->GetBinError(4));

 Tr4EffDetAll->Fill(Testedcamera.at(3),	Tr4EffDet3->GetBinContent(Testedcamera.at(3)+2));
 //Tr4EffDetAll->SetBinError(5,Tr4EffDet3->GetBinError(5));

 Tr4EffDetAll->Fill(Testedcamera.at(4),	Tr4EffDet4->GetBinContent(Testedcamera.at(4)+2));
 //Tr4EffDetAll->SetBinError(6,Tr4EffDet4->GetBinError(6));
 
 Tr4EffDetAllpad->Fill(Testedcamera.at(0),	Tr4EffDet0pad->GetBinContent(Testedcamera.at(0)+2));
 Tr4EffDetAllpad->Fill(Testedcamera.at(1),	Tr4EffDet1pad->GetBinContent(Testedcamera.at(1)+2));
 Tr4EffDetAllpad->Fill(Testedcamera.at(2),	Tr4EffDet2pad->GetBinContent(Testedcamera.at(2)+2));
 Tr4EffDetAllpad->Fill(Testedcamera.at(3),	Tr4EffDet3pad->GetBinContent(Testedcamera.at(3)+2));
 Tr4EffDetAllpad->Fill(Testedcamera.at(4),	Tr4EffDet4pad->GetBinContent(Testedcamera.at(4)+2));

 Tr4EffDetAllstrip->Fill(Testedcamera.at(0),	Tr4EffDet0strip->GetBinContent(Testedcamera.at(0)+2));
 Tr4EffDetAllstrip->Fill(Testedcamera.at(1),	Tr4EffDet1strip->GetBinContent(Testedcamera.at(1)+2));
 Tr4EffDetAllstrip->Fill(Testedcamera.at(2),	Tr4EffDet2strip->GetBinContent(Testedcamera.at(2)+2));
 Tr4EffDetAllstrip->Fill(Testedcamera.at(3),	Tr4EffDet3strip->GetBinContent(Testedcamera.at(3)+2));
 Tr4EffDetAllstrip->Fill(Testedcamera.at(4),	Tr4EffDet4strip->GetBinContent(Testedcamera.at(4)+2));




 //Tr3EffDetAll->Fill(0.,	Tr3EffDet0->GetBinContent(2));
 Tr3EffDetAll->Fill(Testedcamera.at(0),	Tr3EffDet0->GetBinContent(Testedcamera.at(0)+2));
//Tr3EffDetAll->SetBinError(2,	Tr3EffDet0->GetBinError(2));
 Tr3EffDetAll->Fill(Testedcamera.at(1),	Tr3EffDet1->GetBinContent(Testedcamera.at(1)+2));
 //Tr3EffDetAll->SetBinError(3,	Tr3EffDet1->GetBinError(3));
 Tr3EffDetAll->Fill(Testedcamera.at(2),	Tr3EffDet2->GetBinContent(Testedcamera.at(2)+2));
 //Tr3EffDetAll->SetBinError(4,	Tr3EffDet2->GetBinError(4));
 Tr3EffDetAll->Fill(Testedcamera.at(3),	Tr3EffDet3->GetBinContent(Testedcamera.at(3)+2));
 //Tr3EffDetAll->SetBinError(5,	Tr3EffDet3->GetBinError(5));
 Tr3EffDetAll->Fill(Testedcamera.at(4),	Tr3EffDet4->GetBinContent(Testedcamera.at(4)+2));
 //Tr3EffDetAll->SetBinError(6,	Tr3EffDet4->GetBinError(6));




 Tr3EffDetAllpad->Fill(Testedcamera.at(0),	Tr3EffDet0pad->GetBinContent(Testedcamera.at(0)+2));
 Tr3EffDetAllpad->Fill(Testedcamera.at(1),	Tr3EffDet1pad->GetBinContent(Testedcamera.at(1)+2));
 Tr3EffDetAllpad->Fill(Testedcamera.at(2),	Tr3EffDet2pad->GetBinContent(Testedcamera.at(2)+2));
 Tr3EffDetAllpad->Fill(Testedcamera.at(3),	Tr3EffDet3pad->GetBinContent(Testedcamera.at(3)+2));
 Tr3EffDetAllpad->Fill(Testedcamera.at(4),	Tr3EffDet4pad->GetBinContent(Testedcamera.at(4)+2));

 Tr3EffDetAllstrip->Fill(Testedcamera.at(0),	Tr3EffDet0strip->GetBinContent(Testedcamera.at(0)+2));
 Tr3EffDetAllstrip->Fill(Testedcamera.at(1),	Tr3EffDet1strip->GetBinContent(Testedcamera.at(1)+2));
 Tr3EffDetAllstrip->Fill(Testedcamera.at(2),	Tr3EffDet2strip->GetBinContent(Testedcamera.at(2)+2));
 Tr3EffDetAllstrip->Fill(Testedcamera.at(3),	Tr3EffDet3strip->GetBinContent(Testedcamera.at(3)+2));
 Tr3EffDetAllstrip->Fill(Testedcamera.at(4),	Tr3EffDet4strip->GetBinContent(Testedcamera.at(4)+2));

*/




 Tr2EffDetAll->Fill(0.,	Tr2EffDet0->GetBinContent(2));
 Tr2EffDetAll->Fill(1.,	Tr2EffDet1->GetBinContent(3));
 Tr2EffDetAll->Fill(2.,	Tr2EffDet2->GetBinContent(4));
 Tr2EffDetAll->Fill(3.,	Tr2EffDet3->GetBinContent(5));
 Tr2EffDetAll->Fill(4.,	Tr2EffDet4->GetBinContent(6));
 
 Tr2EffDetAllpad->Fill(0.,	Tr2EffDet0pad->GetBinContent(2));
 Tr2EffDetAllpad->Fill(1.,	Tr2EffDet1pad->GetBinContent(3));
 Tr2EffDetAllpad->Fill(2.,	Tr2EffDet2pad->GetBinContent(4));
 Tr2EffDetAllpad->Fill(3.,	Tr2EffDet3pad->GetBinContent(5));
 Tr2EffDetAllpad->Fill(4.,	Tr2EffDet4pad->GetBinContent(6));

 Tr2EffDetAllstrip->Fill(0.,	Tr2EffDet0strip->GetBinContent(2));
 Tr2EffDetAllstrip->Fill(1.,	Tr2EffDet1strip->GetBinContent(3));
 Tr2EffDetAllstrip->Fill(2.,	Tr2EffDet2strip->GetBinContent(4));
 Tr2EffDetAllstrip->Fill(3.,	Tr2EffDet3strip->GetBinContent(5));
 Tr2EffDetAllstrip->Fill(4.,	Tr2EffDet4strip->GetBinContent(6));


 


 /* ----------- Save the output ------------ */
 /* ----------- Save the output ------------ */
 /* ----------- Save the output ------------ */


  
TFile *f = TFile::Open(outputFileName.c_str(), "recreate");
 if( !f || !f->IsWritable() ){
   std::cout << "Output file not opened correctly !!" << std::endl;
 }
 

 gDirectory = f->mkdir("EffiGeometry");

 for(int i=0;i<40;i++){
   HGeometry_StripEfficiency_Num[i]->Write(""); 
   HGeometry_StripEfficiency_Den[i]->Write("");
   HGeometry_PadEfficiency_Num[i]->Write("");
   HGeometry_PadEfficiency_Den[i]->Write("");
 }


 /*------------------------------------------------- ----------- NOISE NOISE ------------ */
 /*------------------------------------------------- ----------- NOISE NOISE ------------ */
 /*------------------------------------------------- ----------- NOISE NOISE ------------ */
                            

 gDirectory = f->mkdir("Noise");
 /*
 for(unsigned int m=0;m<10; m++)
   {
     Class1NoiseDRDPHISep[m]->Write("");
     TrackingHitVsClass1NoiseDRSep[m]->Write("");
     TrackingHitRVsStripCluNoiseDRSep[m]->Write("");
     TrackingHitRVsPadCluNoiseDRSep[m]->Write("");
  
     NoiseHitDet[m]->Write("");
     NoisePadDet[m]->Write("");
     NoiseStripDet[m]->Write("");
     Noise01HitDet[m]->Write("");
     Noise01StripDet[m]->Write("");
     Noise01PadDet[m]->Write("");

   }
 //Noise Histograms  of a selected detector
 PhiNoiseCl1Hit->Write(""); 
 CLSPadNoiseCl1Hit->Write(""); 
 CLSStripNoiseCl1Hit->Write("");
 CLSPhiNoisePad->Write(""); 
 PhiNoisePad->Write(""); 
 RNoiseStrip->Write(""); 
 CLSRNoiseStrip->Write("");
 NoiseEventsStat->Write("");
 RNoiseCl1Hit->Write("");
 DRHit12->Write("");

 Average_Cl1HitNoiseNumb_vsDet->Write(""); 
 STRIPNoiseLevelvsDet->Write(""); 
 PADNoiseLevelvsDet->Write(""); 
 HITNoiseLevelvsDet->Write("");
 */
                                                    
 /*------------------------------------------------- ----------- Resolution ------------ */
 /*------------------------------------------------- ----------- Resolution ------------ */
 /*------------------------------------------------- ----------- Resolution ------------ */

 
 gDirectory = f->mkdir("Resolution");
 /*
 for(unsigned int m=0;m<10; m++)
   {
     DYAlignDet[m]->Write("");
     DXAlignDet[m]->Write("");  
     DRResol[m]->Write("");  
     DPhiResol[m]->Write("");  
     drHit12[m]->Write("");  
     DXResp0[m]->Write("");  
     DYResp0[m]->Write("");  
     DXResp9[m]->Write("");  
     DYResp9[m]->Write("");
   }
     
 AL_DPhiHitNoCutdet0->Write("");
 AL_DRHitNoCutdet0->Write("");   
 AL_DPhiHitNoCutdet1->Write("");
 AL_DRHitNoCutdet1->Write("");
 AL_DPhiHitNoCutdet2->Write("");
 AL_DRHitNoCutdet2->Write("");
 AL_DPhiHitNoCutdet3->Write("");
 AL_DRHitNoCutdet3->Write("");
 AL_DPhiHitNoCutdet4->Write("");
 AL_DRHitNoCutdet4->Write("");
 AL_DPhiHitNoCutdet5->Write("");
 AL_DRHitNoCutdet5->Write("");
 AL_DPhiHitNoCutdet6->Write("");
 AL_DRHitNoCutdet6->Write("");
 AL_DPhiHitNoCutdet7->Write("");
 AL_DRHitNoCutdet7->Write("");
 AL_DPhiHitNoCutdet8->Write("");
 AL_DRHitNoCutdet8->Write("");
 AL_DPhiHitNoCutdet9->Write("");
 AL_DRHitNoCutdet9->Write("");


 */

  gDirectory = f->mkdir("VFAT_efficiency");
  for(unsigned int m=0;m<10; m++)
   {
     VFATEFF[m]->Write("");
     VfatStatistics[m]->Write("");
     VFATEFFNormalized[m]->Write("");
     TrkHitR_vsplane[m]->Write("");
     
   }
  
  for(unsigned int m=0;m<40; m++)
   {
     VfatsCorruptionFrequency[m]->Write("");
     PadClusterSize_AllvsPlane[m]->Write("");
     PadClusterR_AllvsPlane[m]->Write("");  
   } 


  T2vfatinfoErrorMapNOTDEAD->Write("");
  T2vfatinfoErrorMap->SetDrawOption("text");  
  T2vfatinfoErrorMap->Write("");
  T2vfatinfoErrorMapGeoQ1->Write("");
  T2vfatinfoErrorMapFootPr->Write("");
  T2vfatinfoErrorMapCRC->Write("");
  T2vfatinfoErrorNotInData->Write("");
  T2vfatinfoErrorMapNotIncluded->Write("");
  T2vfatinfoErrorErrorMapping->Write("");
  T2vfatinfoErrorMap_DEADWithCorrectFrame->Write("");
  
//T2vfatinfoErrorMapFootPr T2vfatinfoErrorMapCRC T2vfatinfoErrorNotInData T2vfatinfoErrorMapNotIncluded T2vfatinfoErrorMapping
  VfatsCorruptionFrequency2D->Write("");



  
  gDirectory = f->mkdir("VFAT_Monitoring");

  if(VFATMonitoring){

  for(unsigned int m=0;m<40; m++)
    for(unsigned int k=0;k<17; k++)
      {
	vFatCumulative[m][k]->Write();
	vFatMultiplicty[m][k]->Write();
	vFatMultiplictyVsTrk[m][k]->Write();
	TrkMultiplictyWhenEveryThingIsOn[m][k]->Write();
      }
   
  Vfat15_Pl19and17MultCorrel->Write();
  Vfat15_Pl19and11_Pl19MultCorrel->Write();
  Vfat15_Pl19and17MultCorrel_LowMulti->Write();
  //EvtVfat_Strip_IdWhenCompletelyOn_ID->Write();
  EvtVfat_Strip_WhenCompletelyOn->Write();
  NumH1VfatCompletelyOnPerEvt->Write();
  for(unsigned int k=0;k<5; k++){
    Vfat15_Pl10andH1Vfat15_PlDispariMultCorrel[k]->Write();
    Vfat16_Pl10andH1Vfat16_PlDispariMultCorrel[k]->Write();
  }
  
  }  


  
 /*------------------------------------------------- ----------- Global  ------------ */
 /*------------------------------------------------- ----------- Global ------------ */
 /*------------------------------------------------- ----------- Global ------------ */


 gDirectory = f->mkdir("Global");

EventCandidate->Write();
 BunchStatusBit->Write();EventHemisphere->Write();
 Count_t2trackVectorALL_H0->Write();
 Count_t2trackVectorALL_H1->Write();
 Count_t2trackVectorALL_H2->Write();
 Count_t2trackVectorALL_H3->Write();
 
 NumStripCluVsPlaneAll3H0->Write();
 NumPadCluVsPlaneAll3H0->Write();
 NumStripCluVsPlaneAll3H1->Write();
 NumPadCluVsPlaneAll3H1->Write();
 NumStripCluVsPlaneAll3H2->Write();
 NumPadCluVsPlaneAll3H2->Write();
 NumStripCluVsPlaneAll3H3->Write();
 NumPadCluVsPlaneAll3H3->Write();


 NumStripCluVsPlaneAll3H0_Cutted->Write();
 NumStripCluVsPlaneAll3H1_Cutted->Write();
 NumStripCluVsPlaneAll3H2_Cutted->Write();
 NumStripCluVsPlaneAll3H3_Cutted->Write();
 
 NumPadCluVsPlaneAll3H0_Cutted->Write();
 NumPadCluVsPlaneAll3H1_Cutted->Write();
 NumPadCluVsPlaneAll3H2_Cutted->Write();
 NumPadCluVsPlaneAll3H3_Cutted->Write();

 NumPadCluVsPlaneAll3H0_CuttedMinus->Write();
 NumPadCluVsPlaneAll3H1_CuttedMinus->Write();
 NumPadCluVsPlaneAll3H2_CuttedMinus->Write();
 NumPadCluVsPlaneAll3H3_CuttedMinus->Write();
 NumPadCluVsPlaneAll3H0_CuttedPlus->Write();
 NumPadCluVsPlaneAll3H1_CuttedPlus->Write();
 NumPadCluVsPlaneAll3H2_CuttedPlus->Write();
 NumPadCluVsPlaneAll3H3_CuttedPlus->Write();



 NumStripCluVsPlaneAll3H0_CuttedMinus->Write();
 NumStripCluVsPlaneAll3H1_CuttedMinus->Write();
 NumStripCluVsPlaneAll3H2_CuttedMinus->Write();
 NumStripCluVsPlaneAll3H3_CuttedMinus->Write();
 NumStripCluVsPlaneAll3H0_CuttedPlus->Write();
 NumStripCluVsPlaneAll3H1_CuttedPlus->Write();
 NumStripCluVsPlaneAll3H2_CuttedPlus->Write();
 NumStripCluVsPlaneAll3H3_CuttedPlus->Write();
 
 


 NumPadCluVsPlaneAll3H0_Cutted_LOWMultipl->Write();
 NumPadCluVsPlaneAll3H1_Cutted_LOWMultipl->Write();
 NumPadCluVsPlaneAll3H2_Cutted_LOWMultipl->Write();
 NumPadCluVsPlaneAll3H3_Cutted_LOWMultipl->Write();

 CumulativeNumPadCluAll3H0->Write();
 CumulativePadCluSizeAll3H0->Write();
 CumulativeNumStripCluAll3H0->Write();
 CumulativeStripCluSizeAll3H0->Write();
 CumulativePadCluSize_UsedInEffi_HX->Write();
 CumulativeStripCluSize_UsedInEffi_HX->Write();
   
 StripCluSizeVsPlaneAll3H0->Write();
 PadCluSizeVsPlaneAll3H0->Write();
 StripCluSizeVsPlaneAll3H1->Write();
 PadCluSizeVsPlaneAll3H1->Write();
 StripCluSizeVsPlaneAll3H2->Write();
 PadCluSizeVsPlaneAll3H2->Write();
 StripCluSizeVsPlaneAll3H3->Write();
 PadCluSizeVsPlaneAll3H3->Write();


PadCluMultiplVsPlane->Write();
 NumPadCluVsPlane->Write(""); 
 NumStripCluVsPlane->Write(""); 
 NumCl1HitVsPlane->Write(""); 
 PadCluSizeVsPlaneAll->Write(""); 
 StipCluSizeVsPlaneAll->Write("");
 PadCluSizeVsPlaneAll2->Write(""); 

 Class1HitPadStripCLSCorrel->Write("");
 DigiStripOccupancy->Write(""); 
 DigiPadOccupancy->Write("");

 PlaneNumPadClu11_vs_16->Write("");
 for(unsigned int m=0;m<40; m++)
   {
     PadClusterR[m]->Write("");     
     PadClusterXYRadiographyYPlus[m]->Write("");
     PadClusterXYRadiographyYMinus[m]->Write("");
     PlaneHitRadiographyXY[m]->Write("");
     StripCluRadiographyXY[m]->Write("");
     if((m%10)==0)
       {
	 unsigned int hindex=m/10;	 
	 HalfTeleTrkRadiographyXY[hindex]->Write("");
       }

     if(m<10)
       PadClusterXYForEffiEvents[m]->Write("");
     

     //StripClusterXYRadiographyYPlus[m]->Write("");
     //StripClusterXYRadiographyYMinus[m]->Write("");
   }

 
 
 /*
 for(unsigned int m=0;m<10; m++)
   {
     RPhiEvtHit[m]->Write("");
      if(simufile){
	PAdNr[m]->Write(""); 
	StripNr[m]->Write("");      
      }
   } 
 */

 NumTrackALL_EveryQuarter->Write(""); 
 NumTrackALL_ONEQuarter->Write("");

 TrkQuarterIdALL->Write("");
 NumhitinTrackALL->Write("");
 TrkphiALL->Write("");
 TrketaALL->Write("");


 NumhitinTrack->Write("");
 HitMatrixR->Write(""); 


 CluStripentriesGoodHit->Write("");
 CluPadentriesGoodHit->Write("");
 TrackPhi0R0->Write("");
 NumhitinTrackGood->Write("");
 chiRredREF->Write("");
 ActiveSymbDet->Write("");
 Chi2phired->Write("");
 PolarAngles->Write("");
 
 clusterstripentries->Write("");
 clusterpadentries->Write("");

 Chi2red->Write("");

 Tr4EffDetAll->Write("");
 IstoDisp->Write("");
 AllHitR->Write("");
AllHitZ->Write("");
AllHitZCL1->Write(""); 
AllHitRCL1->Write(""); 
AllHitZStripOnly->Write(""); 
AllHitRStripOnly->Write(""); 
AllHitZPadOnly->Write(""); 
AllHitRPadOnly->Write("");

 drvsIddet->Write(""); 
 alldetdr->Write(""); 
 dphivsIddet->Write(""); 
 alldetdphi->Write("");
 
 DetIdvsYShOnlyYfit->Write("");
 DetIdvsXShOnlyXfit->Write("");

 AllDetId->Write("");
 AllDetIdPlaneS->Write("");
 AllDetIdPlane->Write("");
 C12AllDetIdPlaneS->Write("");
 C12AllDetIdPlane->Write("");
 
 Trketagood->Write("");
 Trkphigood->Write("");
 Trketa->Write("");
 Trkphi->Write("");
 TrkphiRZPlus->Write("");
 TrkphiRZMinus->Write("");



 Chi2PhiProb->Write("");
 Chi2RProb->Write(""); 
 
 NumAllHitvsIdAllevt->Write("");
 NumAllHitvsIdGoodevt->Write("");
 
 NumhitinTrackAligned->Write("");

 CumulativeExpXUncert->Write("");
 CumulativeExpYUncert->Write("");
 CumulativePhiResol->Write("");
 CumulativeRResol->Write("");
 RandomGauss->Write("");
 TrkChi2X->Write("");
 TrkChi2Y->Write("");
 ProbXHisto->Write("");
 ProbYHisto->Write("");
 AXHisto->Write("");
 AYHisto->Write("");
 AXError->Write("");
 AYError->Write("");
 BXError->Write("");
 BYError->Write("");
 HitMatrixPhi->Write("");
 
 ProbXHistoAl->Write("");
 ProbYHistoAl->Write("");
 HitXErr->Write("");
 HitYErr->Write("");
 HitX->Write("");
 HitY->Write("");
 Plside0PSimHitXD0->Write("");
 Plside1PSimHitXD0->Write("");
 PSimHitXD1->Write("");
 Plside0PSimHitYD0->Write("");
 Plside1PSimHitYD0->Write("");
 PSimHitYD1->Write("");

 C11Hit_ClsPadVsClsPadCol->Write("");
 C11Hit_ClsPadVsClsStrip->Write("");
 C1Hit_DphiVsPadCLS->Write("");
 C11Hit_ClsPadVsR->Write("");
 Cl1Hit_ClsStripVsR->Write("");
 C1Hit_Dphi->Write("");
 C1Hit_PadClsVsPadNumRow->Write("");  
 C1Hit_PadClsVsPadNumCol->Write("");  
 C1Hit_Pad_PadNumRowVsPadNumCol->Write("");  





 /*------------------------------------------------- ----------- Efficiency  ------------ */
 /*------------------------------------------------- ----------- Efficiency ------------ */
 /*------------------------------------------------- ----------- Efficiency ------------ */



 gDirectory = f->mkdir("Efficiency");

 
 NumEventNonCorruptedUtilized->Write("");
 for(unsigned int m=0;m<4; m++)
   {    
     EFFpadVsIdNormalized_sectorX[m]->Write("");    
     EFFstripVsIdNormalized_sectorX[m]->Write("");
   }

 
 Tracksplanestat3->Write("");
 Tracksplanestat4->Write("");
 TracksplanestatN->Write("");
 TracksplanestatNm1->Write("");
 TrnEffDetAllpad->Write("");
 Trnm1EffDetAllpad->Write("");
 TrnEffDetAllstrip->Write("");
 Trnm1EffDetAllstrip->Write("");
 
 TrnEffDetAll->Write("");
 Trnm1EffDetAll->Write("");
 
 Tr4EffDetAllstrip->Write("");
 
 Tr4EffDetAllpad->Write("");
 
 Tr3EffDetAll->Write("");
 Tr3EffDetAllstrip->Write("");
 Tr3EffDetAllpad->Write("");
 
 Tr2EffDetAll->Write("");
 Tr2EffDetAllstrip->Write("");
 Tr2EffDetAllpad->Write("");

 
 EFFHitVsIdNormalized->Write("");
 EFFstripVsIdNormalized->Write("");  
 EFFpadVsIdNormalized->Write("");



 EFFstripVsIdNoCut->Write("");
 EFFstripVsId->Write("");
 EFFpadVsIdNoCut->Write("");
 EFFpadVsId->Write("");
 EFFHitVsIdNoCut->Write("");
 EFFHitVsId->Write("");
 EFFClusterVsId->Write("");
 
 EFF_DRstripNoCutdet0->Write("");
 EFF_DRstripdet0->Write("");
 EFF_DRHitNoCutdet0->Write("");
 EFF_DRHitdet0->Write("");
 
 EFF_DRstripNoCutdet2->Write("");
 EFF_DRstripdet2->Write("");
 EFF_DRHitNoCutdet2->Write("");
 EFF_DRHitdet2->Write("");
   
 EFF_DRstripNoCutdet4->Write("");
 EFF_DRstripdet4->Write("");
 EFF_DRHitNoCutdet4->Write("");
 EFF_DRHitdet4->Write("");

 EFF_DRstripNoCutdet6->Write("");
 EFF_DRstripdet6->Write("");
 EFF_DRHitNoCutdet6->Write("");
 EFF_DRHitdet6->Write("");
 
 EFF_DRstripNoCutdet8->Write("");
 EFF_DRstripdet8->Write("");
 EFF_DRHitNoCutdet8->Write("");
 EFF_DRHitdet8->Write("");
 
 //--
 EFF_DRstripNoCutdet1->Write("");
 EFF_DRstripdet1->Write("");
 EFF_DRHitNoCutdet1->Write("");
 EFF_DRHitdet1->Write("");
     
    
 EFF_DRstripNoCutdet3->Write("");
 EFF_DRstripdet3->Write("");
 EFF_DRHitNoCutdet3->Write("");
 EFF_DRHitdet3->Write("");
 
   
 EFF_DRstripNoCutdet5->Write("");
 EFF_DRstripdet5->Write("");
 EFF_DRHitNoCutdet5->Write("");
 EFF_DRHitdet5->Write("");

 EFF_DRstripNoCutdet7->Write("");
 EFF_DRstripdet7->Write("");
 EFF_DRHitNoCutdet7->Write("");
 EFF_DRHitdet7->Write("");
 
 EFF_DRstripNoCutdet9->Write("");
 EFF_DRstripdet9->Write("");
 EFF_DRHitNoCutdet9->Write("");
 EFF_DRHitdet9->Write("");
   //--


   EFF_DPhipadNoCutdet0->Write("");
   EFF_DPhipaddet0->Write("");
   EFF_DPhiHitNoCutdet0->Write("");
   EFF_DPhiHitdet0->Write("");

   EFF_DPhipadNoCutdet2->Write("");
   EFF_DPhipaddet2->Write("");
   EFF_DPhiHitNoCutdet2->Write("");
   EFF_DPhiHitdet2->Write("");

   EFF_DPhipadNoCutdet4->Write("");
   EFF_DPhipaddet4->Write("");
   EFF_DPhiHitNoCutdet4->Write("");
   EFF_DPhiHitdet4->Write("");

   EFF_DPhipadNoCutdet6->Write("");
   EFF_DPhipaddet6->Write("");
   EFF_DPhiHitNoCutdet6->Write("");
   EFF_DPhiHitdet6->Write("");

   EFF_DPhipadNoCutdet8->Write("");
   EFF_DPhipaddet8->Write("");
   EFF_DPhiHitNoCutdet8->Write("");
   EFF_DPhiHitdet8->Write("");

   //--
   EFF_DPhipadNoCutdet1->Write("");
   EFF_DPhipaddet1->Write("");
   EFF_DPhiHitNoCutdet1->Write("");
   EFF_DPhiHitdet1->Write("");

   EFF_DPhipadNoCutdet3->Write("");
   EFF_DPhipaddet3->Write("");
   EFF_DPhiHitNoCutdet3->Write("");
   EFF_DPhiHitdet3->Write("");

   EFF_DPhipadNoCutdet5->Write("");
   EFF_DPhipaddet5->Write("");
   EFF_DPhiHitNoCutdet5->Write("");
   EFF_DPhiHitdet5->Write("");

   EFF_DPhipadNoCutdet7->Write("");
   EFF_DPhipaddet7->Write("");
   EFF_DPhiHitNoCutdet7->Write("");
   EFF_DPhiHitdet7->Write("");

   EFF_DPhipadNoCutdet9->Write("");
   EFF_DPhipaddet9->Write("");
   EFF_DPhiHitNoCutdet9->Write("");
   EFF_DPhiHitdet9->Write("");
   //--

   EFF_CLSStripdet0->Write("");
   EFF_CLSPaddet0->Write("");

   EFF_CLSStripdet2->Write("");
   EFF_CLSPaddet2->Write("");

   EFF_CLSStripdet4->Write("");
   EFF_CLSPaddet4->Write("");

   EFF_CLSStripdet6->Write("");
   EFF_CLSPaddet6->Write("");

   EFF_CLSStripdet8->Write("");
   EFF_CLSPaddet8->Write("");
   EFF_CLSStripdet0NoCut->Write("");
   EFF_CLSPaddet0NoCut->Write("");

   EFF_CLSStripdet2NoCut->Write("");
   EFF_CLSPaddet2NoCut->Write("");

   EFF_CLSStripdet4NoCut->Write("");
   EFF_CLSPaddet4NoCut->Write("");

   EFF_CLSStripdet6NoCut->Write("");
   EFF_CLSPaddet6NoCut->Write("");

   EFF_CLSStripdet8NoCut->Write("");
   EFF_CLSPaddet8NoCut->Write("");
     
   //--
   EFF_CLSStripdet1->Write("");
   EFF_CLSPaddet1->Write("");

   EFF_CLSStripdet3->Write("");
   EFF_CLSPaddet3->Write("");

   EFF_CLSStripdet5->Write("");
   EFF_CLSPaddet5->Write("");

   EFF_CLSStripdet7->Write("");
   EFF_CLSPaddet7->Write("");

   EFF_CLSStripdet9->Write("");
   EFF_CLSPaddet9->Write("");


   EFF_CLSStripdet1NoCut->Write("");
   EFF_CLSPaddet1NoCut->Write("");

   EFF_CLSStripdet3NoCut->Write("");
   EFF_CLSPaddet3NoCut->Write("");

   EFF_CLSStripdet5NoCut->Write("");
   EFF_CLSPaddet5NoCut->Write("");

   EFF_CLSStripdet7NoCut->Write("");
   EFF_CLSPaddet7NoCut->Write("");

   EFF_CLSStripdet9NoCut->Write("");
   EFF_CLSPaddet9NoCut->Write("");

   EFF_C1CLSStripdet0->Write("");
   EFF_C1CLSPaddet0->Write("");
   
   EFF_C1CLSStripdet1->Write("");
   EFF_C1CLSPaddet1->Write("");

   EFF_C1CLSStripdet2->Write("");
   EFF_C1CLSPaddet2->Write("");

   EFF_C1CLSStripdet3->Write("");
   EFF_C1CLSPaddet3->Write("");

   EFF_C1CLSStripdet4->Write("");
   EFF_C1CLSPaddet4->Write("");

   EFF_C1CLSStripdet5->Write("");
   EFF_C1CLSPaddet5->Write("");

   EFF_C1CLSStripdet6->Write("");
   EFF_C1CLSPaddet6->Write("");

   EFF_C1CLSStripdet7->Write("");
   EFF_C1CLSPaddet7->Write("");

   EFF_C1CLSStripdet8->Write("");
   EFF_C1CLSPaddet8->Write("");

   EFF_C1CLSStripdet9->Write("");
   EFF_C1CLSPaddet9->Write("");

   RelEffi_ForTracking->Write("");
   
   
   
    
   


  f->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(T2AnalyzerRaw);

// -*- C++ -*-
//
// Original Author:  Mirko Berretti
//         Created:  Wed Dec 12 18:19:44 CET 2010
// $Id: T2HalfQuarterTrkEfficiency.h,v 1.1.2.1 2009/11/08 14:50:48 berretti Exp $
//
//T2HalfQuarterTrkEfficiency.h
#ifndef _TotemAnalysis_T2HalfQuarterTrkEfficiency_T2HalfQuarterTrkEfficiency_H_
#define _TotemAnalysis_T2HalfQuarterTrkEfficiency_T2HalfQuarterTrkEfficiency_H_

// system include files
#include <memory>
#include <math.h>
//#include <TMath.h>
// user include files
#include "TotemT1T2Validation/T2GeometryValidation/interface/DAQInformationSourceXML_a.h"
#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/T2Cluster/interface/T2Cluster.h"
#include "DataFormats/T2Hit/interface/T2Hit.h"
#include "DataFormats/T2Hit/interface/T2HitCollection.h"
#include "DataFormats/T2Cluster/interface/T2PadClusterCollection.h"
#include "DataFormats/T2Cluster/interface/T2StripClusterCollection.h"
#include "DataFormats/T2Digi/interface/T2PadDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2PadDigi.h"
#include "DataFormats/T2Digi/interface/T2StripDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2StripDigi.h"
#include "TotemAnalysis/T2Cuts/interface/T2SelectionCutUtils.h"

#include "DataFormats/T2DetId/interface/T2DetId.h"
#include "DataFormats/T2Road/interface/T2Road.h"
#include "DataFormats/T2Road/interface/T2RoadCollection.h"
#include "DataFormats/T1T2Track/interface/T1T2Track.h"
#include "DataFormats/T1T2Track/interface/T1T2TrackCollection.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "Geometry/TotemGeometry/interface/T2Geometry.h"
#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"
#include <boost/shared_ptr.hpp>
#include "TotemRawDataLibrary/DataFormats/interface/RawEvent.h"
#include "DataFormats/T2DigiVfat/interface/T2VfatInformation.h"


//#include "Geometry/TotemGeometry/interface/T1Geometry.h"
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include "TCanvas.h"
#include "TStyle.h"
#include "TProfile.h"

//#include "FWCore/Framework/interface/GeneratedInputSource.h"


// class decleration
//

class T2HalfQuarterTrkEfficiency : public edm::EDAnalyzer {
public:
  explicit T2HalfQuarterTrkEfficiency(const edm::ParameterSet&);
  ~T2HalfQuarterTrkEfficiency();
  virtual void beginJob();

  virtual void analyze(const edm::Event&,
		       const edm::EventSetup&
		       );

  virtual void endJob(); 
  DAQInformationSourceXML_a* onemap;
  boost::shared_ptr<DAQInformationT2_a> testmap;

  DAQInformationSourceXML_a* MapProducer_FullVFATs;
  boost::shared_ptr<DAQInformationT2_a> Map_FullVfats;
  

  DAQInformationSourceXML_a* MapProducer_NotDeadUsedVfat;
  boost::shared_ptr<DAQInformationT2_a> Map_NotDeadUsedVfat;

  std::vector<double> MyLinearfitX(std::vector<T2Hit> hitvec,unsigned int UseJointProb);
std::vector<double> MyLinearfitY(std::vector<T2Hit> hitvec,unsigned int UseJointProb);

  std::vector<double> TrkPredictedPoint(T1T2Track &thetrk, double planeZ);
  bool RecoTrk_UnfoldingCut_56Division(T1T2Track &atrk);
  bool RecoTrk_UnfoldingCut(T1T2Track &atrk, double &Zimp);
  unsigned int RawtoSymb(uint32_t thedet);
  bool IsVfatMapped(int symbvfat);
  // virtual void beginJob(const edm::EventSetup&);
  
  std::vector<int> AllCorruptedEvents;


  std::vector<bool> vfatsCorrupted;
  std::vector<int> vfatsSymbId;
  
  double Trk_RSeparation;
  double Trk_PhiSeparation;
  double Trk_EtaSeparation;
  
  void CalculateEfficiency(std::vector<T1T2Track>* CleanRefArmXTrk ,std::vector<T1T2Track>* TestArmXTrk,unsigned int refQuarter,unsigned int tested,std::string Up_or_DownY,edm::Handle<T2HitCollection> t2hitcoll_,int numreftrks_raw, double avgHitMultRef, double &investigatethisevent);

  double GetMiniumumPhiDistance(T1T2Track reftrk,T1T2TrackCollection::const_iterator thebegin_, T1T2TrackCollection::const_iterator theend_);

  // void FillBPShadow(edm::Handle<T1T2TrackCollection> trackColl , const T2PadDigiCollection* PadDigiptr, const T2StripDigiCollection* StripDigiptr); //old (2010) def.
  void FillBPShadow(edm::Handle<T1T2TrackCollection> trackColl , const T2PadDigiCollection* PadDigiptr, const T2StripDigiCollection* StripDigiptr,int quarter);


  bool QuarterBombarda(int numplane, int multipadplane, int multstripplane, const T2PadDigiCollection* PadDigiptr, const T2StripDigiCollection* StripDigiptr);

  //  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //virtual void endJob() ;

 
  std::vector<T2Hit> TrkInRightQuarter(T1T2Track &trk, int quarter,bool &rightq);

private:

  edm::InputTag t2PadDigiCollectionLabel;  
  edm::InputTag t2StripDigiCollectionLabel;
  edm::InputTag t2VfatInformationLabel;
  //  std::vector<TPolyLine3D*> All3dtracks;
  //std::vector<TPolyLine3D*> tracks0123;
  //std::vector<TPolyLine3D*> tracks0124;
  //std::vector<TPolyLine3D*> tracks0234;
  //std::vector<TPolyLine3D*> tracks1234;
  //std::vector<TPolyLine3D*> tracks0134;

  std::vector<bool> ExludeThisNoisyplaneForEffMeas(T2StripClusterCollection::const_iterator itstrip, T2PadClusterCollection::const_iterator itpad);
  
  bool PlaneInRightHalf(uint32_t cmsswid);

  bool HitInOverlappingRegion(T2HitCollection::const_iterator ithit);

  // bool MettiTracciaBuona(std::vector<TrackInfo>* matricetracce,T1T2TrackCollection::const_iterator  trk, unsigned int numevento);

  bool HitIndetector(std::vector<double> vpar, unsigned int symbdetid,std::map<unsigned int,double> mapsymbtoZ);


  bool CloseTrksInEvt(T1T2TrackCollection trackColl);

  T2Hit GiveExtrapolatedHit(int plane0_10, std::vector<T2Hit> referenceHits);



  
  bool HitIsInTrackColl(T2Hit trackinghit,T1T2TrackCollection trackColl);

  std::vector<double> ResiduiForStrip(std::vector<T2Hit>refhitv, T2Hit hit);
  std::vector<double> ResiduiForPad(std::vector<T2Hit>refhitv, T2Hit hit);
  std::vector<double> ResiduiForC1HIT(std::vector<T2Hit>refhitv, T2Hit hit);

  bool IsEventCorrupted(const T2VfatInformation* t2vfatinfoptr);
 
  bool CheckEventCorruption();  
  bool TrackInOverlappingRegion(T1T2Track trk);
  unsigned int Howmanyplanes(std::vector<unsigned int> *total_cl1Hit_CloseInArm0Up,unsigned int testedQuarter);

 
  int countegood;
  std::auto_ptr<TH1F> NumhitinTrack;
  std::auto_ptr<TH1F> NumhitinTrackGood;


//Noise Histograms  of a selected detector

std::auto_ptr<TProfile> Average_Cl1HitNoiseNumb_vsDet;
std::auto_ptr<TProfile>   DigiPadOccupancy;
std::auto_ptr<TProfile>   DigiStripOccupancy;
std::auto_ptr<TH1F> PhiNoiseCl1Hit;
std::auto_ptr<TH1F> CLSPadNoiseCl1Hit;
std::auto_ptr<TH1F> CLSStripNoiseCl1Hit;
std::auto_ptr<TH1F> CLSPhiNoisePad;
std::auto_ptr<TH1F> PhiNoisePad;
std::auto_ptr<TH1F> RNoiseStrip; 
std::auto_ptr<TH1F> CLSRNoiseStrip;
std::auto_ptr<TH1F> NoiseEventsStat;
std::auto_ptr<TH1F> RNoiseCl1Hit;
std::auto_ptr<TH2F> C11Hit_ClsPadVsClsPadCol;
std::auto_ptr<TH2F> C11Hit_ClsPadVsClsStrip;
std::auto_ptr<TH2F> C1Hit_DphiVsPadCLS;
std::auto_ptr<TH2F> C1Hit_PadClsVsPadNumRow;
std::auto_ptr<TH2F> C1Hit_PadClsVsPadNumCol;
std::auto_ptr<TH2F> C1Hit_Pad_PadNumRowVsPadNumCol;

std::auto_ptr<TProfile> C11Hit_ClsPadVsR;
std::auto_ptr<TProfile> Cl1Hit_ClsStripVsR;
std::auto_ptr<TH1F> C1Hit_Dphi;
std::auto_ptr<TH1F> DRTrk_HRefHTest;

  std::auto_ptr<TH2D> Class1HitPadStripCLSCorrel;

  std::auto_ptr<TH1F> HitZ;
  std::auto_ptr<TH1F> clusterstripentries;
  std::auto_ptr<TH1F> clusterpadentries;
  std::auto_ptr<TH1F> diffphiCluGun;
  std::auto_ptr<TH1F> RLocHRecoH;
  std::auto_ptr<TH1F> Trketa;
  std::auto_ptr<TH1F> Trkphi;
  
  std::auto_ptr<TH1F> RelEffi_ForTracking;

  std::auto_ptr<TH1F> TrkphiALL;
  std::auto_ptr<TH1F> TrketaALL;
  std::auto_ptr<TH1F> TrketaCandidateForShad;
  std::auto_ptr<TH1F> TrketaInWantedQuarterProb; 
  std::auto_ptr<TH1F> TrketaInWantedQuarter;

  std::auto_ptr<TH1F>  Num_RefHalfTrkH0;
  std::auto_ptr<TH1F>  Num_RefHalfTrkH1;
  std::auto_ptr<TH1F>  Num_RefHalfTrkH2;
  std::auto_ptr<TH1F>  Num_RefHalfTrkH3;

  std::auto_ptr<TH1F> NumhitinTrackALL;
  std::auto_ptr<TH1F> TrkQuarterIdALL;
  std::auto_ptr<TH1F> Reference_RawTrkNumHit_NoCut;

  std::auto_ptr<TH1F> ProbXHisto;
  std::auto_ptr<TH1F> ProbYHisto;
 std::auto_ptr<TH1F> HitRForXShad;

  std::auto_ptr<TH1F> DPhiGoodTrk;
  std::auto_ptr<TH1F> DEtaGoodTrk;
  std::auto_ptr<TH1F> diffphiGUNHIT;
  std::auto_ptr<TH1F> diffphiCluHit;
  std::auto_ptr<TH1F> Trkphigood;
  std::auto_ptr<TH1F> Trketagood;
  std::auto_ptr<TH1F> diffRCluHit;
  std::auto_ptr<TH1F> Chi2RProb;
  std::auto_ptr<TH1F> Chi2PhiProb;
  std::auto_ptr<TH1F> chiRredREF;
  std::auto_ptr<TH1F> HitXErr;
   std::auto_ptr<TH1F> HitYErr;
  std::auto_ptr<TH1F> HitX;
   std::auto_ptr<TH1F> HitY;

  std::auto_ptr<TH1F> NumTrackALL_ONEQuarter;
  std::auto_ptr<TH1F> NumTrackALL_EveryQuarter;

  std::auto_ptr<TH1F> IstoDisp;
  std::auto_ptr<TProfile> TrnEffDetAllpad;
  std::auto_ptr<TProfile> Trnm1EffDetAllpad;
  std::auto_ptr<TProfile> TrnEffDetAllstrip;
  std::auto_ptr<TProfile> Trnm1EffDetAllstrip;

  std::auto_ptr<TProfile>  SingleParticleEfficiency;
  std::auto_ptr<TProfile>  SingleParticleEfficiencyCut;
  std::auto_ptr<TProfile>  MultiParticleEfficiencyCut;
  std::auto_ptr<TProfile>  MultiParticleEfficiencyCutNorm;


  std::auto_ptr<TProfile> Tr4EffDet0;
  std::auto_ptr<TProfile> Tr4EffDet1;
  std::auto_ptr<TProfile> Tr4EffDet2;
  std::auto_ptr<TProfile> Tr4EffDet3;
  std::auto_ptr<TProfile> Tr4EffDet4;
  std::auto_ptr<TProfile> drvsIddet;
  std::auto_ptr<TH1F> alldetdr; 
  std::auto_ptr<TProfile> dphivsIddet; 
  std::auto_ptr<TH1F> alldetdphi; 
  std::auto_ptr<TH2F> TrackPhi0R0;
  

  std::auto_ptr<TProfile>  DRHit12;

  

  
  std::auto_ptr<TH1F> ProbYHistoAl;
  std::auto_ptr<TH1F> ProbXHistoAl;

  std::auto_ptr<TH1F>  AXError;
  std::auto_ptr<TH1F>   AYError;
  std::auto_ptr<TH1F>   BXError;
  std::auto_ptr<TH1F>    BYError;
  std::auto_ptr<TH1F>     AXHisto;
  std::auto_ptr<TH1F>     AYHisto;


  
  std::auto_ptr<TH1F>   NumhitinTrackAligned;
 
  
  std::auto_ptr<TH1F> AL_DPhiHitNoCutdet0;
  std::auto_ptr<TH1F> AL_DPhiHitNoCutdet1;
  std::auto_ptr<TH1F> AL_DPhiHitNoCutdet2;
  std::auto_ptr<TH1F> AL_DPhiHitNoCutdet3;
  std::auto_ptr<TH1F> AL_DPhiHitNoCutdet4;
  std::auto_ptr<TH1F> AL_DPhiHitNoCutdet5;
  std::auto_ptr<TH1F> AL_DPhiHitNoCutdet6;
  std::auto_ptr<TH1F> AL_DPhiHitNoCutdet7;
  std::auto_ptr<TH1F> AL_DPhiHitNoCutdet8;
  std::auto_ptr<TH1F> AL_DPhiHitNoCutdet9;

  std::auto_ptr<TH1F> AL_DRHitNoCutdet0;
  std::auto_ptr<TH1F> AL_DRHitNoCutdet1;
  std::auto_ptr<TH1F> AL_DRHitNoCutdet2;
  std::auto_ptr<TH1F> AL_DRHitNoCutdet3;
  std::auto_ptr<TH1F> AL_DRHitNoCutdet4;
  std::auto_ptr<TH1F> AL_DRHitNoCutdet5;
  std::auto_ptr<TH1F> AL_DRHitNoCutdet6;
  std::auto_ptr<TH1F> AL_DRHitNoCutdet7;
  std::auto_ptr<TH1F> AL_DRHitNoCutdet8;
  std::auto_ptr<TH1F> AL_DRHitNoCutdet9;
  


  std::auto_ptr<TH1F>   CluPadentriesGoodHit; 
  std::auto_ptr<TH1F>   CluStripentriesGoodHit;
  std::auto_ptr<TH1F> HITNoiseLevelvsDet;
  std::auto_ptr<TH1F> PADNoiseLevelvsDet;
  std::auto_ptr<TH1F> STRIPNoiseLevelvsDet;
  std::vector<unsigned int>  HITNoiseEnt;
  std::vector<unsigned int>  PADNoiseEnt;
  std::vector<unsigned int>  STRIPNoiseEnt;

  std::vector<unsigned int>  TotHITNoiseEnt;
  std::vector<unsigned int>  TotPADNoiseEnt;
  std::vector<unsigned int>  TotSTRIPNoiseEnt;


  
  std::auto_ptr<TProfile> Tr4EffDet0strip;
  std::auto_ptr<TProfile> Tr4EffDet1strip;
  std::auto_ptr<TProfile> Tr4EffDet2strip;
  std::auto_ptr<TProfile> Tr4EffDet3strip;
  std::auto_ptr<TProfile> Tr4EffDet4strip;
  std::auto_ptr<TProfile> Tr4EffDetAllstrip;
  std::auto_ptr<TProfile> Tr4EffDet0pad;
  std::auto_ptr<TProfile> Tr4EffDet1pad;
  std::auto_ptr<TProfile> Tr4EffDet2pad;
  std::auto_ptr<TProfile> Tr4EffDet3pad;
  std::auto_ptr<TProfile> Tr4EffDet4pad;
  std::auto_ptr<TProfile> Tr4EffDetAllpad;

  std::auto_ptr<TProfile> Trnm1EffDetAll;
  std::auto_ptr<TProfile> Tr4EffDetAll;
  std::auto_ptr<TProfile> TrnEffDetAll;


  
  std::auto_ptr<TProfile> Tr3EffDet0;
  std::auto_ptr<TProfile> Tr3EffDet1;
  std::auto_ptr<TProfile> Tr3EffDet2;
  std::auto_ptr<TProfile> Tr3EffDet3;
  std::auto_ptr<TProfile> Tr3EffDet4;
  std::auto_ptr<TProfile> Tr3EffDetAll;

  std::auto_ptr<TProfile> Tr3EffDet0strip;
  std::auto_ptr<TProfile> Tr3EffDet1strip;
  std::auto_ptr<TProfile> Tr3EffDet2strip;
  std::auto_ptr<TProfile> Tr3EffDet3strip;
  std::auto_ptr<TProfile> Tr3EffDet4strip;
  std::auto_ptr<TProfile> Tr3EffDetAllstrip;
  std::auto_ptr<TProfile> Tr3EffDet0pad;
  std::auto_ptr<TProfile> Tr3EffDet1pad;
  std::auto_ptr<TProfile> Tr3EffDet2pad;
  std::auto_ptr<TProfile> Tr3EffDet3pad;
  std::auto_ptr<TProfile> Tr3EffDet4pad;
  std::auto_ptr<TProfile> Tr3EffDetAllpad;

  std::auto_ptr<TH2F>  RPhiEvtHit[10];
  std::auto_ptr<TH1F> drHit12[10];
  std::auto_ptr<TH1F>  Noise01HitDet[10];
  std::auto_ptr<TH1F>  NoiseHitDet[10];
  std::auto_ptr<TH1F>  Noise01PadDet[10];
  std::auto_ptr<TH1F>  NoisePadDet[10];
  std::auto_ptr<TH1F>  Noise01StripDet[10];
  std::auto_ptr<TH1F>  NoiseStripDet[10];
  std::auto_ptr<TH1F> DYAlignDet[10];
  std::auto_ptr<TH1F> DXAlignDet[10];
  std::auto_ptr<TH1F>  DRResol[10];
  std::auto_ptr<TH1F>  DPhiResol[10];

  std::auto_ptr<TH1F>  DXResp0[10];
  std::auto_ptr<TH1F>  DYResp0[10];
  std::auto_ptr<TH1F>  DXResp9[10];
  std::auto_ptr<TH1F>  DYResp9[10];

  //Noise correlation Histogr
  std::auto_ptr<TH2F> Class1NoiseDRDPHISep[10];  
  std::auto_ptr<TProfile> TrackingHitVsClass1NoiseDRSep[10];  
  std::auto_ptr<TProfile> TrackingHitRVsStripCluNoiseDRSep[10];  
  std::auto_ptr<TProfile> TrackingHitRVsPadCluNoiseDRSep[10];

  std::auto_ptr<TH1D>  BP_XProjection_HX_singlePlaneFromTrk[10];
  std::auto_ptr<TH1D>  BP_YProjection_HX_singlePlaneFromTrk[10];
  std::auto_ptr<TH1D>  BP_YPositiveProjection_HX_singlePlaneFromTrk[10];
  std::auto_ptr<TH1D>  BP_YNegativeProjection_HX_singlePlaneFromTrk[10];
  std::auto_ptr<TH1D>  BP_XYProjection_HX_singlePlaneZValue[10]; 

  std::auto_ptr<TH1F> PAdNr[10];
  std::auto_ptr<TH1F> StripNr[10];
  std::auto_ptr<TProfile> VFATEFF[10];
  std::auto_ptr<TProfile> VFATEFFNormalized[10];
  std::auto_ptr<TH1F> VfatStatistics[10];
  

  std::auto_ptr<TH2D> T2vfatinfoErrorMap;
  std::auto_ptr<TH2D> T2vfatinfoErrorMapNOTDEAD;
  std::auto_ptr<TH2D> T2vfatinfoErrorMapFootPr;
  std::auto_ptr<TH2D> T2vfatinfoErrorMapCRC;
  std::auto_ptr<TH2D> T2vfatinfoErrorNotInData;
  std::auto_ptr<TH2D> T2vfatinfoErrorMapNotIncluded;
  std::auto_ptr<TH2D> T2vfatinfoErrorErrorMapping;

  std::auto_ptr<TH2D> T2vfatinfoErrorMapGeoQ1;
  std::auto_ptr<TH2D> T2vfatinfoErrorMapGeoQ2;
  std::auto_ptr<TH2D> T2vfatinfoErrorMapGeoQ3;
  std::auto_ptr<TH2D> T2vfatinfoErrorMapGeoQ4;

  std::auto_ptr<TH2D> VfatsCorruptionFrequency2D;
  std::auto_ptr<TH1F> VfatsCorruptionFrequency[40];
  std::auto_ptr<TH1F> CumulativePhiResol;
  std::auto_ptr<TH1F> CumulativeRResol;
  std::auto_ptr<TH1F> CumulativeExpXUncert;
  std::auto_ptr<TH1F> CumulativeExpYUncert;



  std::auto_ptr<TH1D> NumCl1HitInOvRegionArm0_Up;
  std::auto_ptr<TH1D> NumCl1HitInOvRegionArm0_Down;


  
  std::auto_ptr<TProfile> NumCluPadOverlap_perPlane;
  std::auto_ptr<TProfile> Tr2EffDet0;
  std::auto_ptr<TProfile> Tr2EffDet1;
  std::auto_ptr<TProfile> Tr2EffDet2;
  std::auto_ptr<TProfile> Tr2EffDet3;
  std::auto_ptr<TProfile> Tr2EffDet4;
  std::auto_ptr<TProfile> Tr2EffDetAll;

  std::auto_ptr<TProfile> Tr2EffDet0strip;
  std::auto_ptr<TProfile> Tr2EffDet1strip;
  std::auto_ptr<TProfile> Tr2EffDet2strip;
  std::auto_ptr<TProfile> Tr2EffDet3strip;
  std::auto_ptr<TProfile> Tr2EffDet4strip;
  std::auto_ptr<TProfile> Tr2EffDetAllstrip;
  std::auto_ptr<TProfile> Tr2EffDet0pad;
  std::auto_ptr<TProfile> Tr2EffDet1pad;
  std::auto_ptr<TProfile> Tr2EffDet2pad;
  std::auto_ptr<TProfile> Tr2EffDet3pad;
  std::auto_ptr<TProfile> Tr2EffDet4pad;
  std::auto_ptr<TProfile> Tr2EffDetAllpad;
  std::auto_ptr<TH1F> AllHitR;
  
  std::auto_ptr<TH1F> PSimHitXD1;
  std::auto_ptr<TH1F> PSimHitYD1;
  
  std::auto_ptr<TH1F> Plside0PSimHitYD0; 
  std::auto_ptr<TH1F> Plside1PSimHitYD0;
  std::auto_ptr<TH1F> Plside0PSimHitXD0; 
  std::auto_ptr<TH1F> Plside1PSimHitXD0;

  std::auto_ptr<TProfile>  NumAllHitvsIdGoodevt;
  std::auto_ptr<TProfile>  NumAllHitvsIdAllevt;

  std::auto_ptr<TCanvas> Chi2PhiProbLogy;
  std::auto_ptr<TCanvas> Chi2RProbLogy;

  std::auto_ptr<TCanvas> CosmTracks;
 
  std::auto_ptr<TH1F> RandomGauss;
  std::auto_ptr<TH1F> HitMatrixR;
  std::auto_ptr<TH1F> HitMatrixPhi;
  std::auto_ptr<TH1F>  ActiveSymbDet;
  std::auto_ptr<TH1F>  PolarAngles;
  std::auto_ptr<TH1F>  Tracksplanestat3;
  std::auto_ptr<TH1F>  Tracksplanestat4;  
  std::auto_ptr<TH1F>  TracksplanestatN;
  std::auto_ptr<TH1F>  TracksplanestatNm1;
  unsigned int trkcounter;
  std::auto_ptr<TH1F> Chi2red;
  std::auto_ptr<TH1F> Chi2phired;
 

  std::vector<unsigned int> ReferenceQuarters;

  double MaxDphi;
  double maxdphihit;
  double maxdrhit; 
  double Effmaxdphihit; 
  double Effmaxdrhit; 



  std::auto_ptr<TProfile> DetIdvsXShOnlyXfit;
  std::auto_ptr<TProfile> DetIdvsYShOnlyYfit;
  std::auto_ptr<TH1F> AllDetId;
  std::auto_ptr<TH1F> AllDetIdPlane;
  std::auto_ptr<TH1F> AllDetIdPlaneS;
  std::auto_ptr<TH1F> TrkChi2X;
  std::auto_ptr<TH1F> TrkChi2Y;
  std::vector<double> alldetZ;
  unsigned int HitNumb4Align;
 

  std::auto_ptr<TH1F> C12AllDetId;
  std::auto_ptr<TH1F> C12AllDetIdPlane;
  std::auto_ptr<TH1F> C12AllDetIdPlaneS;



  //---------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //-------------------------------------------------------------
  std::auto_ptr<TH2F> ReferenceTrackXY;
  std::auto_ptr<TH1D>    HX_ClusterPhi; 


  std::auto_ptr<TProfile> H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[10];
  std::auto_ptr<TProfile> H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_56Division[10];
  std::auto_ptr<TProfile> H0_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult_ZcutOnGeant[10];
  std::auto_ptr<TProfile> H1_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[10];
  std::auto_ptr<TProfile> H2_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[10];
  std::auto_ptr<TProfile> H3_TrkEtaEfficiencyCumulative_UnfoldCutsBinMult[10];
  




  std::auto_ptr<TProfile>    HXTrkEffi;
  std::auto_ptr<TProfile>    HXTrkEffi_AxisRawTrk;

  std::auto_ptr<TProfile>    H0TrkEffi_AxisRawTrk;
  std::auto_ptr<TProfile>    H1TrkEffi_AxisRawTrk;  
  std::auto_ptr<TProfile>    H2TrkEffi_AxisRawTrk;
  std::auto_ptr<TProfile>    H3TrkEffi_AxisRawTrk;
  

  std::auto_ptr<TProfile>    H0TrkEffi_AxisAvgMult;
  std::auto_ptr<TProfile>    H1TrkEffi_AxisAvgMult;
  std::auto_ptr<TProfile>    H2TrkEffi_AxisAvgMult;
  std::auto_ptr<TProfile>    H3TrkEffi_AxisAvgMult;


  std::auto_ptr<TProfile>     H0TrkEffi_AxisRawTrk_EtaCut;
  std::auto_ptr<TProfile>     H1TrkEffi_AxisRawTrk_EtaCut;
  std::auto_ptr<TProfile>     H2TrkEffi_AxisRawTrk_EtaCut;
  std::auto_ptr<TProfile>     H3TrkEffi_AxisRawTrk_EtaCut;


  std::auto_ptr<TProfile>    H0TrkEffi_AxisAvgMult_EtaCut;
  std::auto_ptr<TProfile>    H1TrkEffi_AxisAvgMult_EtaCut;
  std::auto_ptr<TProfile>    H2TrkEffi_AxisAvgMult_EtaCut;
  std::auto_ptr<TProfile>    H3TrkEffi_AxisAvgMult_EtaCut;


  std::auto_ptr<TProfile>    H0TrkEffi_AxisAvgMult_ZImpactCut;
  std::auto_ptr<TProfile>    H1TrkEffi_AxisAvgMult_ZImpactCut;
  std::auto_ptr<TProfile>    H2TrkEffi_AxisAvgMult_ZImpactCut;
  std::auto_ptr<TProfile>    H3TrkEffi_AxisAvgMult_ZImpactCut;

 

  std::auto_ptr<TProfile>    HXTrkEffivsR; 
 
  std::auto_ptr<TProfile>    HXTrkEffivsRplus; 
  std::auto_ptr<TProfile>    HXTrkEffi_Yplus; 
  std::auto_ptr<TProfile>    HXTrkEffi_Yminus; 
  std::auto_ptr<TProfile>    HXTrkEffivsRminus;


  std::auto_ptr<TProfile>   H0TrkEffivsR;
  std::auto_ptr<TProfile>   H1TrkEffivsR;
  std::auto_ptr<TProfile>   H2TrkEffivsR;
  std::auto_ptr<TProfile>   H3TrkEffivsR;
  std::auto_ptr<TProfile>   H0TrkEffi;
  std::auto_ptr<TProfile>   H1TrkEffi;
  std::auto_ptr<TProfile>   H2TrkEffi;
  std::auto_ptr<TProfile>   H3TrkEffi;

  std::auto_ptr<TProfile>   H0TrkEffi_Yplus;
  std::auto_ptr<TProfile>   H1TrkEffi_Yplus;
  std::auto_ptr<TProfile>   H2TrkEffi_Yplus;
  std::auto_ptr<TProfile>   H3TrkEffi_Yplus;
  
  std::auto_ptr<TProfile>   H0TrkEffi_Yminus;
  std::auto_ptr<TProfile>   H1TrkEffi_Yminus;
  std::auto_ptr<TProfile>   H2TrkEffi_Yminus;
  std::auto_ptr<TProfile>   H3TrkEffi_Yminus;

  std::auto_ptr<TProfile>   NumRefTrk_vs_RefAvgHitMultH0;
  std::auto_ptr<TProfile>   NumRefTrk_vs_RefAvgHitMultH1;
  std::auto_ptr<TProfile>   NumRefTrk_vs_RefAvgHitMultH2;
  std::auto_ptr<TProfile>   NumRefTrk_vs_RefAvgHitMultH3;
  std::auto_ptr<TH1D>   UnassHitFractH0;
  std::auto_ptr<TH1D>   UnassHitFractH1;
  std::auto_ptr<TH1D>   UnassHitFractH2;
  std::auto_ptr<TH1D>   UnassHitFractH3;
  std::auto_ptr<TH1D>  AvgRefMultH0;
  std::auto_ptr<TH1D>   AvgRefMultH1;
  std::auto_ptr<TH1D>   AvgRefMultH2;
  std::auto_ptr<TH1D>   AvgRefMultH3;
  std::auto_ptr<TH1D>   H0RefHitActivePlane;
  std::auto_ptr<TH1F>  ZIMPdistrRefH0; 
  std::auto_ptr<TH1F>  ZIMPdistrRefH1;
  std::auto_ptr<TH1F>  ZIMPdistrRefH2; 
  std::auto_ptr<TH1F>  ZIMPdistrRefH3;
 
  std::auto_ptr<TProfile>   DPhiTrkHit_HRefHTestCloseVsZTest;
  std::auto_ptr<TProfile>   DRTrkHit_HRefHTestCloseVsZTest;
  std::auto_ptr<TProfile>   DXTrkHit_HRefHTestCloseVsZTest;
  std::auto_ptr<TProfile>   DYTrkHit_HRefHTestCloseVsZTest;

  std::auto_ptr<TH1D>   QuarterQuarterDXEntry;
  std::auto_ptr<TH1D>   QuarterQuarterDXExit;
  std::auto_ptr<TH1D>   QuarterQuarterDYEntry;
  std::auto_ptr<TH1D>   QuarterQuarterDYExit;


  std::auto_ptr<TH1D>   DPhiTrkHit_HRefHTestClose;
  std::auto_ptr<TH1D>   DRTrkHit_HRefHTestClose;
  std::auto_ptr<TH1D>   DXTrkHit_HRefHTestClose;
  std::auto_ptr<TH1D>   DYTrkHit_HRefHTestClose;

  std::auto_ptr<TH1D> TrkRinRefQuarter;
  std::auto_ptr<TH1D> TrkPhiinRefQuarter;
  std::auto_ptr<TH1D> TrkXinRefQuarter;
  std::auto_ptr<TH1D> TrkYinRefQuarter;
  std::auto_ptr<TH1F> H0TrackPhiALL;
  std::auto_ptr<TH2F> ReferenceZImpVsEta2;
  std::auto_ptr<TH2F> ScatterPlotPhi_ROverlaps;  std::auto_ptr<TH2F> ScatterPlotPhi_R_fullPhiH0;
  std::auto_ptr<TH2F> ScatterPlotPhi_ROverlapsH0;std::auto_ptr<TH2F> ScatterPlotPhi_ROverlapsH1;
  std::auto_ptr<TH2F> ScatterPlotPhi_ROverlapsH2;std::auto_ptr<TH2F> ScatterPlotPhi_ROverlapsH3;
  
  std::auto_ptr<TH1F> HitInOVERLTrk_ALL_Phi_Residual;
  std::auto_ptr<TH1F> HitInOVERLTrk_ALL_R_Residual;
  
  std::auto_ptr<TH1F> DPhiTrk_HRefHTest;
  std::auto_ptr<TH2F> DPhiDRTrk_HRefHTest;
  
  std::auto_ptr<TH1D>   RefTrackXall;
  
  std::auto_ptr<TH1F> DEtaMinRefTrkVsTestQuarterTrkArm0;
  std::auto_ptr<TH1F> DPhiMinRefTrkVsTestQuarterTrkArm0;


  std::auto_ptr<TH2F> HitH0vsHitH1_OvrlTrk;
  std::auto_ptr<TH2F> HitH2vsHitH3_OvrlTrk;

  std::auto_ptr<TH1F> HowManyPadCluCloseToRef_PluNear;
  std::auto_ptr<TH1F> HowManyPadCluCloseToRef_PluFar;
  
  std::auto_ptr<TH1F> HowManyPadCluCloseToRef_ArmXFar; 
  std::auto_ptr<TH1F> HowManyPadCluCloseToRef_ArmXNear;
  
  std::auto_ptr<TH1F> ReasonOfMissingTrkArm0PluNear;
  std::auto_ptr<TH1F> ReasonOfMissingTrkArm0PluFar;
   // 
 std::auto_ptr<TH1F>  ReasonOfMissingTrkArmXFar;  
 std::auto_ptr<TH1F> ReasonOfMissingTrkArmXNear;
  
  
  std::auto_ptr<TH1D> TrkPhiinTestedQuarterCloseToRef; 
  std::auto_ptr<TH1D> TrkRinTestedQuarterCloseToRef ;
  std::auto_ptr<TH1D> TrkXinTestedQuarterCloseToRef ;
  std::auto_ptr<TH1D> TrkYinTestedQuarterCloseToRef ;
  std::auto_ptr<TH1D> DPhiTrk_HRefHTest2CloseToRef;

 
  unsigned int RefTrkHitMult;
  double UnassHitThr;

  double Chi2ProbRefTrk ;
  bool UseRestrictedOverlap ;
  bool OldReco;

  unsigned int  MaxPadCluInOverlapUporDown ;
  unsigned int NumHitGood ;
  unsigned int numevent;
 
  double TrkEtamin;
  double TrkEtaMAX;
  double PhiChiProbCut;
  double RChiProbCut;
  double Rwind;
  double Phiwind;
  int numhittrkeff;
  double chiRredCut,chiPhiredCut;
 
  bool ExcludeNoisyplane;
  unsigned int CommonNoiseClSize;
 
  unsigned int MaxPad, MaxStrip;
  unsigned int EffMaxPad, EffMaxStrip;

  double ShadowTan;	  
  double AllowedDRTrackDistance;
  unsigned int  MaxEvents;
  unsigned int RefTrkMultiplicity;
  unsigned int MaxTrkInProcess;
  std::string outputFileName, CluLabel,RefCluLabel, HitLabel, RoadLabel,TrackLabel;
  std::string Refstring;
 
  std::string TrackLabelLastPlanes;
  std::string TrackLabelFirstPlanes;
  std::string PadCluLabelFirstPlanes;
  std::string PadCluLabelLastPlanes;

  std::string RoadLabelH0, TrackLabelH0, RoadLabelH1, TrackLabelH1;
  std::string RoadLabelH2, TrackLabelH2, RoadLabelH3, TrackLabelH3;
  std::string HitLabelTestedQ;


  std::string xmlfilenameFull;
  std::string xmlfilenameUsed_NotDead;
  std::string AnType; //def=INCLUSIVE, DOUBLE SINGLE
  bool produceVfatEffiFile;
  bool simufile;
  bool verbosity;
  
  bool skipSelectedEvents;
  std::string skipEventFileName;

  bool requiregoodChi;
  bool LookToRawEvent;
  bool UseUncorrupetdEventMap; //made corruption analysis
  bool OnlyShadowAnalysis;  


  bool trackCollectionH0FirstsValid;
  bool trackCollectionH1FirstsValid;
  bool trackCollectionH2FirstsValid;
  bool trackCollectionH3FirstsValid;
  bool trackCollectionH0LastsValid;
  bool trackCollectionH1LastsValid;
  bool trackCollectionH2LastsValid;
  bool trackCollectionH3LastsValid;



  unsigned int MaxplaneCorruptInquarter;
};

#endif 

// -*- C++ -*-
//
// Original Author:  Mirko Berretti
//         Created:  Wed Dec 12 18:19:44 CET 2007
// $Id: T2AnalyzerRaw.h,v 1.1.2.1 2009/11/08 14:50:48 berretti Exp $
//

#ifndef _TotemT1T2ValidationT2RecoValidationT2RecoAnalyzerRaw_H_
#define _TotemT1T2ValidationT2RecoValidationT2RecoAnalyzerRaw_H_

// system include files
#include <memory>
#include <math.h>
#include <TMath.h>
// user include files
#include "TotemT1T2Validation/T2GeometryValidation/interface/DAQInformationSourceXML_a.h"
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include "TCanvas.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TPad.h"
#include "TView.h"
#include "TAxis3D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPolyLine3D.h"
//#include <THnSparse.h>
#include "THnSparse.h"
//#include "TotemBackground/BeamGas/interface/BeamGasBckgStudy.h"
//#include "TotemBackground/BeamGas/interface/THnSparse.h"
#include "TotemAnalysis/T2ValidRawData/interface/TrackInfo.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2D.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2DCollection.h"
#include "DataFormats/T2Cluster/interface/T2Cluster.h"
#include "DataFormats/T2Hit/interface/T2Hit.h"
#include "DataFormats/T2Hit/interface/T2HitCollection.h"
#include "DataFormats/T2Cluster/interface/T2PadClusterCollection.h"
#include "DataFormats/T2Cluster/interface/T2StripClusterCollection.h"
#include "DataFormats/T2DetId/interface/T2DetId.h"
#include "DataFormats/T2Road/interface/T2Road.h"
#include "DataFormats/T2Road/interface/T2RoadCollection.h"
#include "DataFormats/T1T2Track/interface/T1T2Track.h"
#include "DataFormats/T1T2Track/interface/T1T2TrackCollection.h"
#include "TotemAnalysis/T2Cuts/interface/T2SelectionCutUtils.h"

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
#include "DataFormats/T2DigiVfat/interface/T2DigiVfat.h"
#include "DataFormats/T2DigiVfat/interface/T2DigiVfatCollection.h"



// class decleration
//

class T2AnalyzerRaw : public edm::EDAnalyzer {
public:
  explicit T2AnalyzerRaw(const edm::ParameterSet&);
  ~T2AnalyzerRaw();
  //unsigned int T2AnalyzerRaw::RawtoSymb(uint32_t thedet);
  //std::vector<double> MyLinearfitX(std::vector<T2Hit> hitvec,bool UseJointProb);
  //std::vector<double> MyLinearfitY(std::vector<T2Hit> hitvec,bool UseJointProb);
  // std::vector<double> MyLinearfit(std::vector<T2Hit> hitvec,bool UseJointProb);
  DAQInformationSourceXML_a* onemap;
  boost::shared_ptr<DAQInformationT2_a> Map_NOTDEAD;

  DAQInformationSourceXML_a* MapProducer_FullVFATs;
  boost::shared_ptr<DAQInformationT2_a> Map_FullVfats;
  

  DAQInformationSourceXML_a* MapProducer_NotDeadUsedVfat;
  boost::shared_ptr<DAQInformationT2_a> Map_NotDeadUsedVfat;



  bool IsVfatMapped(int symbvfat);
  virtual void beginJob();
  
  std::vector<int> AllCorruptedEvents;


  std::vector<bool> vfatsCorrupted;
  std::vector<int> vfatsSymbId;




private:

  edm::InputTag t2PadDigiCollectionLabel;
  edm::InputTag t2StripDigiCollectionLabel;
  edm::InputTag t2DigiVfatCollectionLabel;
  edm::InputTag rawEventLabel;

  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  std::vector<TPolyLine3D*> All3dtracks;
  std::vector<TPolyLine3D*> tracks0123;
  std::vector<TPolyLine3D*> tracks0124;
  std::vector<TPolyLine3D*> tracks0234;
  std::vector<TPolyLine3D*> tracks1234;
  std::vector<TPolyLine3D*> tracks0134;

  std::vector<bool> ExludeThisNoisyplaneForEffMeas(T2StripClusterCollection::const_iterator itstrip, T2PadClusterCollection::const_iterator itpad);
  
  bool PlaneInRightHalf(uint32_t cmsswid);

  bool MettiTracciaBuona(std::vector<TrackInfo>* matricetracce,T1T2TrackCollection::const_iterator  trk, unsigned int numevento);

  bool HitIndetector(std::vector<double> vpar, unsigned int symbdetid/*,std::map<unsigned int,double> mapsymbtoZ*/);

  bool HitInQuarterAnalyzed(T2Hit &aHit);

  bool CloseTrksInEvt(T1T2TrackCollection trackColl);

  T2Hit GiveExtrapolatedHit(int plane0_10, std::vector<T2Hit> referenceHits);
  
  bool HitIsInTrackColl(T2Hit trackinghit,T1T2TrackCollection trackColl);

  std::vector<double> ResiduiForStrip(std::vector<T2Hit>refhitv, T2Hit hit);
  //std::vector<double> ResiduiForPad(std::vector<T2Hit>refhitv, T2Hit hit);
  std::vector<double> ResiduiForPad(std::vector<T2Hit>refhitv, T2Hit hit,double &expPadXErr,double &expPadYErr);

  bool PadIsNeighbour(std::vector<T2Hit>refhitv,T2Hit clusterchecked);


  std::vector<double> ResiduiForC1HIT(std::vector<T2Hit>refhitv, T2Hit hit);

  bool CheckEventCorruption();
  

  //  void LoadDeadSectorFile();
  
  bool HitInDeadSector(unsigned int Rsector,unsigned int plane0_10);

  void ConvertCluRPhiInEffiGeoCell(double Rcl, double Phicl, uint32_t thedet, int &rCell, int &phiCell);



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
std::auto_ptr<TProfile>  PadCluMultiplVsPlane;

  std::auto_ptr<TH2D> Class1HitPadStripCLSCorrel;
  
   std::auto_ptr<TH1F> EventCandidate;
  std::auto_ptr<TH1F> clusterstripentries;
  std::auto_ptr<TH1F> clusterpadentries;
  std::auto_ptr<TH1F> diffphiCluGun;
  std::auto_ptr<TH1F> RLocHRecoH;
  std::auto_ptr<TH1F> Trketa;
  std::auto_ptr<TH1F> Trkphi;
  std::auto_ptr<TH1F> TrkphiRZPlus;
  std::auto_ptr<TH1F> TrkphiRZMinus;

  std::auto_ptr<TH1F> RelEffi_ForTracking;

  std::auto_ptr<TH1F> TrkphiALL;
  std::auto_ptr<TH1F> TrketaALL;
  std::auto_ptr<TH1F> NumhitinTrackALL;
  std::auto_ptr<TH1F> TrkQuarterIdALL;
  

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

  

  std::auto_ptr<TProfile> EFFstripVsIdNoCut;
  std::auto_ptr<TProfile> EFFstripVsId;
  std::auto_ptr<TProfile> EFFpadVsIdNoCut;
  std::auto_ptr<TProfile> EFFpadVsId;
  std::auto_ptr<TProfile> EFFHitVsIdNoCut;
  std::auto_ptr<TProfile> EFFHitVsId;
  std::auto_ptr<TProfile> EFFClusterVsId;

//R assumed :43-145 spacer: 50, 92.5 121.5 pm 1 -> 30 means  each cell=3.4 mm: ~2 points for sect 1 ~12 for central, ~8 for outer
  std::auto_ptr<TProfile> EFFpadVsIdNormalized;
  std::auto_ptr<TProfile> EFFstripVsIdNormalized;
  std::auto_ptr<TProfile> EFFHitVsIdNormalized;


  std::auto_ptr<TProfile> EFFpadVsIdNormalized_sectorX[4];
  std::auto_ptr<TProfile> EFFstripVsIdNormalized_sectorX[4];
  
  
  std::auto_ptr<TH1F> EFF_DRstripNoCutdet0;
  std::auto_ptr<TH1F> EFF_DRstripdet0;
  std::auto_ptr<TH1F> EFF_DRpadNoCutdet0;
  std::auto_ptr<TH1F> EFF_DRpaddet0;
  std::auto_ptr<TH1F> EFF_DRHitNoCutdet0;
  std::auto_ptr<TH1F> EFF_DRHitdet0;
  std::auto_ptr<TH1F> EFF_DRstripNoCutdet2;
  std::auto_ptr<TH1F> EFF_DRstripdet2;
  std::auto_ptr<TH1F> EFF_DRpadNoCutdet2;
  std::auto_ptr<TH1F> EFF_DRpaddet2;
  std::auto_ptr<TH1F> EFF_DRHitNoCutdet2;
  std::auto_ptr<TH1F> EFF_DRHitdet2;
  std::auto_ptr<TH1F> EFF_DRstripNoCutdet4;
  std::auto_ptr<TH1F> EFF_DRstripdet4;
  std::auto_ptr<TH1F> EFF_DRpadNoCutdet4;
  std::auto_ptr<TH1F> EFF_DRpaddet4;
  std::auto_ptr<TH1F> EFF_DRHitNoCutdet4;
  std::auto_ptr<TH1F> EFF_DRHitdet4;
  std::auto_ptr<TH1F> EFF_DRstripNoCutdet6;
  std::auto_ptr<TH1F> EFF_DRstripdet6;
  std::auto_ptr<TH1F> EFF_DRpadNoCutdet6;
  std::auto_ptr<TH1F> EFF_DRpaddet6;
  std::auto_ptr<TH1F> EFF_DRHitNoCutdet6;
  std::auto_ptr<TH1F> EFF_DRHitdet6;
  std::auto_ptr<TH1F> EFF_DRstripNoCutdet8;
  std::auto_ptr<TH1F> EFF_DRstripdet8;
  std::auto_ptr<TH1F> EFF_DRpadNoCutdet8;
  std::auto_ptr<TH1F> EFF_DRpaddet8;
  std::auto_ptr<TH1F> EFF_DRHitNoCutdet8;
  std::auto_ptr<TH1F> EFF_DRHitdet8;



  std::auto_ptr<TH1F> EFF_DPhistripNoCutdet0;
  std::auto_ptr<TH1F> EFF_DPhistripdet0;
  std::auto_ptr<TH1F> EFF_DPhipadNoCutdet0;
  std::auto_ptr<TH1F> EFF_DPhipaddet0;
  std::auto_ptr<TH1F> EFF_DPhiHitNoCutdet0;
  std::auto_ptr<TH1F> EFF_DPhiHitdet0;
  std::auto_ptr<TH1F> EFF_DPhistripNoCutdet2;
  std::auto_ptr<TH1F> EFF_DPhistripdet2;
  std::auto_ptr<TH1F> EFF_DPhipadNoCutdet2;
  std::auto_ptr<TH1F> EFF_DPhipaddet2;
  std::auto_ptr<TH1F> EFF_DPhiHitNoCutdet2;
  std::auto_ptr<TH1F> EFF_DPhiHitdet2;
  std::auto_ptr<TH1F> EFF_DPhistripNoCutdet4;
  std::auto_ptr<TH1F> EFF_DPhistripdet4;
  std::auto_ptr<TH1F> EFF_DPhipadNoCutdet4;
  std::auto_ptr<TH1F> EFF_DPhipaddet4;
  std::auto_ptr<TH1F> EFF_DPhiHitNoCutdet4;
  std::auto_ptr<TH1F> EFF_DPhiHitdet4;
  std::auto_ptr<TH1F> EFF_DPhistripNoCutdet6;
  std::auto_ptr<TH1F> EFF_DPhistripdet6;
  std::auto_ptr<TH1F> EFF_DPhipadNoCutdet6;
  std::auto_ptr<TH1F> EFF_DPhipaddet6;
  std::auto_ptr<TH1F> EFF_DPhiHitNoCutdet6;
  std::auto_ptr<TH1F> EFF_DPhiHitdet6;
  std::auto_ptr<TH1F> EFF_DPhistripNoCutdet8;
  std::auto_ptr<TH1F> EFF_DPhistripdet8;
  std::auto_ptr<TH1F> EFF_DPhipadNoCutdet8;
  std::auto_ptr<TH1F> EFF_DPhipaddet8;
  std::auto_ptr<TH1F> EFF_DPhiHitNoCutdet8;
  std::auto_ptr<TH1F> EFF_DPhiHitdet8;


  std::auto_ptr<TH1F> ProbYHisto;
  std::auto_ptr<TH1F> ProbXHisto;
  std::auto_ptr<TH1F> ProbYHistoAl;
  std::auto_ptr<TH1F> ProbXHistoAl;

  std::auto_ptr<TH1F>  AXError;
  std::auto_ptr<TH1F>   AYError;
  std::auto_ptr<TH1F>   BXError;
  std::auto_ptr<TH1F>    BYError;
  std::auto_ptr<TH1F>     AXHisto;
  std::auto_ptr<TH1F>     AYHisto;

  std::auto_ptr<TH1F> EFF_CLSHitdet0;
  std::auto_ptr<TH1F> EFF_CLSStripdet0;
  std::auto_ptr<TH1F> EFF_CLSPaddet0;
  std::auto_ptr<TH1F> EFF_CLSHitdet2;
  std::auto_ptr<TH1F> EFF_CLSStripdet2;
  std::auto_ptr<TH1F> EFF_CLSPaddet2;
  std::auto_ptr<TH1F> EFF_CLSHitdet4;
  std::auto_ptr<TH1F> EFF_CLSStripdet4;
  std::auto_ptr<TH1F> EFF_CLSPaddet4;
  std::auto_ptr<TH1F> EFF_CLSHitdet6;
  std::auto_ptr<TH1F> EFF_CLSStripdet6;
  std::auto_ptr<TH1F> EFF_CLSPaddet6;
  std::auto_ptr<TH1F> EFF_CLSHitdet8;
  std::auto_ptr<TH1F> EFF_CLSStripdet8;
  std::auto_ptr<TH1F> EFF_CLSPaddet8;
  std::auto_ptr<TH1F> EFF_CLSHitdet0NoCut;
  std::auto_ptr<TH1F> EFF_CLSStripdet0NoCut;
  std::auto_ptr<TH1F> EFF_CLSPaddet0NoCut;
  std::auto_ptr<TH1F> EFF_CLSHitdet2NoCut;
  std::auto_ptr<TH1F> EFF_CLSStripdet2NoCut;
  std::auto_ptr<TH1F> EFF_CLSPaddet2NoCut;
  std::auto_ptr<TH1F> EFF_CLSHitdet4NoCut;
  std::auto_ptr<TH1F> EFF_CLSStripdet4NoCut;
  std::auto_ptr<TH1F> EFF_CLSPaddet4NoCut;
  std::auto_ptr<TH1F> EFF_CLSHitdet6NoCut;
  std::auto_ptr<TH1F> EFF_CLSStripdet6NoCut;
  std::auto_ptr<TH1F> EFF_CLSPaddet6NoCut;
  std::auto_ptr<TH1F> EFF_CLSHitdet8NoCut;
  std::auto_ptr<TH1F> EFF_CLSStripdet8NoCut;
  std::auto_ptr<TH1F> EFF_CLSPaddet8NoCut;
  
  

  
  std::auto_ptr<TH1F>   NumhitinTrackAligned;
  std::auto_ptr<TH1F> EFF_DRstripNoCutdet1;
  std::auto_ptr<TH1F> EFF_DRstripdet1;
  std::auto_ptr<TH1F> EFF_DRpadNoCutdet1;
  std::auto_ptr<TH1F> EFF_DRpaddet1;
  std::auto_ptr<TH1F> EFF_DRHitNoCutdet1;
  std::auto_ptr<TH1F> EFF_DRHitdet1;
  std::auto_ptr<TH1F> EFF_DRstripNoCutdet3;
  std::auto_ptr<TH1F> EFF_DRstripdet3;
  std::auto_ptr<TH1F> EFF_DRpadNoCutdet3;
  std::auto_ptr<TH1F> EFF_DRpaddet3;
  std::auto_ptr<TH1F> EFF_DRHitNoCutdet3;
  std::auto_ptr<TH1F> EFF_DRHitdet3;
  std::auto_ptr<TH1F> EFF_DRstripNoCutdet5;
  std::auto_ptr<TH1F> EFF_DRstripdet5;
  std::auto_ptr<TH1F> EFF_DRpadNoCutdet5;
  std::auto_ptr<TH1F> EFF_DRpaddet5;
  std::auto_ptr<TH1F> EFF_DRHitNoCutdet5;
  std::auto_ptr<TH1F> EFF_DRHitdet5;
  std::auto_ptr<TH1F> EFF_DRstripNoCutdet7;
  std::auto_ptr<TH1F> EFF_DRstripdet7;
  std::auto_ptr<TH1F> EFF_DRpadNoCutdet7;
  std::auto_ptr<TH1F> EFF_DRpaddet7;
  std::auto_ptr<TH1F> EFF_DRHitNoCutdet7;
  std::auto_ptr<TH1F> EFF_DRHitdet7;
  std::auto_ptr<TH1F> EFF_DRstripNoCutdet9;
  std::auto_ptr<TH1F> EFF_DRstripdet9;
  std::auto_ptr<TH1F> EFF_DRpadNoCutdet9;
  std::auto_ptr<TH1F> EFF_DRpaddet9;
  std::auto_ptr<TH1F> EFF_DRHitNoCutdet9;
  std::auto_ptr<TH1F> EFF_DRHitdet9;
  



  std::auto_ptr<TH1F> EFF_DPhistripNoCutdet1;
  std::auto_ptr<TH1F> EFF_DPhistripdet1;
  std::auto_ptr<TH1F> EFF_DPhipadNoCutdet1;
  std::auto_ptr<TH1F> EFF_DPhipaddet1;
  std::auto_ptr<TH1F> EFF_DPhiHitNoCutdet1;
  std::auto_ptr<TH1F> EFF_DPhiHitdet1;
  std::auto_ptr<TH1F> EFF_DPhistripNoCutdet3;
  std::auto_ptr<TH1F> EFF_DPhistripdet3;
  std::auto_ptr<TH1F> EFF_DPhipadNoCutdet3;
  std::auto_ptr<TH1F> EFF_DPhipaddet3;
  std::auto_ptr<TH1F> EFF_DPhiHitNoCutdet3;
  std::auto_ptr<TH1F> EFF_DPhiHitdet3;
  std::auto_ptr<TH1F> EFF_DPhistripNoCutdet5;
  std::auto_ptr<TH1F> EFF_DPhistripdet5;
  std::auto_ptr<TH1F> EFF_DPhipadNoCutdet5;
  std::auto_ptr<TH1F> EFF_DPhipaddet5;
  std::auto_ptr<TH1F> EFF_DPhiHitNoCutdet5;
  std::auto_ptr<TH1F> EFF_DPhiHitdet5;
  std::auto_ptr<TH1F> EFF_DPhistripNoCutdet7;
  std::auto_ptr<TH1F> EFF_DPhistripdet7;
  std::auto_ptr<TH1F> EFF_DPhipadNoCutdet7;
  std::auto_ptr<TH1F> EFF_DPhipaddet7;
  std::auto_ptr<TH1F> EFF_DPhiHitNoCutdet7;
  std::auto_ptr<TH1F> EFF_DPhiHitdet7;
  std::auto_ptr<TH1F> EFF_DPhistripNoCutdet9;
  std::auto_ptr<TH1F> EFF_DPhistripdet9;
  std::auto_ptr<TH1F> EFF_DPhipadNoCutdet9;
  std::auto_ptr<TH1F> EFF_DPhipaddet9;
  std::auto_ptr<TH1F> EFF_DPhiHitNoCutdet9;
  std::auto_ptr<TH1F> EFF_DPhiHitdet9;


  std::auto_ptr<TH1F> EFF_C1CLSStripdet0;
  std::auto_ptr<TH1F> EFF_C1CLSPaddet0;
  std::auto_ptr<TH1F> EFF_C1CLSStripdet1;
  std::auto_ptr<TH1F> EFF_C1CLSPaddet1;
  std::auto_ptr<TH1F> EFF_C1CLSStripdet2;
  std::auto_ptr<TH1F> EFF_C1CLSPaddet2;
  std::auto_ptr<TH1F> EFF_C1CLSStripdet3;
  std::auto_ptr<TH1F> EFF_C1CLSPaddet3;
  std::auto_ptr<TH1F> EFF_C1CLSStripdet4;
  std::auto_ptr<TH1F> EFF_C1CLSPaddet4;
  std::auto_ptr<TH1F> EFF_C1CLSStripdet5;
  std::auto_ptr<TH1F> EFF_C1CLSPaddet5;

  std::auto_ptr<TH1F> EFF_C1CLSStripdet6;
  std::auto_ptr<TH1F> EFF_C1CLSPaddet6;
  std::auto_ptr<TH1F> EFF_C1CLSStripdet7;
  std::auto_ptr<TH1F> EFF_C1CLSPaddet7;
  std::auto_ptr<TH1F> EFF_C1CLSStripdet8;
  std::auto_ptr<TH1F> EFF_C1CLSPaddet8;
  std::auto_ptr<TH1F> EFF_C1CLSStripdet9;
  std::auto_ptr<TH1F> EFF_C1CLSPaddet9;






  
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
  


  std::auto_ptr<TH1F> EFF_CLSHitdet1;
  std::auto_ptr<TH1F> EFF_CLSStripdet1;
  std::auto_ptr<TH1F> EFF_CLSPaddet1;
  std::auto_ptr<TH1F> EFF_CLSHitdet3;
  std::auto_ptr<TH1F> EFF_CLSStripdet3;
  std::auto_ptr<TH1F> EFF_CLSPaddet3;
  std::auto_ptr<TH1F> EFF_CLSHitdet5;
  std::auto_ptr<TH1F> EFF_CLSStripdet5;
  std::auto_ptr<TH1F> EFF_CLSPaddet5;
  std::auto_ptr<TH1F> EFF_CLSHitdet7;
  std::auto_ptr<TH1F> EFF_CLSStripdet7;
  std::auto_ptr<TH1F> EFF_CLSPaddet7;
  std::auto_ptr<TH1F> EFF_CLSHitdet9;
  std::auto_ptr<TH1F> EFF_CLSStripdet9;
  std::auto_ptr<TH1F> EFF_CLSPaddet9;
  std::auto_ptr<TH1F> EFF_CLSHitdet1NoCut;
  std::auto_ptr<TH1F> EFF_CLSStripdet1NoCut;
  std::auto_ptr<TH1F> EFF_CLSPaddet1NoCut;
  std::auto_ptr<TH1F> EFF_CLSHitdet3NoCut;
  std::auto_ptr<TH1F> EFF_CLSStripdet3NoCut;
  std::auto_ptr<TH1F> EFF_CLSPaddet3NoCut;
  std::auto_ptr<TH1F> EFF_CLSHitdet5NoCut;
  std::auto_ptr<TH1F> EFF_CLSStripdet5NoCut;
  std::auto_ptr<TH1F> EFF_CLSPaddet5NoCut;
  std::auto_ptr<TH1F> EFF_CLSHitdet7NoCut;
  std::auto_ptr<TH1F> EFF_CLSStripdet7NoCut;
  std::auto_ptr<TH1F> EFF_CLSPaddet7NoCut;
  std::auto_ptr<TH1F> EFF_CLSHitdet9NoCut;
  std::auto_ptr<TH1F> EFF_CLSStripdet9NoCut;
  std::auto_ptr<TH1F> EFF_CLSPaddet9NoCut;


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


  std::auto_ptr<TProfile> NumPadCluVsPlane; 
  std::auto_ptr<TProfile> NumStripCluVsPlane; 
  std::auto_ptr<TProfile> NumCl1HitVsPlane; 
  std::auto_ptr<TProfile> PadCluSizeVsPlaneAll; 
  std::auto_ptr<TProfile> PadCluSizeVsPlaneAll2;


  std::auto_ptr<TProfile> PadCluSizeVsPlaneAll3H0;
  std::auto_ptr<TProfile> PadCluSizeVsPlaneAll3H1;
  std::auto_ptr<TProfile> PadCluSizeVsPlaneAll3H2;
  std::auto_ptr<TProfile> PadCluSizeVsPlaneAll3H3;

  std::auto_ptr<TProfile> StripCluSizeVsPlaneAll3H0;
  std::auto_ptr<TProfile> StripCluSizeVsPlaneAll3H1;
  std::auto_ptr<TProfile> StripCluSizeVsPlaneAll3H2;
  std::auto_ptr<TProfile> StripCluSizeVsPlaneAll3H3;

  std::auto_ptr<TProfile> NumPadCluVsPlaneAll3H0; 
  std::auto_ptr<TProfile> NumStripCluVsPlaneAll3H0;
  std::auto_ptr<TProfile> NumPadCluVsPlaneAll3H1; 
  std::auto_ptr<TProfile> NumStripCluVsPlaneAll3H1;
  std::auto_ptr<TProfile> NumPadCluVsPlaneAll3H2; 
  std::auto_ptr<TProfile> NumStripCluVsPlaneAll3H2;
  std::auto_ptr<TProfile> NumPadCluVsPlaneAll3H3; 
  std::auto_ptr<TProfile> NumStripCluVsPlaneAll3H3;

  std::auto_ptr<TProfile>  NumStripCluVsPlaneAll3H0_Cutted;
  std::auto_ptr<TProfile>  NumStripCluVsPlaneAll3H1_Cutted;
  std::auto_ptr<TProfile>  NumStripCluVsPlaneAll3H2_Cutted;
  std::auto_ptr<TProfile>  NumStripCluVsPlaneAll3H3_Cutted;
 
  std::auto_ptr<TProfile>  NumPadCluVsPlaneAll3H0_Cutted;
  std::auto_ptr<TProfile>  NumPadCluVsPlaneAll3H1_Cutted;
  std::auto_ptr<TProfile>  NumPadCluVsPlaneAll3H2_Cutted;
  std::auto_ptr<TProfile>  NumPadCluVsPlaneAll3H3_Cutted;

  std::auto_ptr<TProfile> NumPadCluVsPlaneAll3H0_CuttedMinus;
  std::auto_ptr<TProfile> NumPadCluVsPlaneAll3H1_CuttedMinus;
  std::auto_ptr<TProfile> NumPadCluVsPlaneAll3H2_CuttedMinus;
  std::auto_ptr<TProfile> NumPadCluVsPlaneAll3H3_CuttedMinus;
  std::auto_ptr<TProfile> NumStripCluVsPlaneAll3H0_CuttedMinus;
  std::auto_ptr<TProfile> NumStripCluVsPlaneAll3H1_CuttedMinus;
  std::auto_ptr<TProfile> NumStripCluVsPlaneAll3H2_CuttedMinus;
  std::auto_ptr<TProfile> NumStripCluVsPlaneAll3H3_CuttedMinus;

  std::auto_ptr<TProfile> NumPadCluVsPlaneAll3H0_CuttedPlus;
  std::auto_ptr<TProfile> NumPadCluVsPlaneAll3H1_CuttedPlus;
  std::auto_ptr<TProfile> NumPadCluVsPlaneAll3H2_CuttedPlus;
  std::auto_ptr<TProfile> NumPadCluVsPlaneAll3H3_CuttedPlus;
  std::auto_ptr<TProfile> NumStripCluVsPlaneAll3H0_CuttedPlus;
  std::auto_ptr<TProfile> NumStripCluVsPlaneAll3H1_CuttedPlus;
  std::auto_ptr<TProfile> NumStripCluVsPlaneAll3H2_CuttedPlus;
  std::auto_ptr<TProfile> NumStripCluVsPlaneAll3H3_CuttedPlus;






  std::auto_ptr<TProfile>  NumPadCluVsPlaneAll3H0_Cutted_LOWMultipl;
  std::auto_ptr<TProfile>  NumPadCluVsPlaneAll3H1_Cutted_LOWMultipl;
  std::auto_ptr<TProfile>  NumPadCluVsPlaneAll3H2_Cutted_LOWMultipl;
  std::auto_ptr<TProfile>  NumPadCluVsPlaneAll3H3_Cutted_LOWMultipl;
 




  std::auto_ptr<TH1F> CumulativeNumPadCluAll3H0;
  std::auto_ptr<TH1F> CumulativePadCluSizeAll3H0;
  std::auto_ptr<TH1F> CumulativeNumStripCluAll3H0;
  std::auto_ptr<TH1F> CumulativeStripCluSizeAll3H0;

  
  std::auto_ptr<TH1F> CumulativePadCluSize_UsedInEffi_HX;
  std::auto_ptr<TH1F> CumulativeStripCluSize_UsedInEffi_HX;


  std::auto_ptr<TProfile> StipCluSizeVsPlaneAll;


  std::auto_ptr<TH2D> PlaneNumPadClu11_vs_16;

  std::auto_ptr<TH2D> PlaneHitRadiographyXY[40];
  std::auto_ptr<TH2D> PadClusterXYRadiographyYPlus[40];
  std::auto_ptr<TH2D> PadClusterXYRadiographyYMinus[40];
  std::auto_ptr<TH1D> PadClusterR[40];
  
  std::auto_ptr<TH2D> StripCluRadiographyXY[40];
  std::auto_ptr<TH2D> PadClusterXYForEffiEvents[10];
  std::auto_ptr<TH2D> HalfTeleTrkRadiographyXY[4];
  
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

  std::auto_ptr<TH1D> vFatMultiplicty[40][17];
  std::auto_ptr<TH1D> vFatCumulative[40][17];
  std::auto_ptr<TProfile> vFatMultiplictyVsTrk[40][17];
  std::auto_ptr<TH1D> TrkMultiplictyWhenEveryThingIsOn[40][17];
  std::auto_ptr<TH2D> Vfat15_Pl19and17MultCorrel;  
  std::auto_ptr<TH2D> Vfat15_Pl19and11_Pl19MultCorrel;
  std::auto_ptr<TH2D> Vfat15_Pl19and17MultCorrel_LowMulti;
  std::auto_ptr<THnSparseD> EvtVfat_Strip_WhenCompletelyOn; 
  std::auto_ptr<THnSparseS> EvtVfat_Strip_WhenCompletelyOnS; 
  
  std::vector<std::pair<int,int> > NoisyStripVfatH1_PlIId;
  std::auto_ptr<TH1D> EvtVfat_Strip_IdWhenCompletelyOn_ID;

  std::auto_ptr<TH1D> NumH1VfatCompletelyOnPerEvt;
  std::auto_ptr<TH2D> Vfat15_Pl10andH1Vfat15_PlDispariMultCorrel[5];  //11-15 13-15 15-15 17-15 19-15 
  std::auto_ptr<TH2D> Vfat16_Pl10andH1Vfat16_PlDispariMultCorrel[5];   //11-16 13-16 15-16 17-16 19-16
  
  //Noise correlation Histogr
  std::auto_ptr<TH2F> Class1NoiseDRDPHISep[10];  
  std::auto_ptr<TProfile> TrackingHitVsClass1NoiseDRSep[10];  
  std::auto_ptr<TProfile> TrackingHitRVsStripCluNoiseDRSep[10];  
  std::auto_ptr<TProfile> TrackingHitRVsPadCluNoiseDRSep[10];

  // std::auto_ptr<TH1F> PadClusterR_AllvsPlane_LowMulti[40]; 
  std::auto_ptr<TH1F> TrkHitR_vsplane[10];
  std::auto_ptr<TH1F> PadClusterR_AllvsPlane[40]; 
  std::auto_ptr<TH1F> PadClusterSize_AllvsPlane[40];
  std::auto_ptr<TH1F> PAdNr[10];
  std::auto_ptr<TH1F> StripNr[10];
  std::auto_ptr<TProfile> VFATEFF[10];
  std::auto_ptr<TProfile> VFATEFFNormalized[10];
  std::auto_ptr<TH1F> VfatStatistics[10];


  std::vector<std::vector<std::vector<double> > > Geometry_PadEfficiency_Num;
  std::vector<std::vector<std::vector<double> > > Geometry_PadEfficiency_Den;
  std::vector<std::vector<std::vector<double> > > Geometry_PadEfficiency;

  std::vector<std::vector<std::vector<double> > > Geometry_StripEfficiency_Num;
  std::vector<std::vector<std::vector<double> > > Geometry_StripEfficiency_Den;
  std::vector<std::vector<std::vector<double> > > Geometry_StripEfficiency;
  int numRsectEffi;
  int numPhisectEffi;
  
  std::auto_ptr<TH2D> HGeometry_StripEfficiency_Num[40]; 
  std::auto_ptr<TH2D> HGeometry_StripEfficiency_Den[40]; 
  std::auto_ptr<TH2D> HGeometry_PadEfficiency_Num[40]; 
  std::auto_ptr<TH2D> HGeometry_PadEfficiency_Den[40]; 

  //std::auto_ptr<TH2D> Geometry_PadEfficiency[8][24]; //[num_r_sector][num_pad_sector]
  //std::auto_ptr<TH2D> Geometry_StripEfficiency[8][24];

  //192/24= 8 gradi, 24 phi sector
  //8 Rsector 
  std::auto_ptr<TH2D> Geometry_HitEfficiency_caseON[30][15]; //[num_r_sector][num_pad_sector] 
//R assumed :43-145 spacer: 50, 92.5 121.5 pm 1 -> 30 means  each cell=3.4 mm: ~2 points for sect 1 ~12 for central, ~8 for outer
 
  std::auto_ptr<TH2D> T2vfatinfoErrorMap;         
  std::auto_ptr<TH2D> T2vfatinfoErrorMapNOTDEAD;
  std::auto_ptr<TH2D> T2vfatinfoErrorMap_DEADWithCorrectFrame;

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
  
 


  std::auto_ptr<TH1D> EventHemisphere;
  std::auto_ptr<TH1F> AllHitZ;
  std::auto_ptr<TH1F> AllHitR;
   
  std::auto_ptr<TH1F> AllHitRPadOnly; 
  std::auto_ptr<TH1F> AllHitZPadOnly; 
  std::auto_ptr<TH1F> AllHitRStripOnly; 
  std::auto_ptr<TH1F> AllHitZStripOnly; 
  std::auto_ptr<TH1F> AllHitZCL1; 
  std::auto_ptr<TH1F> AllHitRCL1;

  std::auto_ptr<TH1F> PSimHitXD1;
  std::auto_ptr<TH1F> PSimHitYD1;
  
  std::auto_ptr<TH1F> Plside0PSimHitYD0; 
  std::auto_ptr<TH1F> Plside1PSimHitYD0;
  std::auto_ptr<TH1F> Plside0PSimHitXD0; 
  std::auto_ptr<TH1F> Plside1PSimHitXD0;
  std::auto_ptr<TH1F> NumEventNonCorruptedUtilized;

  std::auto_ptr<TH1D> Count_t2trackVectorALL_H0; 
  std::auto_ptr<TH1D> Count_t2trackVectorALL_H1;
  std::auto_ptr<TH1D> Count_t2trackVectorALL_H2; 
  std::auto_ptr<TH1D> Count_t2trackVectorALL_H3;


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
  
  std::map<unsigned int,double> mapsymbtoZ; 
  unsigned int matrixentries;
  std::vector<unsigned int> Testedcamera;
  std::vector<TrackInfo> MatriceTracce;
  double MaxDphi;
  double maxdphihit;
  double maxdrhit; 
  double Effmaxdphihit; 
  double Effmaxdrhit; 












  // ALGNMENT PART
  std::vector<std::vector<T2Hit> > roadXfit;
  std::vector<std::vector<T2Hit> > roadYfit;
  std::vector<std::vector<T2Hit> > origroadXfit;
  std::vector<std::vector<T2Hit> > origroadYfit;
  std::vector<uint32_t> alldetid;
  std::vector<std::vector<double> > startingtracksDX;
  std::vector<std::vector<double> > originaltracksDX;
  std::vector<std::vector<double> > startingtracksDY;
  std::vector<std::vector<double> > originaltracksDY;

  std::auto_ptr<TProfile> DetIdvsXShOnlyXfit;
  std::auto_ptr<TProfile> DetIdvsYShOnlyYfit;
  std::auto_ptr<TH1F> AllDetId;
  std::auto_ptr<TH1F> AllDetIdPlane;
  std::auto_ptr<TH1F> AllDetIdPlaneS;
  std::auto_ptr<TH1F> TrkChi2X;
  std::auto_ptr<TH1F> TrkChi2Y;
  std::vector<double> alldetZ;
  unsigned int HitNumb4Align;
  double MeasuredXYResol; 
  double  SHIFTprescale;
  unsigned int MaxStepalignstep;
   unsigned int Idreferencedet;
  double centrograv;
  double  meanz;
  double AlignmentHitRMax;
  //bool DoALign;
  unsigned int UseJointProb;
  double FitgravcenterZ;
  unsigned int  useRZforResol;
  //END ALIGNMENT PART

  double PhiMinForNoiseStudies;
  double PhiMaxForNoiseStudies;
  unsigned int DetForNoiseStudies; 

  std::vector<unsigned int> VectDeadSect_Plane;
  std::vector<unsigned int> VectDeadSect_Sector;
  
  std::auto_ptr<TH1F> C12AllDetId;
  std::auto_ptr<TH1F> C12AllDetIdPlane;
  std::auto_ptr<TH1F> C12AllDetIdPlaneS;



  //bool MakeVfatCumulativePlots;
  std::auto_ptr<TH2F> BunchStatusBit;

  unsigned int SelectedHalf;
  unsigned int NumHitGood ;
  unsigned int numevent;
  double Chicut;
  double energy;
  double DZScale;
  double tracketamin;
  double tracketamax;
  double PhiChiProbCut;
  double RChiProbCut;
  double Rwind;
  double Phiwind;
  int numhittrkeff;
  double chiRredCut,chiPhiredCut;
 
  bool ExcludeNoisyplane;
  unsigned int CommonNoiseClSize;

  unsigned int MaxPadAllowedInQuarter;
  unsigned int MaxPad, MaxStrip;
  unsigned int EffMaxPad, EffMaxStrip;
  double NoiseDphiMAX;
  double NoiseDrMAX;
  double AllowedDRTrackDistance;
  unsigned int  MaxEvents;
  unsigned int Effgoodhitnumber;
  unsigned int MaxTrkInQuarter;
  unsigned int MinTrkInQuarter;
  std::string outputFileName, CluLabel, HitLabel, RoadLabel,TrackLabel;
  std::string xmlfilenameFull;
  std::string xmlfilenameUsed_NotDead;
  std::string DeadSectFileName;
 

  bool produceVfatEffiFile;
  bool simufile;
  bool verbosity; 

  bool OnlycorruptionAnalysis;
  bool OnlyClusterAnalysis;
  bool  VFATMonitoring;
  bool skipSelectedEvents;
  std::string skipEventFileName;
  bool DispVtx;
  bool requiregoodChi;
  bool LookToRawEvent;
  bool UseUncorrupetdEventMap; //made corruption analysis
 std::vector<unsigned int> bunchesToAnalyse;
 
 
};

#endif // _TotemT1T2ValidationT2RecoValidationT2RecoValidation_H_

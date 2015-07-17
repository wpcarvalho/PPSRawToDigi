// -*- C++ -*-
//
// Original Author:  Mirko Berretti
//         Created:  Wed Dec 12 18:19:44 CET 2007
// $Id: T2RecoValidation.h,v 1.4 2008/11/17 14:50:08 berretti Exp $
//

#ifndef _TotemT1T2ValidationT2GeometryValidationT2RecoValidation_H_
#define _TotemT1T2ValidationT2GeometryValidationT2RecoValidation_H_

// system include files
#include <memory>
#include <math.h>
#include <TMath.h>
// user include files

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include "TCanvas.h"
#include "TStyle.h"
#include "TProfile.h"

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
#include "DataFormats/T2Cluster/interface/T2ROGeometry.h"
#include <boost/shared_ptr.hpp>
#include "TotemT1T2Validation/T2GeometryValidation/interface/T2DetClustReconsT.h"
#include "DataFormats/T2Digi/interface/T2PadDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2StripDigiCollection.h"
//
// class decleration
//




class T2GeometryAnalyzer : public edm::EDAnalyzer {
public:
  explicit T2GeometryAnalyzer(const edm::ParameterSet&);
  ~T2GeometryAnalyzer();
 
  T2DetClustReconsT theT2Clusterizer; 
  auto_ptr<T2PadClusterCollection> thePadClusters;
  auto_ptr<T2StripClusterCollection> theStripClusters;
  void InfoForVfatController(int* CCuAddrInt, std::string* Fec,  std::string* CCuAddr, std::string* I2CChan, int absvfatiid);
  void MakeEvCluster(std::auto_ptr<T2PadDigiCollection> PadDigiptr,std::auto_ptr<T2StripDigiCollection> StripDigiptr);
  
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
private:
  edm::InputTag T2PadDigiCollectionLabel;
  edm::InputTag T2StripDigiCollectionLabel;
  std::auto_ptr<TH1F> NumhitinTrack;
  std::auto_ptr<TH1F> NumhitinTrackGood;

  std::auto_ptr<TH1F> clusterstripentries;
  std::auto_ptr<TH1F> clusterpadentries;
  std::auto_ptr<TH1F> diffphiCluGun;
  std::auto_ptr<TH1F> RLocHRecoH;
  std::auto_ptr<TH1F> Trketa;
  std::auto_ptr<TH1F> Trkphi;
  std::auto_ptr<TH1F> DPhiGoodTrk;
  std::auto_ptr<TH1F> DEtaGoodTrk;
  std::auto_ptr<TH1F> diffphiGUNHIT;
  std::auto_ptr<TH1F> diffphiCluHit;
  std::auto_ptr<TH1F> Trkphigood;
  std::auto_ptr<TH1F> Trketagood;
  std::auto_ptr<TH1F> diffRCluHit;
  std::auto_ptr<TH1F> Chi2RProb;
  std::auto_ptr<TH1F> Chi2PhiProb;
  std::auto_ptr<TH1F> TrKPhi0degstudyg;
  std::auto_ptr<TH1F> SourceofSecondary;
  std::auto_ptr<TH1F> SourceofReco;
  std::auto_ptr<TH1F> Z0distro;
  std::auto_ptr<TH1F> R0distro;
  std::auto_ptr<TH1F> tantheta;
  std::auto_ptr<TH1F> tanthetam0;
  std::auto_ptr<TH1F> RecoZHit;
  std::auto_ptr<TH1F> SimuZHit;


  std::auto_ptr<TProfile>  MeanT2PadvsEta; 
  std::auto_ptr<TProfile>  MeanT2StripvsEta;
  std::auto_ptr<TProfile>  MeanT2RoadEntriesvsEta;
  std::auto_ptr<TProfile>  MeanT2HitCl1vsEta;
  std::auto_ptr<TProfile>  SingleParticleEfficiency;
  std::auto_ptr<TProfile>  SingleParticleEfficiencyAllCuts;
  std::auto_ptr<TProfile>  SingleParticleEfficiencyCutEta;
  std::auto_ptr<TProfile>  SingleParticleEfficiencyCutChi;
  std::auto_ptr<TProfile>  SingleParticleEfficiencyCutDZ;
  std::auto_ptr<TProfile>  MultiParticleEfficiencyCut;
  std::auto_ptr<TProfile>  MultiParticleEfficiencyCutNorm;
  std::auto_ptr<TProfile>  MeanT2HitvsEta;
  std::auto_ptr<TProfile>  MeanT2GeantHitvsEta;
  std::auto_ptr<TProfile>  P0NumEtacorr;
  std::auto_ptr<TProfile>  P0NumEnergycorr;

   std::auto_ptr<TH1F> TrkEtaGoodH0;
   std::auto_ptr<TH1F> TrkEtaGoodH1;
     std::auto_ptr<TH1F> TrkEtaH0;
   std::auto_ptr<TH1F> TrkEtaH1;

  // std::auto_ptr<TProfile>   SingleParticleEfficiencyAllCutsH0; 
  // std::auto_ptr<TProfile>   SingleParticleEfficiencyAllCutsH1;
  std::auto_ptr<TH1F> EtaParticle;  
  std::auto_ptr<TH1F> EtaParticleAll;  
  std::auto_ptr<TProfile>  MeanEtaResvsPhiHem2; 
  std::auto_ptr<TProfile>  MeanEtaResvsPhiHem1; 
  std::auto_ptr<TProfile>  MeanEtaResvsPhi;

  std::auto_ptr<TProfile>  MeanEtaResvsPhiHem2Part; 
  std::auto_ptr<TProfile>  MeanEtaResvsPhiHem1Part; 
  std::auto_ptr<TProfile>  MeanEtaResvsPhiPart;
  std::auto_ptr<TProfile>  P0EtaEnergycorr;
  std::auto_ptr<TProfile>  MeanEtaResvsPhiCut;
  std::auto_ptr<TProfile>  MeanEtaResvsPhiHem2Cut; 
  std::auto_ptr<TProfile>  MeanEtaResvsPhiHem1Cut;
  std::auto_ptr<TProfile>   MeanT2StripDigivsEta;
  std::auto_ptr<TProfile>   MeanT2PadDigivsEta;

  std::auto_ptr<TProfile>  NumChPartVsNumP0; 
  std::auto_ptr<TProfile>  NumP0VsNumTracks;
  std::auto_ptr<TProfile>  NumChPartInVsNumChPartOut;

  std::auto_ptr<TH1F> HNumrecotrackcutmag0;
  std::auto_ptr<TH1F> HTrackcounter;
  std::auto_ptr<TH1F> chPartinT2;

  std::auto_ptr<TH2F>  P0EtaEnergy;
  std::auto_ptr<TH2F>  muEtaEnergy;
  std::auto_ptr<TH2F>  ChOutEtaEnergy;
  std::auto_ptr<TH1F>  ChEnergyinT2;
  std::auto_ptr<TH1F> HTrackcounterNonCrit;
  std::auto_ptr<TH1F> HTrackcounterCrit;
  std::auto_ptr<TH1F> stablepdg;
  std::auto_ptr<TCanvas> Chi2PhiProbLogy;
  std::auto_ptr<TCanvas> Chi2RProbLogy;

  std::auto_ptr<TH2F>  PimenoEtaEnergy[15];

  std::auto_ptr<TH1F> chiPhidistro;
  std::auto_ptr<TH1F> chiRdistro;
 
  std::auto_ptr<TH1F> SymbIdHT0;
  std::auto_ptr<TH1F> Trketachicut;
  std::auto_ptr<TH1F> Trketaetachicut; 
  std::auto_ptr<TH1F> Trketaetacut;
  std::auto_ptr<TH1F> Trketadzcut;
  std::auto_ptr<TH1F> Trketadzetacut;
  std::auto_ptr<TH1F> Trketadzetachicut;
  std::auto_ptr<TH1F> DPhiChiCutOneTrk;
  std::auto_ptr<TH1F> DEtaChiCutOneTrk;
  std::auto_ptr<TH1F> MinDz2Plane;
  std::auto_ptr<TH1F>  energyP04575;
  std::auto_ptr<TH1F>  etaP04575;
  bool singleparticle;

  unsigned int numevent;
  double Chicut;
  double energy;
  double DZScale;
  double tracketamin;
  double tracketamax;
  double PhiChiProbCut;
  double RChiProbCut;
  //bool onlyMultiParticle;
  bool DoClusterFromVFatChannel;
  bool PrintOnlyFinalOutput;
  bool PrintChannelNumerationAsInVfatController;
  bool OnlyVfatANDchannel_Finalprinting;
  std::string outputFileName, HepMCProductLabel, CluLabel, HitLabel, RoadLabel,TrackLabel;

};

#endif

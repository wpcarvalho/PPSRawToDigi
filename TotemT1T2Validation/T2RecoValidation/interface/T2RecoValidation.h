// -*- C++ -*-
//
// Original Author:  Mirko Berretti
//         Created:  Wed Dec 12 18:19:44 CET 2007
// $Id: T2RecoValidation.h,v 1.6.2.2 2009/08/07 16:02:41 zhangz Exp $
//

#ifndef _TotemT1T2ValidationT2RecoValidationT2RecoValidation_H_
#define _TotemT1T2ValidationT2RecoValidationT2RecoValidation_H_

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
#include <boost/shared_ptr.hpp>
//
// class decleration
//

class T2RecoAnalyzer : public edm::EDAnalyzer {
public:
  explicit T2RecoAnalyzer(const edm::ParameterSet&);
  ~T2RecoAnalyzer();
  unsigned int RawtoSymb(uint32_t thedet);

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
edm::InputTag SimTrackContainerLabel;
std::auto_ptr<TH1D> NumRecoTrkPerEvent;
 std::auto_ptr<TH1D> NumRecoTrkWhenAtLeastAstableInT2;
  std::auto_ptr<TH2D> EventEffiXChpartInT2YNumTrk;
std::auto_ptr<TH1F> TrackHemisphereWhenChargedParticle;
std::auto_ptr<TH1F> RoadHemisphereNorequirement;

std::auto_ptr<TH1F> TrackHemisphereNorequirement;
std::auto_ptr<TProfile>  PadCluMultiplVsPlane;
std::auto_ptr<TProfile> StripCluMultiplVsPlane;
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
  std::auto_ptr<TH1D> NumTotRoadH0;
  std::auto_ptr<TH1D> NumTotRoadH1;
  std::auto_ptr<TH1D> NumTotRoadH2;
  std::auto_ptr<TH1D> NumTotRoadH3;
  std::auto_ptr<TH1F>    NonTriggeredEnergyPlus;
std::auto_ptr<TH1D> EtaMaxAroundLeftProton;
   std::auto_ptr<TH1D> NumEvtSel2Proton;
   std::auto_ptr<TH1D> AllRecoMassHist;
   std::auto_ptr<TH1D> forwInvMassHist;
   std::auto_ptr<TH1D> CentInvMassHist;
   std::auto_ptr<TH1D> PomeronInvMassHist;
   std::auto_ptr<TH1D> CentralEnergyWhenHigh005Csi;
   std::auto_ptr<TH1F> NonTriggeredEnergyInPlusOnly_When100GeVDiscr;
   std::auto_ptr<TH1F>  NonTriggeredEnergyInPlusAndMinus_When100GeVDiscr;




  std::auto_ptr<TH1D>     Pom_minus_Cent_003_008;
   std::auto_ptr<TH1D>   CentralEnergyWhenfscT1T2OFF_An005Cut;
   std::auto_ptr<TH1D>   ForwEnergyWhenForwOFF_005Cut_Discr100GeV;
   std::auto_ptr<TH1D>   ForwEnergyWhenForwOFF_005Cut_Discr200GeV;
   std::auto_ptr<TH1D>   ForwEnergyWhenForwOFF_005Cut_Discr300GeV;
   std::auto_ptr<TH1D>   ForwEnergyWhenForwOFF_005Cut_Discr400GeV;
   std::auto_ptr<TH1D>   ForwEnergyWhenForwOFF_003_008_Cut_Discr100GeV;
   std::auto_ptr<TH1D>   ForwEnergyWhenForwOFF_003_008_Cut_Discr200GeV;
   std::auto_ptr<TH1D>   ForwEnergyWhenForwOFF_003_008_Cut_Discr300GeV;
   std::auto_ptr<TH1D>   ForwEnergyWhenForwOFF_003_008_Cut_Discr400GeV;
   std::auto_ptr<TH1D>   ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC100GeV;
   std::auto_ptr<TH1D>   ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC200GeV;
   std::auto_ptr<TH1D>   ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC300GeV;
   std::auto_ptr<TH1D>   ForwEnergyWhenForwOFF_003_008_Cut_Discr_NoFSC400GeV;







std::auto_ptr<TH1D> ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr100GeV;
std::auto_ptr<TH1D> ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr200GeV;
std::auto_ptr<TH1D> ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr300GeV;
std::auto_ptr<TH1D>  ForwEnergyWhenForwOFFCMSOFF_005Cut_Discr400GeV;


  std::auto_ptr<TH1D> T1chargeWhen_CentralEnergyWhenHigh005Csi_30;
  std::auto_ptr<TH1D> FSCchargeWhen_CentralEnergyWhenHigh005Csi_30;  
  std::auto_ptr<TH1D> T2TrkWhen_CentralEnergyWhenHigh005Csi_30;
  std::auto_ptr<TH1D> T1chargeWhen_CentralEnergyWhenHigh005Csi_100;
  std::auto_ptr<TH1D> FSCchargeWhen_CentralEnergyWhenHigh005Csi_100;  
  std::auto_ptr<TH1D> T2TrkWhen_CentralEnergyWhenHigh005Csi_100;  

  std::auto_ptr<TH1D> ForwardWhen_CentralEnergyWhenHigh005Csi_30;  
  std::auto_ptr<TH1D> ForwardWhen_CentralEnergyWhenHigh005Csi_100;  

  std::auto_ptr<TH2D>  csi1_vs_csi2;
  std::auto_ptr<TH1D> ForwRegionMassWhenT2OFF_Discr100GeV;
  std::auto_ptr<TH1D> ForwRegionMassWhenT2OFF_Discr200GeV;
  std::auto_ptr<TH1D> ForwRegionMassWhenT2OFF_Discr300GeV;
  std::auto_ptr<TH1D> ForwRegionMassWhenT2OFF_Discr400GeV;

  std::auto_ptr<TH1D> ForwRegionMassWhenT2OFF_Discr100GeVNolargeEta;
  std::auto_ptr<TH1D> ForwRegionMassWhenT2OFF_Discr200GeVNolargeEta;
  std::auto_ptr<TH1D> ForwRegionMassWhenT2OFF_Discr300GeVNolargeEta;
  std::auto_ptr<TH1D> ForwRegionMassWhenT2OFF_Discr400GeVNolargeEta;

   std::auto_ptr<TH1D> NonCentralMassFromPom; std::auto_ptr<TH1D> csi1Histo;
  std::auto_ptr<TProfile> EfficiencyTrigT2PlusVSEplus;
	
   std::auto_ptr<TH1D>  LeakT1T2;
   std::auto_ptr<TH1D>  NonTriggeredEnergy;
 
 std::auto_ptr<TH1D> CumulativeEtaGeneratorStableParticles_WhenACHARGEInT2;
 std::auto_ptr<TH1D> CumulativeEtaGeneratorStableParticles;
   std::auto_ptr<TH1D> CumulativeEtaGeneratorStableParticlesRightOnly;
   std::auto_ptr<TH1D> CumulativeEtaGeneratorStableParticlesLeftOnly;
 std::auto_ptr<TH1D> CumulativeEtaGeneratorStableCHParticles;
 std::auto_ptr<TH1D> CumulativeEtaGeneratorALLParticles;
  std::auto_ptr<TH1D> TrackPlusWhenGenRapGapPlus;
  std::auto_ptr<TH1D>   GeneratorEvtHemisp;
  std::auto_ptr<TH1D>  ProbGetLowMass;
  std::auto_ptr<TH1D>  GeneratorEvtHemispCH;
  std::auto_ptr<TH1D>  MultiplicityWhenT2ShouldHAVE;
  std::auto_ptr<TH1D>  DifferenceMassFromProtons_andAllOtherParticles;

 
  std::auto_ptr<TH1D> CentralEnergyWhenfscT1T2OFF;
  std::auto_ptr<TH2D> PomeronInvMassWhenfscT1T2OFFvsCentralE;
  std::auto_ptr<TH1D> CentralEnergyWhenHigh005Csi_AndfscT1T2OFF;
  std::auto_ptr<TProfile>  SingleParticleEfficiency;
  std::auto_ptr<TProfile>  SingleParticleEfficiencyCut;
  std::auto_ptr<TProfile>  MultiParticleEfficiencyCut;
  std::auto_ptr<TProfile>  MultiParticleEfficiencyCutNorm;

  std::auto_ptr<TCanvas> Chi2PhiProbLogy;
  std::auto_ptr<TCanvas> Chi2RProbLogy;
  std::auto_ptr<TH1D>  ProtonRapidity;
  unsigned int numevent;
  double Chicut;
  double energy;
  double DZScale;
  double tracketamin;
  double tracketamax;
  double PhiChiProbCut;
  double RChiProbCut;
  bool singleparticle;  
 
  std::string outputFileName, HepMCProductLabel, CluLabel, HitLabel, RoadLabel,TrackLabel;

};

#endif // _TotemT1T2ValidationT2RecoValidationT2RecoValidation_H_

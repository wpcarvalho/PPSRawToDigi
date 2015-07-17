// -*- C++ -*-
//
// Original Author:  Mirko Berretti
//         Created:  Wed Dec 12 18:19:44 CET 2007
// $Id: T2GeneratorAnalyzer.h,v 1.6.2.2 2009/08/07 16:02:41 zhangz Exp $
//T2GeneratorAnalyzer

#ifndef _TotemT1T2ValidationGeneratorAnalysisT2GeneratorAnalyzer_H_
#define _TotemT1T2ValidationGeneratorAnalysisT2GeneratorAnalyzer_H_

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


#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"




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

class T2GeneratorAnalyzer : public edm::EDAnalyzer {
public:
  explicit T2GeneratorAnalyzer(const edm::ParameterSet&);
  ~T2GeneratorAnalyzer();


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

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
  std::auto_ptr<TH1F>   EnergyPi0InT2;
  std::auto_ptr<TH1F>   StatusPi0InT2;
  std::auto_ptr<TH1F>     EnergyPiMenoPlusInT2;
   std::auto_ptr<TH2F>   P0EtaEnergy;
  std::auto_ptr<TH2F>   PmenoEtaEnergy2;
  std::auto_ptr<TProfile>  SingleParticleEfficiency;
  std::auto_ptr<TProfile>  SingleParticleEfficiencyCut;
  std::auto_ptr<TProfile>  MultiParticleEfficiencyCut;
  std::auto_ptr<TProfile>  MultiParticleEfficiencyCutNorm;
  std::auto_ptr<TProfile>  DNDetaMBALLMCGenerator_GeneratorTriggered;
  std::auto_ptr<TProfile>  DNDetaMBALLMCGenerator_NoEvtRequirement;
  std::auto_ptr<TProfile>  DNDetaMBALLMCGenerator_LHCbRequirement;	
  std::auto_ptr<TProfile>  DNDetaMBALLMCGenerator_T2Requirement;
  std::auto_ptr<TProfile>  DNDetaMBALLMCGenerator_K0s_NoEvtRequirement;
  std::auto_ptr<TProfile>   DNDE_gamma_NoEvtRequirement;
  std::auto_ptr<TProfile>  DNDetaMBALLMCGenerator_K0s_CMS;
  std::auto_ptr<TProfile>  DNDetaMBALLMCGenerator_CMS;
  std::auto_ptr<TProfile>  DNDetaMBALLMCGenerator_ALICE;
  std::auto_ptr<TProfile>  EtaParticleVsNumTrk;

  std::auto_ptr<TH1D> StableChPartilceEnergy;
  std::auto_ptr<TProfile> Eta2_vs_Energy;

  
  std::auto_ptr<TCanvas> Chi2PhiProbLogy;
  std::auto_ptr<TCanvas> Chi2RProbLogy;



std::auto_ptr<TProfile>  NumPadCluVsPlaneAll3H0;
 std::auto_ptr<TProfile>  PadCluSizeVsPlaneAll3H0;
std::auto_ptr<TProfile>  NumPadCluVsPlaneAll3H1;
std::auto_ptr<TProfile>   PadCluSizeVsPlaneAll3H1;
std::auto_ptr<TProfile>  NumPadCluVsPlaneAll3H2;
std::auto_ptr<TProfile>   PadCluSizeVsPlaneAll3H2;
std::auto_ptr<TProfile>  NumPadCluVsPlaneAll3H3;
std::auto_ptr<TProfile>   PadCluSizeVsPlaneAll3H3;
std::auto_ptr<TProfile>  NumStripCluVsPlaneAll3H0;
 std::auto_ptr<TProfile>  StripCluSizeVsPlaneAll3H0;
std::auto_ptr<TProfile>  NumStripCluVsPlaneAll3H1;
std::auto_ptr<TProfile>   StripCluSizeVsPlaneAll3H1;
std::auto_ptr<TProfile>  NumStripCluVsPlaneAll3H2;
std::auto_ptr<TProfile>   StripCluSizeVsPlaneAll3H2;
std::auto_ptr<TProfile>  NumStripCluVsPlaneAll3H3;
std::auto_ptr<TProfile>   StripCluSizeVsPlaneAll3H3;







  std::vector<unsigned int> vectorMBALLMCGenerator_gammaEInLhcf;
  std::vector<unsigned int> vectorMBALLMCGenerator;
  std::vector<unsigned int> vectorMBALLMCGenerator_K0s;
  std::vector<unsigned int> vectorMBALLMCGeneratorPtCut;
  unsigned int numevent;
  double Chicut;
  double energy;
  double DZScale;
  double tracketamin;
  double tracketamax;
  double PhiChiProbCut;
  double RChiProbCut;
  bool singleparticle;  

  unsigned int totnumberbinfordndeta;
  double maxetafordndetahisto;
  double etabinsize;
  bool FullInfos;
  std::string outputFileName, HepMCProductLabel, CluLabel, HitLabel, RoadLabel,TrackLabel;

};

#endif // _TotemT1T2ValidationT2GeneratorAnalyzerT2GeneratorAnalyzer_H_

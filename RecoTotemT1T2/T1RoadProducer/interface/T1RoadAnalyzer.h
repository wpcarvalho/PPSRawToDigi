/*
  Created by Fabrizio Ferro - INFN Genova for TOTEM
*/
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/T1Road/interface/T1Road.h"
#include "DataFormats/T1Road/interface/T1RecHitGlobal.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/Provenance/interface/Provenance.h"
#include "DataFormats/Provenance/interface/BranchDescription.h"

#include "Geometry/TotemGeometry/interface/T1Geometry.h"
//
// class decleration
//

#include "TROOT.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

//class TFile;
//class TH1D;

class T1RoadAnalyzer : public edm::EDAnalyzer 
{

  enum {maxHits = 10000};  

 public:
  explicit T1RoadAnalyzer(const edm::ParameterSet&);
  ~T1RoadAnalyzer();
  float Eta(float x,float y,float z);
  float Phi(float, float);

 private:
  edm::InputTag t1RecHit2DCollectionMyRecoCollLabel;
  edm::InputTag t1RoadCollectionLabel;
  edm::InputTag simVertexContainerLabel;
  edm::InputTag simTrackContainerLabel;
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;


  TFile* theFile; 

  TTree* tree;

  unsigned int nHits;
  double rx[maxHits],  ry[maxHits],  rz[maxHits], ee[maxHits], ff[maxHits];
  int evto[maxHits];

  TH1D * hDeltaEta;
  TH1D * hDeltaPhi;
  TH1D * hDeltaEtaPri;
  TH1D * hDeltaPhiPri;

  TH1D * hRoadNumber;
  TH1D * hRoadSize;
  TH1D * hRoadPlaneSize;

  TH1D * hPrimaryRoadNumber;
  TH1D * hPrimaryTracks;
  TH1D * hPrimaryRoadTrackRatio;
  TH1D * hRoadTrackRatio;
  TH1D * hPrimaryRoadRoadRatio;
  TH1D * hLostPrimaryTracks;
  TH1D * hLostPrimaryTrackRatio;

  T1Geometry *layer;
  // ----------member data ---------------------------
};

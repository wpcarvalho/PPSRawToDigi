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
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"

//class TFile;
//class TH1D;

class T1RecHitAnalyzer : public edm::EDAnalyzer 
{

  enum {maxHits = 10000};  

 public:
  explicit T1RecHitAnalyzer(const edm::ParameterSet&);
  ~T1RecHitAnalyzer();
  float Eta(float x,float y,float z);
  float Phi(float, float);

  void writeHistograms();
 private:
  edm::InputTag t1RecHit2DCollectionLabel;
  edm::InputTag t1DigiWireCollectionLabel;
  edm::InputTag t1ClusterCollectionLabel;
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();


  TFile* theFile; 

  std::string outputFileName;


  std::auto_ptr<TH2>hCSC000;
  std::auto_ptr<TH2>hCSC004;
  std::auto_ptr<TH2>hCSC005;

  std::auto_ptr<TH2>hCSC010;
  std::auto_ptr<TH2>hCSC014;
  std::auto_ptr<TH2>hCSC015;

  std::auto_ptr<TH2>hCSC020;
  std::auto_ptr<TH2>hCSC024;
  std::auto_ptr<TH2>hCSC025;

  std::auto_ptr<TH2>hCSC030;
  std::auto_ptr<TH2>hCSC034;
  std::auto_ptr<TH2>hCSC035;

  std::auto_ptr<TH2>hCSC040;
  std::auto_ptr<TH2>hCSC044;
  std::auto_ptr<TH2>hCSC045;

  std::auto_ptr<TH2>hCSC001;
  std::auto_ptr<TH2>hCSC002;
  std::auto_ptr<TH2>hCSC003;

  std::auto_ptr<TH2>hCSC011;
  std::auto_ptr<TH2>hCSC012;
  std::auto_ptr<TH2>hCSC013;

  std::auto_ptr<TH2>hCSC021;
  std::auto_ptr<TH2>hCSC022;
  std::auto_ptr<TH2>hCSC023;

  std::auto_ptr<TH2>hCSC031;
  std::auto_ptr<TH2>hCSC032;
  std::auto_ptr<TH2>hCSC033;

  std::auto_ptr<TH2>hCSC041;
  std::auto_ptr<TH2>hCSC042;
  std::auto_ptr<TH2>hCSC043;

//=================================================


  std::auto_ptr<TH1>hCSCcsc000;
  std::auto_ptr<TH1>hCSCcsc004;
  std::auto_ptr<TH1>hCSCcsc005;

  std::auto_ptr<TH1>hCSCcsc010;
  std::auto_ptr<TH1>hCSCcsc014;
  std::auto_ptr<TH1>hCSCcsc015;

  std::auto_ptr<TH1>hCSCcsc020;
  std::auto_ptr<TH1>hCSCcsc024;
  std::auto_ptr<TH1>hCSCcsc025;

  std::auto_ptr<TH1>hCSCcsc030;
  std::auto_ptr<TH1>hCSCcsc034;
  std::auto_ptr<TH1>hCSCcsc035;

  std::auto_ptr<TH1>hCSCcsc040;
  std::auto_ptr<TH1>hCSCcsc044;
  std::auto_ptr<TH1>hCSCcsc045;

  std::auto_ptr<TH1>hCSCcsc001;
  std::auto_ptr<TH1>hCSCcsc002;
  std::auto_ptr<TH1>hCSCcsc003;

  std::auto_ptr<TH1>hCSCcsc011;
  std::auto_ptr<TH1>hCSCcsc012;
  std::auto_ptr<TH1>hCSCcsc013;

  std::auto_ptr<TH1>hCSCcsc021;
  std::auto_ptr<TH1>hCSCcsc022;
  std::auto_ptr<TH1>hCSCcsc023;

  std::auto_ptr<TH1>hCSCcsc031;
  std::auto_ptr<TH1>hCSCcsc032;
  std::auto_ptr<TH1>hCSCcsc033;

  std::auto_ptr<TH1>hCSCcsc041;
  std::auto_ptr<TH1>hCSCcsc042;
  std::auto_ptr<TH1>hCSCcsc043;






  T1Geometry *layer;
  // ----------member data ---------------------------
};

// -*- C++ -*-
//
// Package:    T1RecHitAnalyzer2
// Class:      T1RecHitAnalyzer2
// 
/**\class T1RecHitAnalyzer2 T1RecHitAnalyzer2.cc RecoTotemT1T2/T1RecHitAnalyzer2/src/T1RecHitAnalyzer2.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Enrico Robutti
//         Created:  Wed Mar 28 10:38:39 CEST 2012
// $Id$
//
//


// system include files
#include <memory>
#include <string>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/T1RecHit/interface/T1RecHit2DCollection.h"
#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "TH2F.h"
#include "TFile.h"

using namespace std;
using namespace edm;

//
// class declaration
//

class T1RecHitAnalyzer2 : public edm::EDAnalyzer {
public:
  explicit T1RecHitAnalyzer2(const edm::ParameterSet&);
  ~T1RecHitAnalyzer2();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  // ----------member data ---------------------------
  unsigned int _verbosity;
  string _t1RecHitModuleLabel;
  string _t1RecHitProductLabel;
  string _outputHistoFileName;
  unsigned int _printProgressFrequency;
  TH2F* _hHitOccupancy[2][5][6];
  double _binSize;
  int _evtCount;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
T1RecHitAnalyzer2::T1RecHitAnalyzer2(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getUntrackedParameter<unsigned int>("verbosity")),
  _t1RecHitModuleLabel(iConfig.getParameter<string>("t1RecHitModuleLabel")),
  _t1RecHitProductLabel(iConfig.getParameter<string>("t1RecHitProductLabel")),
  _outputHistoFileName(iConfig.getUntrackedParameter<string>("outputHistoFileName")),
  _printProgressFrequency(iConfig.getUntrackedParameter<unsigned int>("printProgressFrequency")),
  _binSize(iConfig.getUntrackedParameter<double>("binSize", 10.))
{
  //now do what ever initialization is needed

}


T1RecHitAnalyzer2::~T1RecHitAnalyzer2()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  for (int iArm = 0; iArm < 2; iArm++)
    for (int iLyr = 0; iLyr < 5; iLyr++)
      for (int iSext = 0; iSext < 6; iSext++)
	delete _hHitOccupancy[iArm][iLyr][iSext];
}


//
// member functions
//

// ------------ method called for each event  ------------
void
T1RecHitAnalyzer2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Print progress frequency
  _evtCount++;
  if (_evtCount%_printProgressFrequency == 0)
    cout << "Event n. " << _evtCount << endl;

  // Get reconstructed hits
  Handle<T1RecHit2DCollection> hit2DCollection;
  iEvent.getByLabel(_t1RecHitModuleLabel, _t1RecHitProductLabel, hit2DCollection);
  for (T1RecHit2DCollection::const_iterator itH = hit2DCollection->begin(); itH != hit2DCollection->end(); itH++) {
    const T1DetId& csc = (*itH).t1DetId();
    int arm = csc.Arm();
    int layer = csc.Plane();
    int sextant = csc.CSC();
    _hHitOccupancy[arm][layer][sextant]->Fill((*itH).localPosition().x(), (*itH).localPosition().y());
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
T1RecHitAnalyzer2::beginJob()
{
  _evtCount = -1;
  // Book histograms
  for (int iArm = 0; iArm < 2; iArm++)
    for (int iLyr = 0; iLyr < 5; iLyr++)
      for (int iSext = 0; iSext < 6; iSext++) {
	TString hName = TString("hHitOcc_") + (unsigned long)iArm + "_" +
	                (unsigned long)iLyr + "_" + (unsigned long)iSext;
	TString hTitle = TString("Hit occupancy, Arm ") + ((iArm == 0) ? "+" : "-") +
	                 ", Layer " + (unsigned long)iLyr + ", Sextant " + (ULong_t)iSext;
	unsigned int nBins = 1600./_binSize;
	_hHitOccupancy[iArm][iLyr][iSext] = new TH2F(hName, hTitle, nBins, -800., 800., nBins, -800., 800.);
      }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
T1RecHitAnalyzer2::endJob() 
{
  // Write histograms to ROOT file
  TFile outFile(_outputHistoFileName.c_str(), "RECREATE");
  for (int iArm = 0; iArm < 2; iArm++)
    for (int iLyr = 0; iLyr < 5; iLyr++)
      for (int iSext = 0; iSext < 6; iSext++)
	_hHitOccupancy[iArm][iLyr][iSext]->Write();
  outFile.Close();
}

// ------------ method called when starting to processes a run  ------------
void 
T1RecHitAnalyzer2::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
T1RecHitAnalyzer2::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
T1RecHitAnalyzer2::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
T1RecHitAnalyzer2::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
T1RecHitAnalyzer2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(T1RecHitAnalyzer2);

// -*- C++ -*-
//
// Package:    T1Analysis2
// Class:      T1Analysis2
// 
/**\class T1Analysis2 T1Analysis2.cc TotemAnalysis/T1Analysis2/src/T1Analysis2.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  fabrizio ferro
//         Created:  Mon Nov 28 09:33:34 CET 2011
// $Id$
//
//


// system include files
#include <memory>
#include <map>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "TotemRawDataLibrary/DataFormats/interface/RawEvent.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/T1Road/interface/T1RecHitGlobal.h"
#include "DataFormats/T1T2Track/interface/T1T2TrackCollection.h"

#include <boost/shared_ptr.hpp>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

#define PI 3.141592654

//
// class declaration
//

class T1Zetas : public edm::EDAnalyzer {
public:
  explicit T1Zetas(const edm::ParameterSet&);
  ~T1Zetas();

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



  std::auto_ptr<TH1D> hZetaAtRmin;
  std::auto_ptr<TH1D> hZImpact;
  std::auto_ptr<TH2D> hZetas;



  TFile* theFile;

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
T1Zetas::T1Zetas(const edm::ParameterSet& iConfig)

{
//now do what ever initialization is needed

}


T1Zetas::~T1Zetas()
{
 
// do anything here that needs to be done at desctruction time
// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
T1Zetas::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;



//

// T1 data

  edm::Handle<T1T2TrackCollection> T1trackCollection;
  iEvent.getByLabel("t1tracks2","T1TrackColl",T1trackCollection);

  T1T2TrackCollection::const_iterator T1TC_it;

// tagli su numero di tracce e chiquadro? (50)

//  if(T1trackCollection->size() < 40)
  for(T1TC_it=T1trackCollection->begin(); T1TC_it!=T1trackCollection->end(); T1TC_it++){



    if( (*T1TC_it).Eta() > 0 
// && (*T1TC_it).Rmin()<100

 ){

    hZetaAtRmin->Fill( (*T1TC_it).Z_at_Rmin() );
    double C0  =    (*T1TC_it).GetHitT1(0).GlobalPosition().perp2()/ ( (*T1TC_it).GetTx()* (*T1TC_it).GetHitT1(0).GlobalPosition().x() + (*T1TC_it).GetTy()* (*T1TC_it).GetHitT1(0).GlobalPosition().y()  );
    double ZImpact = (*T1TC_it).GetHitT1(0).GlobalPosition().z() - C0;
/*
  (atrkGeant.GetHitT2(0).GetHitX()*atrkGeant.GetHitT2(0).GetHitX()+atrkGeant.GetHitT2(0).GetHitY()*atrkGeant.GetHitT2(0).GetHitY())/(atrkGeant.GetTx()*atrkGeant.GetHitT2(0).GetHitX()+atrkGeant.GetTy()*atrkGeant.GetHitT2(0).GetHitY());
  double Z0impact=atrkGeant.GetHitT2(0).GetHitZ()-C0Geant;
*/
    hZImpact->Fill( ZImpact );
    hZetas->Fill((*T1TC_it).Z_at_Rmin(), ZImpact);
    }


  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
T1Zetas::beginJob()
{


 
 

  hZetaAtRmin = std::auto_ptr<TH1D>(new TH1D("hZetaAtRmin","hZetaAtRmin",1000,-20000,20000));
  hZetaAtRmin->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hZetaAtRmin","hZetaAtRmin",100,-0.5,99.5)

  hZImpact = std::auto_ptr<TH1D>(new TH1D("hZImpact","hZImpact",1000,-20000,20000));
  hZImpact->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hZImpact","hDphi",100,-0.5,99.5)

  hZetas = std::auto_ptr<TH2D>(new TH2D("hZetas","hZetas",1000,-20000,20000,1000,-20000,20000));
  hZetas->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hZetas","hDphi",100,-0.5,99.5)

 




}
// ------------ method called once each job just after ending the event loop  ------------
void 
T1Zetas::endJob() 
{

  theFile = TFile::Open("T1Zetas.root", "recreate");
  if(!theFile || !theFile->IsWritable())
    {
      std::cout<<"Output file not opened correctly!!"<<std::endl;
    }

  hZetaAtRmin->Write(); 
  hZImpact->Write(); 
  hZetas->Write(); 






  theFile->Close();





}

// ------------ method called when starting to processes a run  ------------
void 
T1Zetas::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
T1Zetas::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
T1Zetas::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
T1Zetas::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
T1Zetas::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
//The following says we do not know what parameters are allowed so do no validation
// Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(T1Zetas);

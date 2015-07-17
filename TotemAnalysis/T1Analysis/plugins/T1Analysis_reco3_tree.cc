// -*- C++ -*-
//
// Package:    T1Analysis
// Class:      T1Analysis
// 
/**\class T1Analysis T1Analysis.cc TotemAnalysis/T1Analysis/src/T1Analysis.cc

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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "TotemRawDataLibrary/DataFormats/interface/RawEvent.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2D.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2DCollection.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/T1Road/interface/T1RecHitGlobal.h"
#include "DataFormats/T1T2Track/interface/T1T2TrackCollection.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"

#include <boost/shared_ptr.hpp>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

//
// class declaration
//

class T1Analysis_reco3_tree : public edm::EDAnalyzer {
public:
  explicit T1Analysis_reco3_tree(const edm::ParameterSet&);
  ~T1Analysis_reco3_tree();
  
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

  edm::InputTag rawEventLabel;
  edm::InputTag t1RecHit2DCollectionLabel;
  unsigned int NSelectEvent;
  unsigned int NAtLeastOneTrackEvent;
  unsigned int NEventTot;

  std::auto_ptr<T1Geometry> layer;
  TFile* theFile;


  int Verbosity;
  std::string outputFileName;
  bool DumpOnly0TrackEvents;
  int MaxTracks;
  double ChiProbCut;
  std::string trackLabel;
  double EtaMin, EtaMax;


  std::auto_ptr<TH1D> hZatRmin;
  std::auto_ptr<TH1D> hRmin;
  std::auto_ptr<TH1D> hZImpact;
  std::auto_ptr<TH2D> hZetas;
  std::auto_ptr<TH2D> hZR;
  std::auto_ptr<TH1D> hChiOverN;
  std::auto_ptr<TH1D> hChiProb;

 std::auto_ptr<TTree> _tree;
  unsigned int _event;
  unsigned int _tracks;
  unsigned int _recoHits;
  vector<double>_ZatRmin;
  vector<double>_ZImpact;
  vector<double>_Rmin;
  vector<double>_Deta;
  vector<double>_Dphi;
  vector<double>_ChiProb;
  vector<int>_Hits;
  vector<double>_Eta;
  vector<double>_Phi;
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
T1Analysis_reco3_tree::T1Analysis_reco3_tree(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  rawEventLabel = iConfig.getParameter<edm::InputTag>("RawEventLabel");
  t1RecHit2DCollectionLabel = iConfig.getParameter<edm::InputTag>("T1RecHit2DCollectionLabel");

  std::string thegeometryfile = "Geometry/TotemGeometry/data/T1_data_geometry.dat";
  layer = std::auto_ptr<T1Geometry>(new T1Geometry(thegeometryfile));

  Verbosity = iConfig.getParameter<int>("Verbosity");
  outputFileName = iConfig.getParameter<std::string>("OutputFile");
  DumpOnly0TrackEvents = iConfig.getParameter<bool>("DumpOnly0TrackEvents");
  MaxTracks = iConfig.getParameter<int>("MaxTracks");
  ChiProbCut = iConfig.getParameter<double>("ChiProbCut");

  EtaMin = iConfig.getParameter<double>("EtaMin");
  EtaMax = iConfig.getParameter<double>("EtaMax");
  trackLabel = iConfig.getParameter<std::string>("TrackLabel");
}

T1Analysis_reco3_tree::~T1Analysis_reco3_tree()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
T1Analysis_reco3_tree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

  _event = iEvent.id().event();
  _tracks = 0;
  _recoHits = 0;
  _ZatRmin.clear();
  _ZImpact.clear();
  _Rmin.clear();
  _Deta.clear();
  _Dphi.clear();
  _ChiProb.clear();
  _Hits.clear();
  _Eta.clear();
  _Phi.clear();
// Trigger data

  unsigned int TriggerData_bunch_num;
  unsigned int TriggerData_run_num;
  unsigned int TriggerData_input_status_bits;

  Handle< Totem::RawEvent > input;
  iEvent.getByLabel(rawEventLabel, input);

  TriggerData_bunch_num = input->triggerData.bunch_num;
  TriggerData_run_num = input->triggerData.run_num;
  TriggerData_input_status_bits = input->triggerData.input_status_bits;
 
// T1 data

    edm::Handle<T1T2TrackCollection> T1trackCollection;
    iEvent.getByLabel(trackLabel,"T1TrackColl",T1trackCollection);

    T1T2TrackCollection::const_iterator T1TC_it;

    edm::Handle<T1RecHit2DCollection>myRecoColl;
    iEvent.getByLabel(t1RecHit2DCollectionLabel, myRecoColl);
    T1RecHit2DCollection::const_iterator T1Reco_it;

// T2 data

  Handle<T1T2TrackCollection> T2trackCollection;
  iEvent.getByLabel("T2TrackColl3","T2TrackColl",T2trackCollection);
 
  T1T2TrackCollection::const_iterator T2TC_it;


    unsigned int T2TBIT = TriggerData_input_status_bits & 64;
    unsigned int T2TBIT_HM = TriggerData_input_status_bits & 128;

    bool SelectEvent=false;


// conditions to select event to be analyzed

    bool SelectBunchNumber = false;   
    assert(TriggerData_run_num == 6917 || TriggerData_run_num == 8535 || TriggerData_run_num == 8561);
    if(TriggerData_run_num == 6917 && TriggerData_bunch_num == 0 )SelectBunchNumber=true;
    if( (TriggerData_run_num == 8535 || TriggerData_run_num == 8561) && (TriggerData_bunch_num < 800 || 
TriggerData_bunch_num > 1000)  
)SelectBunchNumber=true;

//    if(TriggerData_bunch_num == 0 && (T2TBIT || T2TBIT_HM))SelectEvent=true; // Big Bunches and T2 trigger
    if(SelectBunchNumber && (T2TBIT || T2TBIT_HM))SelectEvent=true; // Big Bunches and T2 trigger

    NEventTot++;

    if(SelectEvent){
      int RecoHitsPlus=0;
      int RecoHitsMinus=0;
      for(T1Reco_it=myRecoColl->begin();T1Reco_it!=myRecoColl->end();T1Reco_it++){
	
	int lArm=(*T1Reco_it).t1DetId().Arm();
//      _recoHits = myRecoColl->size();
	if(lArm == 0)RecoHitsPlus++;
	else 
	  RecoHitsMinus++;
      }
      _recoHits = RecoHitsPlus;
      NSelectEvent++;
      if(T1trackCollection->size()>0)NAtLeastOneTrackEvent++;
      for(T1TC_it=T1trackCollection->begin(); T1TC_it!=T1trackCollection->end(); T1TC_it++){
      //    if((*T1TC_it).Eta() > 0 && TMath::Prob(  (*T1TC_it).ChiSquared(), 2*(  (*T1TC_it).GetHitEntries() - 2))>0.01 && T1trackCollection->size()<30 ){
	if((*T1TC_it).Eta() > EtaMin &&  (*T1TC_it).Eta()<EtaMax && T1trackCollection->size()< (unsigned int)MaxTracks && TMath::Prob(  (*T1TC_it).ChiSquared(), 2*(  (*T1TC_it).GetHitEntries() - 2))>ChiProbCut){
	  
	  _tracks = T1trackCollection->size();
	  double X0 = (*T1TC_it).GetHitT1(0).GlobalPosition().x();
	  double Y0 = (*T1TC_it).GetHitT1(0).GlobalPosition().y();
	  double Z0 = (*T1TC_it).GetHitT1(0).GlobalPosition().z();
	  double R02 = X0*X0 + Y0*Y0;
	  
	  double C0 = R02/ ((*T1TC_it).GetTx()*X0 + (*T1TC_it).GetTy()*Y0);
	  double ZImpact = Z0 - C0;
	  
	  double Dphi=0;
	  double Deta=0;
    

	  hZatRmin->Fill((*T1TC_it).Z_at_Rmin());
	  hRmin->Fill((*T1TC_it).Rmin());
	  hZImpact->Fill(ZImpact);
	  hZetas->Fill((*T1TC_it).Z_at_Rmin(),ZImpact);
	  hZR->Fill((*T1TC_it).Z_at_Rmin(),(*T1TC_it).Rmin());
	  hChiOverN->Fill( (*T1TC_it).ChiSquaredOverN() );
	  hChiProb->Fill(    TMath::Prob(  (*T1TC_it).ChiSquared(), 2*(  (*T1TC_it).GetHitEntries() - 2)) );
	      for(unsigned int iii = 0; iii<(*T1TC_it).GetHitEntries(); iii++){
 
		Dphi +=  (*T1TC_it).GetHitT1(iii).phi()/(double)(*T1TC_it).GetHitEntries();  
		Deta +=  (*T1TC_it).GetHitT1(iii).eta()/(double)(*T1TC_it).GetHitEntries();  
	      }
	      Dphi -= (*T1TC_it).Phi();
	      Deta -= (*T1TC_it).Eta();

    _ZatRmin.push_back((*T1TC_it).Z_at_Rmin());
    _ZImpact.push_back(ZImpact);
    _ChiProb.push_back(TMath::Prob(  (*T1TC_it).ChiSquared(), 2*(  (*T1TC_it).GetHitEntries() - 2)));
    _Rmin.push_back((*T1TC_it).Rmin());
    _Deta.push_back(Deta);
    _Dphi.push_back(Dphi);
    _Hits.push_back((*T1TC_it).GetHitEntries() );
    _Eta.push_back((*T1TC_it).Eta() );
    _Phi.push_back((*T1TC_it).Phi() );
	}
      }

  _tree->Fill();

    }
    


}






// ------------ method called once each job just before starting event loop  ------------
void 
T1Analysis_reco3_tree::beginJob()
{
 _tree.reset(new TTree("tree","T1tracks"));

  _tree->Branch("event",&_event);
  _tree->Branch("tracks",&_tracks);
  _tree->Branch("recoHits",&_recoHits);
  _tree->Branch("ZatRmin",&_ZatRmin);
  _tree->Branch("ZImpact",&_ZImpact);
  _tree->Branch("Rmin",&_Rmin);
  _tree->Branch("ChiProb",&_ChiProb);
  _tree->Branch("Deta",&_Deta);
  _tree->Branch("Dphi",&_Dphi);
  _tree->Branch("Hits",&_Hits);
  _tree->Branch("Eta",&_Eta);
  _tree->Branch("Phi",&_Phi);


   NSelectEvent=0;
   NAtLeastOneTrackEvent=0;
   NEventTot=0;


  hZatRmin = std::auto_ptr<TH1D>(new TH1D("ZatRmin","ZatRmin",1000,-50000,50000));
  hRmin = std::auto_ptr<TH1D>(new TH1D("Rmin","Rmin",1000,0,2000));
  hZImpact = std::auto_ptr<TH1D>(new TH1D("ZImpact","ZImpact",1000,-50000,50000));
  hZetas   = std::auto_ptr<TH2D>(new TH2D("Zetas  ","Zetas  ",1000,-50000,50000,1000,-50000,50000));
  hZR   = std::auto_ptr<TH2D>(new TH2D("ZR  ","ZR  ",200,-50000,50000,200,0,1000));
  hChiOverN = std::auto_ptr<TH1D>(new TH1D("ChiOverN","ChiOverN",1000,0,100));
  hChiProb = std::auto_ptr<TH1D>(new TH1D("ChiProb","ChiProb",1000,0,1));


  hRmin->SetDirectory(0);
  hZR->SetDirectory(0);
  hZatRmin->SetDirectory(0);
  hZImpact->SetDirectory(0);
  hZetas->SetDirectory(0);
  hChiOverN->SetDirectory(0);
  hChiProb->SetDirectory(0);

 _tree->SetDirectory(0);

}
// ------------ method called once each job just after ending the event loop  ------------
void 
T1Analysis_reco3_tree::endJob() 
{

  theFile = TFile::Open(outputFileName.c_str(), "recreate");
  if(!theFile || !theFile->IsWritable())
    {
      std::cout<<"Output file not opened correctly!!"<<std::endl;
    }
  hZatRmin->Write();
  hZImpact->Write();
  hZetas->Write();
  hRmin->Write();
  hZR->Write();
  hChiOverN->Write();
  hChiProb->Write();



  cout << "Events: Total="<<NEventTot << " Selected="<<NSelectEvent <<" >1Track="<<NAtLeastOneTrackEvent<<endl;





 _tree->Write();





  theFile->Close();





 }

// ------------ method called when starting to processes a run  ------------
void 
T1Analysis_reco3_tree::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
T1Analysis_reco3_tree::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
T1Analysis_reco3_tree::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
T1Analysis_reco3_tree::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
T1Analysis_reco3_tree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}




//define this as a plug-in
DEFINE_FWK_MODULE(T1Analysis_reco3_tree);

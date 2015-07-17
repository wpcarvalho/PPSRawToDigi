#include "IOMC/DPEProtons/interface/DPEProtons.h"
#include <ostream>
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <iostream>
#include <TMath.h>



DPEProtons::DPEProtons(const edm::ParameterSet& pset)
{
  file_names = pset.getUntrackedParameter< std::vector<std::string>  >("FileNames");
  verbosity_ = pset.getUntrackedParameter<int>("Verbosity");
  
  dpe_file_access = new h101(file_names);
  max_entries = dpe_file_access->GetEntries();
  
  produces<edm::HepMCProduct>();
  
  std::cout << "Internal DPEProtons is initialzed" << std::endl ;
  
  edm::Service<edm::RandomNumberGenerator> rng;
  long seed = (long)(rng->mySeed());
  
  fRandomEngine = new CLHEP::HepJamesRandom(seed);
  fRandomGenerator = new CLHEP::RandFlat(fRandomEngine);
  
  //  if(verbosity_)
}

DPEProtons::~DPEProtons()
{
  if(fRandomGenerator)
    delete fRandomGenerator;
  if(dpe_file_access)
    delete dpe_file_access;
}

void DPEProtons::beginRun(Run& r, EventSetup const& es){
  es.getData(fPDGTable);
  PartID_ = 2212;
  const HepPDT::ParticleData* PData = fPDGTable->particle(HepPDT::ParticleID(PartID_));
  m0_ = PData->mass().value();
  
  return ;
}

void DPEProtons::produce(edm::Event& e, const edm::EventSetup& es)
{
  if(verbosity_)
    std::cout<<" DPEProtons : Begin New Event Generation"<<std::endl;
  
  fEvt_ = new HepMC::GenEvent() ;
  HepMC::GenVertex* Vtx = new HepMC::GenVertex(HepMC::FourVector(0.,0.,0.));
  
  long ev_number = e.id().event();
  if(ev_number<0 || ev_number>=max_entries)
  {
    long ev_number = fRandomGenerator->fireInt(0, max_entries);
    std::cout<<"DPEProtons : Data input set size exceeded. Too many events generated. "
    << ev_number << " out of " << max_entries <<std::endl;
  }
  MADProtonPair pair = dpe_file_access->GetProtonPair(ev_number);
  
  int barcode = 1;
  
  double px1, py1, pz1, E1;
  px1 = pair.l.px;
  py1 = pair.l.py;
  pz1 = pair.l.pz;
  E1 = TMath::Sqrt(px1*px1 + py1*py1 + pz1*pz1 + m0_*m0_);
  
  double px2, py2, pz2, E2;
  px2 = pair.r.px;
  py2 = pair.r.py;
  pz2 = pair.r.pz;
  E2 = TMath::Sqrt(px2*px2 + py2*py2 + pz2*pz2 + m0_*m0_);
  
  HepMC::ThreeVector p1(px1,py1,pz1);
  HepMC::GenParticle* Part1 =
    new HepMC::GenParticle(HepMC::FourVector(px1,py1,pz1,E1),2212,1);
  Part1->suggest_barcode(barcode);
  barcode++ ;
  Vtx->add_particle_out(Part1);
  
  HepMC::ThreeVector p2(px2,py2,pz2);
  HepMC::GenParticle* Part2 =
    new HepMC::GenParticle(HepMC::FourVector(px2,py2,pz2,E2),2212,1);
  Part2->suggest_barcode(barcode);
  barcode++ ;
  Vtx->add_particle_out(Part2);
  
  fEvt_->add_vertex(Vtx);
  fEvt_->set_event_number(e.id().event());
  fEvt_->set_signal_process_id(20);
  
  if(verbosity_)
  {
    fEvt_->print();
  }
  
  std::auto_ptr<edm::HepMCProduct> BProduct(new edm::HepMCProduct()) ;
  BProduct->addHepMCData(fEvt_);
  e.put(BProduct);
  
  if(verbosity_)
  {
    std::cout << "DPEProtons : Event Generation Done, "<<std::endl;
  }
  
}

DEFINE_FWK_MODULE(DPEProtons);

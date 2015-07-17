#include "RecoTotemRP/CentralMCJetReconstruction/interface/CentralMCJetReconstruction.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include <string>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "RecoTotemRP/CentralMCJetReconstruction/interface/BuildClusters.h"
#include <iostream>
#include "RecoTotemRP/RPRecoDataFormats/interface/CentralMassInfo.h"
#include "RecoTotemRP/CentralMCJetReconstruction/interface/PseudoJetVectorWrapper.h"


using namespace std;
//using namespace fastjet;

CentralMCJetReconstruction::CentralMCJetReconstruction(const edm::ParameterSet& conf)
{
  hepmc_instance_ = conf.getParameter<string>("hepmcInstance");
  hepmc_label_ = conf.getParameter<string>("hepmcLabel");
  jet_output_label_ = conf.getParameter<string>("jetOutputLabel");
  verbosity_ = conf.getParameter<int>("verbosity");
  R_ = conf.getParameter<double>("R");
  min_eta_ = conf.getParameter<double>("minEta");
  max_eta_ = conf.getParameter<double>("maxEta");
  min_pt_ = conf.getParameter<double>("minPt");
  beamEnergy_ = conf.getParameter<double>("beamEnergy");
  findJets_ = conf.getParameter<bool>("findJets");
  
  string jet_alg_name = conf.getParameter<string>("jetAlgName");

  if(jet_alg_name == "kt_algorithm")
    jet_alg_ = fastjet::kt_algorithm;
  else if(jet_alg_name == "cambridge_algorithm")
    jet_alg_ = fastjet::cambridge_algorithm;
  else if(jet_alg_name == "antikt_algorithm")
    jet_alg_ = fastjet::antikt_algorithm;
  else if(jet_alg_name == "genkt_algorithm")
    jet_alg_ = fastjet::genkt_algorithm;
  else if(jet_alg_name == "ee_kt_algorithm")
    jet_alg_ = fastjet::ee_kt_algorithm;
  else if(jet_alg_name == "ee_genkt_algorithm")
    jet_alg_ = fastjet::ee_genkt_algorithm;
  else
  {
    cout<<"CentralMCJetReconstruction :: wrongly defined jet algorithm, exit(0)"<<endl;
  }
  
  produces< PseudoJetVectorWrapper >(jet_output_label_);
  produces<CentralMassInfo>(jet_output_label_);
}

CentralMCJetReconstruction::~CentralMCJetReconstruction()
{
}

void CentralMCJetReconstruction::beginJob()
{
}

void CentralMCJetReconstruction::produce(edm::Event& e, const edm::EventSetup& c)
{
  edm::Handle<edm::HepMCProduct> mcProd;
  e.getByLabel(hepmc_label_, mcProd);  
  HepMC::GenEvent *evt = (HepMC::GenEvent *) mcProd->GetEvent();
  
  BuildClusters cluster_builder(beamEnergy_);
  cluster_builder.AddParticles(evt);
  vector<fastjet::PseudoJet> rec_jets;
  if(findJets_)
  {
    rec_jets = cluster_builder.GetInclusiveJets(jet_alg_, R_, min_pt_, min_eta_, max_eta_);
  }

  PseudoJetVectorWrapper * jetWrapper = new PseudoJetVectorWrapper();
  jetWrapper->pseudoJetsVector = vector<fastjet::PseudoJet>(rec_jets);
  auto_ptr< PseudoJetVectorWrapper > output(jetWrapper);
  e.put(output, jet_output_label_);
  
  auto_ptr< CentralMassInfo > centr_mass_info_output(new CentralMassInfo(cluster_builder.GetCentralMassInfo()));
  e.put(centr_mass_info_output, jet_output_label_);
}





DEFINE_FWK_MODULE(CentralMCJetReconstruction);

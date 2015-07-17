#ifndef RecoTotemRPCentralMCJetReconstruction_h__
#define RecoTotemRPCentralMCJetReconstruction_h__


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "fastjet/JetDefinition.hh"
#include <string>

using namespace std;
//using namespace fastjet;

class CentralMCJetReconstruction : public edm::EDProducer
{
  public:
    explicit CentralMCJetReconstruction(const edm::ParameterSet& conf);  
    virtual ~CentralMCJetReconstruction();
    virtual void beginJob();  
    virtual void produce(edm::Event& e, const edm::EventSetup& c);
  
  private:
//    edm::ParameterSet conf_;
    string hepmc_instance_;
    string hepmc_label_;
    string jet_output_label_;
    int verbosity_;
    double R_;
    double min_eta_, max_eta_;
    fastjet::JetAlgorithm jet_alg_;
    double min_pt_;
    double beamEnergy_;
    bool findJets_;
};


#endif

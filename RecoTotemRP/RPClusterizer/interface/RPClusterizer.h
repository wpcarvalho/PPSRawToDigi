/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*   Hubert Niewiadomski, CERN
*
****************************************************************************/

#ifndef RPClusterizer_h
#define RPClusterizer_h

#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/TotemRPDigi/interface/TotemRPDigi.h"

#include "RecoTotemRP/RPClusterizer/interface/RPClusterizerAlgorithm.h"
 
#include <iostream>
#include <memory>
#include <string>
 
/**
 * Merges neighbouring active TOTEM RP strips into clusters.
 **/
class RPClusterizer : public edm::one::EDProducer<>
{
  public:
  
    explicit RPClusterizer(const edm::ParameterSet& conf);
  
    virtual ~RPClusterizer();
  
    virtual void beginJob();
  
    virtual void produce(edm::Event& e, const edm::EventSetup& c);
  
  private:
    edm::ParameterSet conf_;
    int verbosity_;
    edm::InputTag digiInputTag_;
    edm::EDGetTokenT<edm::DetSetVector<TotemRPDigi> >digiInputTagToken_;

    RPClusterizerAlgorithm RPClusterizerAlgorithm_;

    void run(const edm::DetSetVector<TotemRPDigi> &input, edm::DetSetVector<TotemRPCluster> &output);
};
  
#endif

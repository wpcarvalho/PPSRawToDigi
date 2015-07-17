//--------------------------------------------------------------------------------------------------
// $Id: ProducerEvtSelData.h,v 1.5 2010/01/07 17:07:55 loizides Exp $
//
// ProducerEvtSelData
//
// Produce event selection data.
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITEDM_PRODUCERS_PRODUCEREVTSELDATA_H
#define MITEDM_PRODUCERS_PRODUCEREVTSELDATA_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TrackerGeometry;

namespace mitedm
{
  class ProducerEvtSelData : public edm::EDProducer {
    public:
      explicit ProducerEvtSelData(const edm::ParameterSet &cfg);
      ~ProducerEvtSelData();
    
    private:
      struct VertexHit {
        float z;
        float r;
        float w;
      };

      void produce(edm::Event &evt, const edm::EventSetup &setup);
      int  getContainedHits(const std::vector<VertexHit> &hits, double z0, double &chi);

      std::string srcHF_;     //hf rec hits
      std::string srcHBHE_;   //hbhe rec hits
      std::string srcCastor_; //castor rec hits
      std::string srcZDC_;    //zdc rec hits
      std::string srcPixels_; //pixel rec hits
      std::string srcVertex_; //vertex (if not set will use pixel counting vertex)
      std::string srcTowers_; //calo towers
      double      hfEthresh_; //hf hit energy threshold
      double      hfETowerh_; //hf calo tower energy threshold
      std::string srcTrk_;    //track collection
  };
}
#endif

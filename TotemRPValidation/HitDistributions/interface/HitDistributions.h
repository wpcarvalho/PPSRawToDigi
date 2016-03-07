#ifndef _TotemRPValidationHitDistributionsHitDistributions_H_
#define _TotemRPValidationHitDistributionsHitDistributions_H_

// -*- C++ -*-
//
// Package:    HitDistributions
// Class:      HitDistributions
//
// Author: Leszek Grzanka

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "TotemRPValidation/HitDistributions/interface/HitDistributionsLibrary.h"

#include <TGraph.h>
#include <map>
#include <boost/shared_ptr.hpp>

namespace edm {
  class ParameterSet;
  class EventSetup;
  class Event;
}


class HitDistributions : public edm::EDAnalyzer
{
public:
  explicit HitDistributions(const edm::ParameterSet&);
  ~HitDistributions();

private:
  unsigned int verbosity;
  unsigned int validTracksOnly;
  std::auto_ptr<HitDistributionsLibrary> library;

  std::string outputFile;
  std::map< unsigned int, TGraph *> hitDists;

  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
};


#endif

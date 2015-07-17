// -*- C++ -*-
//
// Class:      RPCoincidenceAnalyzer
//
// Author: Leszek Grzanka

#ifndef _L1TriggerTotemCoincidenceChipRPCoincidenceAnalyzer_H_
#define _L1TriggerTotemCoincidenceChipRPCoincidenceAnalyzer_H_

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDetTrigger.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/TotemL1Trigger/interface/RPCCBits.h"

#include "TotemRawDataLibrary/DataFormats/interface/RawEvent.h"
#include "TotemRawDataLibrary/DataFormats/interface/VFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/VFATFrameCollection.h"

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"

#include <map>
#include <vector>
#include <cstdlib>
#include <boost/shared_ptr.hpp>

namespace edm {
  class ParameterSet;
  class EventSetup;
  class Event;
}


class RPCoincidenceAnalyzer : public edm::EDAnalyzer
{
public:
  explicit RPCoincidenceAnalyzer(const edm::ParameterSet&);
  ~RPCoincidenceAnalyzer();

private:
  unsigned int verbose_;
  std::string outputFile;

  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
  std::auto_ptr<TH1I> tracksSelectionCategories;

  std::auto_ptr<TH1D> trigSectorsRawEven;
  std::auto_ptr<TH1D> trigSectorsRawOdd;
  std::auto_ptr<TH1D> trigSectorsSimuEven;
  std::auto_ptr<TH1D> trigSectorsSimuOdd;

  std::auto_ptr<TH1I> differenceSimuVsRawEven;
  std::auto_ptr<TH1I> differenceSimuVsRawOdd;
  std::auto_ptr<TH1I> differenceSimuVsRawTotal;

  std::auto_ptr<TH1I> differenceCategories;

  std::auto_ptr<TH1I> goodTracksCategories;
  std::auto_ptr<TH1I> activeChipSim;
  std::auto_ptr<TH1I> activeChipRaw;

  std::auto_ptr<TH1I> goodEventsCategories;

  std::auto_ptr<TH1I> totalEfficiencyRaw;
  std::auto_ptr<TH1I> totalEfficiencySimu;

  std::vector<unsigned int> getEvenSectorsVFAT( std::vector<RPStripDigi>  );
  std::vector<unsigned int> getOddSectorsVFAT( std::vector<RPStripDigi>  );
  
  std::string  modulLabelRaw;
  std::string  productLabelRaw;
  std::string  modulLabelSimu;
  std::string  productLabelSimu; 
  
  edm::InputTag trackCandidateCollectionLabel;
  edm::InputTag detTriggerLabel;

};


#endif

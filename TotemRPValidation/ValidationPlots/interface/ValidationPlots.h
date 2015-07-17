/*
 * ValidationPlots.h
 *
 *  Created on: Oct 29, 2008
 *      Author: Leszek Grzanka
 */

#ifndef _TotemRPValidationValidationPlotsValidationPlots_H_
#define _TotemRPValidationValidationPlotsValidationPlots_H_


#include "FWCore/Framework/interface/EDAnalyzer.h"
#include <TGraph.h>
#include <TCanvas.h>
#include <map>

#include "TotemRPValidation/ValidationPlots/interface/InelasticPlots.h"
#include "TotemRPValidation/ElasticReconstruction/interface/ElasticRecoValLibrary.h"
#include "TotemRPValidation/HitDistributions/interface/HitDistributionsLibrary.h"


namespace edm {
  class ParameterSet;
  class EventSetup;
  class Event;
}


class ValidationPlots : public edm::EDAnalyzer
{
 public:
  explicit ValidationPlots(const edm::ParameterSet&);
  ~ValidationPlots();

 private:
  unsigned char verbosity;
  std::string outputFile;
  ElasticRecoValLibrary * elasticHists;
  HitDistributionsLibrary * hitDistributions;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
};

#endif /* _TotemRPValidationValidationPlotsValidationPlots_H_ */

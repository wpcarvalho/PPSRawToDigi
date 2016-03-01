/****************************************************************************
*
* This is a part of the TOTEM offline software.
* Authors: 
*	Jan Kašpar (jan.kaspar@gmail.com)
*    
* $Id: ElasticRecoVal.h 9977 2015-01-12 14:00:26Z tsodzawi $
* $Revision: 9977 $
* $Date: 2015-01-12 15:00:26 +0100 (Mon, 12 Jan 2015) $
*
****************************************************************************/


#ifndef _TotemRPValidationElasticReconstructionElasticRecoVal_H_
#define _TotemRPValidationElasticReconstructionElasticRecoVal_H_

#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "TotemRPValidation/ElasticReconstruction/interface/ElasticRecoValLibrary.h"
#include "TotemRPValidation/ElasticReconstruction/interface/ElasticAcceptanceLibrary.h"

#include <TH1D.h>

#include <boost/shared_ptr.hpp>

namespace edm {
  class ParameterSet;
  class EventSetup;
  class Event;
}

/**
 *\brief Plugin for elastic reconstruction validation.
 **/
class ElasticRecoVal : public edm::EDAnalyzer
{
 public:
  explicit ElasticRecoVal(const edm::ParameterSet&);
  ~ElasticRecoVal();

 private:

  std::auto_ptr<ElasticRecoValLibrary> libVal;
  std::auto_ptr<ElasticAcceptanceLibrary> libAcc;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
};

#endif 


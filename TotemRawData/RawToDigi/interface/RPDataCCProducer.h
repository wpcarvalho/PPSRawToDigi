/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*   Jan Kašpar (jan.kaspar@gmail.com)
*   Leszek Grzanka 
*
****************************************************************************/

#ifndef DataCCProducer_h
#define DataCCProducer_h

#include "FWCore/Framework/interface/EDProducer.h"

namespace edm {
  class ParameterSet;
  class EventSetup;
  class Event;
  class EventID;
}

namespace Totem {
  class FramePosition;
}


/**
 * \brief Converts raw event data to RP CC information.
**/
class RPDataCCProducer : public edm::EDProducer
{
  public:

    RPDataCCProducer(const edm::ParameterSet& conf);
    virtual ~RPDataCCProducer();

    virtual void beginJob(const edm::EventSetup&);
    virtual void endJob();
    virtual void produce(edm::Event&, const edm::EventSetup&);

  private:
    unsigned char verbosity;
	std::string productLabelRaw;
};

#endif


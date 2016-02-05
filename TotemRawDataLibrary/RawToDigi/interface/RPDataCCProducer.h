/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*   Jan Kašpar (jan.kaspar@gmail.com)
*   Leszek Grzanka 
*
* $$RCSfile: RPDataCCProducer.h,v $: $
* $Revision: 9953 $
* $Date: 2014-12-16 12:03:45 +0100 (Tue, 16 Dec 2014) $
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


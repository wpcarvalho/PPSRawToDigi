/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#ifndef T1XMLDataDigiProducer_h
#define T1XMLDataDigiProducer_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "TotemRawData/RawToDigi/interface/T1FisChannel.h"
#include "DataFormats/T1DigiVfat/interface/T1DigiVfat.h"
#include "DataFormats/T1DigiVfat/interface/T1DigiVfatCollection.h"
#include "DataFormats/T1DigiWire/interface/T1DigiWire.h"
#include "DataFormats/T1DigiWire/interface/T1DigiWireCollection.h"

#include "TotemRawDataLibrary/DataFormats/interface/VFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/RawEvent.h"

#include "TotemRawDataLibrary/DataFormats/interface/CommonDef.h"
#include "TotemCondFormats/DAQInformation/interface/DAQInformationT1.h"
#include "TotemCondFormats/DataRecord/interface/TotemDAQRecord.h"
#include "TotemRawData/RawToDigi/interface/T1Cfecs.h"
#include <math.h>

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
 * \brief Converts raw event data to RP digi information.
**/
class T1XMLDataDigiProducer : public edm::EDProducer
{
  public:

    T1XMLDataDigiProducer(const edm::ParameterSet& conf);
    virtual  ~T1XMLDataDigiProducer();

    virtual void beginJob(const edm::EventSetup&);
    virtual void endJob();
    virtual void produce(edm::Event&, const edm::EventSetup&);
//      int FisKey(int,int,int,int,int,int);
    //   virtual void generateMap( const edm::EventSetup&);
// virtual void generateMap ( map<FramePosition,unsigned int> );

  private:
    unsigned char verbosity;
    std::string _NOISY_MAP;
    unsigned int _GET_NOISY_MAP;
	  bool starting;
    bool eventError, positionError;   ///< whether error found in event/position
 map<int,T1FisChannel> T1Map;
    ///< prints error header
    void PrintErrorHeader(const edm::EventID&, const Totem::FramePosition&, unsigned short, unsigned short);
};





#endif


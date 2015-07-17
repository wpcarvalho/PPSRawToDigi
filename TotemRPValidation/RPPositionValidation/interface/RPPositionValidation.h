#ifndef TotemRPValidation_RPAngleValidation_RPAngleValidation_h
#define TotemRPValidation_RPAngleValidation_RPAngleValidation_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"

namespace edm {
  class ParameterSet;
  class EventSetup;
  class Event;
}

class RPPositionValidation : public edm::EDAnalyzer
{
  public:
    explicit RPPositionValidation(const edm::ParameterSet&);
    ~RPPositionValidation();
    
  private:
    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

	#define NSTATIONS 4
	#define NPOTS 6
    static const int rpids[NSTATIONS][NPOTS];
    edm::ParameterSet conf_;
    unsigned int verbosity_;
};


#endif


#include <iostream>
#include <iomanip>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"
#include "Geometry/TotemRecords/interface/RealGeometryRecord.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "TotemRPValidation/RPPositionValidation/interface/RPPositionValidation.h"

using namespace std;

const int RPPositionValidation::rpids[NSTATIONS][NPOTS] = {
       {   0,   1,   2,   3,   4,   5 },
       {  20,  21,  22,  23,  24,  25 },
       { 100, 101, 102, 103, 104, 105 },
       { 120, 121, 122, 123, 124, 125 }
       };


RPPositionValidation::RPPositionValidation(const edm::ParameterSet& conf)
 : conf_(conf)
{
    verbosity_ = conf.getUntrackedParameter<unsigned int>("Verbosity");
}

RPPositionValidation::~RPPositionValidation()
{
}

void RPPositionValidation::beginRun(edm::Run const& r, edm::EventSetup const& es)
{
  edm::ESHandle<TotemRPGeometry> geometry;
  es.get<RealGeometryRecord>().get(geometry);
  
  if (verbosity_ > 0) {
  	cout << "Thin foil and detector positions and surface normals" << endl;

  for (int s=0; s<NSTATIONS; s++)
    for (int p=0; p<NPOTS; p++) {
      int RPId = rpids[s][p];

      CLHEP::Hep3Vector fp = geometry->GetRPThinFoilPosition(RPId);
      CLHEP::Hep3Vector fv = geometry->GetRPThinFoilNormalVector(RPId);

      cout << setfill('0') << endl << "RP " << setw(3) << RPId << ": thin foil at " << fp << ", surface normal " << fv  << endl;

      for (int d = 0; d < 10; d++) {
		unsigned int detId = RPId * 10 + d;
        unsigned int rawDetId = TotRPDetId::DecToRawId(detId);

        CLHEP::Hep3Vector ep = geometry->GetDetEdgePosition(rawDetId);
        CLHEP::Hep3Vector ev = geometry->GetDetEdgeNormalVector(rawDetId);
        CLHEP::Hep3Vector center = geometry->LocalToGlobal(rawDetId, CLHEP::Hep3Vector(0, 0, 0));

        cout << setw(4) << detId << ": det_e=" << ep << ", det_e_n=" << ev    << endl;
        cout << "      foil - det_e=" << fp-ep << ", foil_n - det_e_n=" << fv-ev << endl;
        cout << "      det_c=" << center << ", det_e - det_c " << (ep - center) << endl;
        cout << "      foil - det_c=" << (fp - center) << endl;
      }
    }
  }
}

void RPPositionValidation::analyze(const edm::Event& e, const edm::EventSetup& es)
{
}

void RPPositionValidation::endJob()
{
}

DEFINE_FWK_MODULE(RPPositionValidation);

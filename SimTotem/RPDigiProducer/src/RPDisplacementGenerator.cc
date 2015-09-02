#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TotemRecords/interface/MisalignedGeometryRecord.h"
#include "Geometry/TotemRecords/interface/RealGeometryRecord.h"
#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"
#include "TotemAlignment/RPDataFormats/interface/RPAlignmentCorrections.h"

#include "SimTotem/RPDigiProducer/interface/RPDisplacementGenerator.h"

using namespace std;
using namespace edm;

//#define DEBUG 1


RPDisplacementGenerator::RPDisplacementGenerator(const edm::ParameterSet &ps, RPDetId _detId, const edm::EventSetup &iSetup) : detId(_detId)
{
  isOn = ps.getParameter<bool>("RPDisplacementOn");

  // read the alignment correction
  ESHandle<RPAlignmentCorrections> alignments;
  try {
    iSetup.get<MisalignedGeometryRecord>().get(alignments);
  }
  catch (...) {
  }

  unsigned int decId = TotRPDetId::RawToDecId(detId);

  DDTranslation S_m;        // zero translation by default
  DDRotationMatrix R_m;     // identity rotation by default

  if (alignments.isValid()) {
    const RPAlignmentCorrection& ac = alignments->GetFullSensorCorrection(decId);
    S_m = ac.Translation(); 
    R_m = ac.RotationMatrix(); 
  } else
    isOn = false;

  // transform shift and rotation to the local coordinate frame
  ESHandle<TotemRPGeometry> geom;
  iSetup.get<RealGeometryRecord>().get(geom);
  DetGeomDesc *g = geom->GetDetector(detId);

  //const DDTranslation& S_l = g->translation();
  const DDRotationMatrix& R_l = g->rotation();

  rotation = R_l.Inverse() * R_m.Inverse() * R_l;
  shift = R_l.Inverse() * R_m.Inverse() * S_m;

#ifdef DEBUG
  cout << ">> RPDisplacementGenerator::RPDisplacementGenerator, det id = " << decId << ", isOn = " << isOn << endl;
  if (isOn) {
    cout << " shift = " << shift << endl;
    cout << " rotation = " << rotation << endl;
  }
#endif
}

//----------------------------------------------------------------------------------------------------

Local3DPoint RPDisplacementGenerator::DisplacePoint(const Local3DPoint &p)
{
  /// input is in mm, shifts are in mm too
  
  //printf("DisplacePoint:\n");
  //printf("\t%E\t%E\t%E\n", p.x(), p.y(), p.z());

  DDTranslation v(p.x(), p.y(), p.z());
  v = rotation * v - shift;

  //printf("\t%E\t%E\t%E\n", v.x(), v.y(), v.z());

  return Local3DPoint(v.x(), v.y(), v.z());
}

//----------------------------------------------------------------------------------------------------

PSimHit RPDisplacementGenerator::Displace(const PSimHit &input)
{
  if (!isOn)
    return input;

  const Local3DPoint &ep = input.entryPoint(), &xp = input.exitPoint();
  const Local3DPoint &dep = DisplacePoint(ep), &dxp = DisplacePoint(xp);

#ifdef DEBUG
  printf(">> RPDisplacementGenerator::Displace\n");
  cout << " entry point: " << ep << " -> " << dep << endl;
  cout << " exit point : " << xp << " -> " << dxp << endl;
#endif


  return PSimHit(dep, dxp, input.pabs(), input.tof(), input.energyLoss(), input.particleType(),
    input.detUnitId(), input.trackId(), input.thetaAtEntry(), input.phiAtEntry(), input.processType());
}


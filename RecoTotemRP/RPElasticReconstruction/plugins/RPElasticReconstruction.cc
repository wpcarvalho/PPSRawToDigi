/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/ProtonTransportFunctions/interface/ProtonTransportFunctions.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include "TotemCondFormats/DataRecord/interface/ProtonTransportRcd.h"
#include "SimG4CMS/TotemRPProtTranspPar/interface/LHCOpticsApproximator.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecoElasticEvent.h"

#include "TGraph.h"

#include <iostream>
#include <map>
#include <vector>



///\defgroup ElasticReconstruction
/**
 * \ingroup ElasticReconstruction
 * \brief Class for reconstruction of elastic events.
 *
 * \b Note: Current implementation treats hits with fixed uncertainty 66/sqrt(12) um = 19.05 um, see recoHit_type
 *
 * \b Note 2: In this module, the angle theta is considered as a PARAMETER of an elastic PROCESS (i.e. theta = sqrt(|t|) / p). 
 * Its components theta_x and theta_y are defined as 
 * (1) theta_x = theta * cos(phi), theta_y = theta * sin(phi),
 * where phi is the azimuthal angle (the second parameter describing the elastic process). This definition differs from another
 * common one, where theta_x and theta_y are related to the direction of a PARTICLE. They're defined
 * (2) theta_x = p_x / |p|, theta_y = p_y / |p|,
 * where p denotes the momentum of the particle.
 * Both definitions concide for the elastic proton with positive p_z, but they differ in sign for the elastic proton with negative p_z.
 * This effects the way how the proton transport is handled. The common way is give the following equation
 * (3) x_hit = Lx * theta_x + vx * x_star
 * This holds also for the proton with positive p_z when the definition (1) is used. But for the proton with negative p_z, one has 
 * to add a minus sign:
 * (4) x_hit = Lx * (-theta_x) + vx * x_star
 * Equations (3) and (4) can be written in a common form
 * (5) x_hit = Lx' * theta_x + vx * x_star,
 * where the Lx' is a kind of "modified" effective length, which becomes negative for RP at the left arm. This modified effective lengths 
 * are used throughout this module. See optics_type(LHCOpticsApproximator *, bool invert = false) method for more information.
 **/
class RPElasticReconstruction : public edm::EDProducer
{
  public:
    /// structure holding information on optical function for a RP
    struct optics_type {
      double Cx, Cy;  ///< in m
      double Lx, Ly;  ///< in m
      double vx, vy;  ///< dimensionless
      optics_type(double _Lx = 0., double _Ly = 0., double _vx = 0., double _vy = 0.) : Lx(_Lx), Ly(_Ly), vx(_vx), vy(_vy) {}
      optics_type(LHCOpticsApproximator *, const BeamOpticsParams &optPar, bool invert = false);
    };

    /// utility structure for a hit (= pair of hit position and optical function), used for fitting
    struct recoHit_type {
      double x, si_x, y, si_y;  ///< hit information; in m
      optics_type *of;      	///< information on optical function
      
      /// constructor; expects x, si_x, y, si_y in m
      recoHit_type(double _x = 0., double _si_x = 0., double _y = 0., double _si_y = 0., optics_type *_of = NULL)
        : x(_x), si_x(_si_x), y(_y), si_y(_si_y), of(_of) {}

      /// constructor; converts mm to m, removes the effect by crossing angle (in 1st order) 
      recoHit_type(const RPFittedTrack &ft, optics_type *_of);
    };

    explicit RPElasticReconstruction(const edm::ParameterSet& conf);
    virtual ~RPElasticReconstruction();

    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void produce(edm::Event& e, const edm::EventSetup& c);
  
  protected:
    /// verbosity level
    unsigned int verbosity;

    /// what to do with IP information
    /// 0 ... do not use it at all, 1 ... add it only to left and right fits, 2 ... add to all fits
    unsigned int addIPHit;

    /// roadsize limits
    double roadsize_x, roadsize_y;

    /// beam spread at IP
    double beamsize_x, beamsize_y;

    /// tolerance for angular difference for left, right and global fit
    double angleTolerance_x, angleTolerance_y;

    /// tolerance for vertex position difference for left and right fit
    double vertexTolerance_x, vertexTolerance_y;

    // parameters of edge effect correction
    bool eec_enabled;
    double eec_si_bd, eec_th_L;

    /// graph of th_y_reco vs. th_y_orig
    TGraph *eec_graph;  

    /// loads file with optics and builds `optics' variable
    void InitializeOpticalFunctions(const ProtonTransportFunctions &optFun, const BeamOpticsParams &optPar);

    /// performs road search
    unsigned char MakeRoads(const RPFittedTrackCollection *, RPRecoElasticEvent *, const BeamOpticsParams&);

    /// performs fits (left, right and global)
    unsigned char MakeFits(const RPFittedTrackCollection *, RPRecoElasticEvent *, const BeamOpticsParams&);

    /// fits collection of recoHit_type
    RPRecoElasticEvent::fit_type DoFit(std::vector<recoHit_type> &);

    /// merges fit information (left, right and global)
    unsigned char MakeResult(const RPFittedTrackCollection *, RPRecoElasticEvent *, const BeamOpticsParams&);

    // builds the correction graph
    void InitializeEdgeEffectCorrection();

    // applies edge effet correction by using the correction graph
    void ApplyEdgeEffectCorrection(RPRecoElasticEvent::fit_type&);

    /// collection of optics information
    std::map<unsigned int, optics_type> optics;
    optics_type opticsIP;

    /// safe method to get optical functions for a RP
    optics_type& GetOptics(unsigned int RPID)
    {
      std::map<unsigned int, optics_type>::iterator it = optics.find(RPID);
      if (it != optics.end())
        return it->second;
      else
        throw cms::Exception("RPElasticReconstruction") << "RP ID=" << RPID << " has no optical functions assigned. Check your CFG file." << std::endl;
    }

    /// Is RP with ID `RPID' at left arm?
    bool IsRPLeft(unsigned int RPID)
    {
      return ((RPID / 100) == 0);
    }
private:
    edm::InputTag fittedTrackCollectionInputLabel;
};

//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

RPElasticReconstruction::recoHit_type::recoHit_type(const RPFittedTrack &ft, optics_type *_of) : of(_of)
{
  x = ft.X0()*1E-3 - of->Cx; 
  y = ft.Y0()*1E-3 - of->Cy; 

  /// \b error is \b fixed to 66/sqrt(12) um
  //si_x = ft.X0Sigma()*1E-3; si_y = ft.Y0Sigma()*1E-3;
  si_x = 19.05E-6; 
  si_y = 19.05E-6;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

RPElasticReconstruction::optics_type::optics_type(LHCOpticsApproximator *o, const BeamOpticsParams &optPar, const bool invert)
{
  double D = 0.;
  double in_null[5] = {optPar.GetBeamDisplacementX(), optPar.GetCrossingAngleX(), optPar.GetBeamDisplacementY(), optPar.GetCrossingAngleY(), 0.};
  o->GetLinearApproximation(in_null, Cx, Lx, vx, Cy, Ly, vy, D);
  if (invert) {
    Lx = -Lx;
    Ly = -Ly;
  }
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

RPElasticReconstruction::RPElasticReconstruction(const edm::ParameterSet& conf) : eec_graph(NULL), opticsIP(0., 0., 1., 1.)
{
  verbosity = conf.getUntrackedParameter<unsigned int>("verbosity", 1);
  addIPHit = conf.getParameter<unsigned int>("addIPHit");

  fittedTrackCollectionInputLabel = conf.getParameter<edm::InputTag>("RPFittedTrackCollLabel");

  roadsize_x = conf.getParameter<double>("roadsize_x");
  roadsize_y = conf.getParameter<double>("roadsize_y");
  angleTolerance_x = conf.getParameter<double>("angleTolerance_x");
  angleTolerance_y = conf.getParameter<double>("angleTolerance_y");
  vertexTolerance_x = conf.getParameter<double>("vertexTolerance_x") * 1E-3;  // mm -> m
  vertexTolerance_y = conf.getParameter<double>("vertexTolerance_y") * 1E-3;

  try {
    const ParameterSet& eec = conf.getParameter<ParameterSet>("edgeEffectCorrection");
    eec_enabled = eec.getParameter<bool>("enabled");
    eec_si_bd = eec.getParameter<double>("si_bd");
    eec_th_L = eec.getParameter<double>("th_L");
  }
  catch (...) {
    eec_enabled = false;
  }
  if (eec_enabled) InitializeEdgeEffectCorrection();

  produces<RPRecoElasticEvent>();
}

//----------------------------------------------------------------------------------------------------

RPElasticReconstruction::~RPElasticReconstruction()
{
  if (eec_graph) delete eec_graph;
}

//----------------------------------------------------------------------------------------------------

void RPElasticReconstruction::InitializeOpticalFunctions(const ProtonTransportFunctions &optFun, const BeamOpticsParams &optPar)
{
  for (ProtonTransportFunctions::MapType::const_iterator it = optFun.GetFunctionMap().begin(); it != optFun.GetFunctionMap().end(); ++it) {
      bool invert = (TotRPDetId::ArmOfRP(it->first) == 0);
        optics[it->first] = optics_type(it->second.real, optPar, invert);
  }

  if (verbosity > 0) {
    printf(">> RPElasticReconstruction::InitializeOpticalFunctions: optics.size = %lu\n", optics.size());
    for (map<unsigned int, optics_type>::iterator it = optics.begin(); it != optics.end(); ++it) {
      printf("ID %i: Cx = %.0f um, Lx = %.2f m, vx = %.2f, Cy = %.0f um, Ly = %.2f m, vy = %.2f\n",
        it->first, it->second.Cx*1E6, it->second.Lx, it->second.vx, it->second.Cy*1E6, it->second.Ly, it->second.vy);
    }
  }
}

//----------------------------------------------------------------------------------------------------

void RPElasticReconstruction::beginRun(edm::Run const&, edm::EventSetup const& es)
{
  // get optics and beam parameters
  edm::ESHandle<BeamOpticsParams> optPar;
  es.get<BeamOpticsParamsRcd>().get(optPar);

  edm::ESHandle<ProtonTransportFunctions> optFun;
  es.get<ProtonTransportRcd>().get(optFun);

  InitializeOpticalFunctions(*optFun, *optPar);
}

//----------------------------------------------------------------------------------------------------

void RPElasticReconstruction::produce(edm::Event& e, const edm::EventSetup& c)
{
  // get optics and beam parameters
  edm::ESHandle<BeamOpticsParams> optPar;
  c.get<BeamOpticsParamsRcd>().get(optPar);

  // get intput
  edm::Handle< RPFittedTrackCollection > input; 
  e.getByLabel(fittedTrackCollectionInputLabel, input);

  if (verbosity > 5 || ((e.id().event() - 1) % 100) == 0)
    cout << endl << "EVENT " << e.id().event() << endl;

  // prepare output
  std::auto_ptr<RPRecoElasticEvent> output(new RPRecoElasticEvent());

  // processing
  if (MakeRoads(input.product(), output.get(), *optPar) == 0)
    if (MakeFits(input.product(), output.get(), *optPar) == 0)
      if (MakeResult(input.product(), output.get(), *optPar) == 0)
        output->status = output->sOK;

  if (verbosity > 0 && verbosity <= 5 && !output->isValid() && output->preferredRoad >= 0)
    cout << "Event " << (e.id().event() - 1) << " rejected." << endl;

  if (verbosity > 5) {
    cout << "** valid = " << output->isValid() << endl;
    if (output->isValid())
      printf("** product:\n\ttheta*_x = %.2E\n\ttheta*_y = %.2E\n\tx* = %.2E\n\ty* = %.2E\n",
          output->result.th_x, output->result.th_y, output->result.x, output->result.y);
  }


  // save output
  e.put(output);
}

//----------------------------------------------------------------------------------------------------

unsigned char RPElasticReconstruction::MakeRoads(const RPFittedTrackCollection *input, RPRecoElasticEvent *output, const BeamOpticsParams &optPar)
{
  if (verbosity > 5)
    cout << "*** MakeRoads (x and y in mm, al in urad)" << endl;

  // go through all hits
  for (map<RPId, RPFittedTrack>::const_iterator it = input->begin(); it != input->end(); ++it) {
    // check validity
    if (!it->second.IsValid())
      continue;

    // compute rough angle (alpha)
    optics_type &of = GetOptics(it->first);
    double al_x = 0., al_y = 0.;
    if (of.Lx != 0.)
      al_x = (it->second.X0()*1E-3 - of.Cx) / of.Lx;  // from mm to m
    if (of.Ly != 0.)
      al_y = (it->second.Y0()*1E-3 - of.Cy) / of.Ly;
    if (verbosity > 5)
      printf("\tID = %03i: x = %.3f, y = %.3f, al_x = %.2E, al_y = %.2E,\n", it->first, it->second.X0(), it->second.Y0(), al_x*1E6, al_y*1E6);

    // look for a road
    bool found = false;
    for (unsigned int i = 0; i < output->roads.size(); i++) {
      if ( fabs(output->roads[i].centerX() - al_x) <= roadsize_x && fabs(output->roads[i].centerY() - al_y) <= roadsize_y ) {
        if (verbosity > 5)
          printf("\t\tadding to road #%i\n", i);
        found = true;
        output->roads[i].sumx += al_x;
        output->roads[i].sumy += al_y;
        output->roads[i].minx = min(output->roads[i].minx, al_x);
        output->roads[i].maxx = max(output->roads[i].maxx, al_x);
        output->roads[i].miny = min(output->roads[i].miny, al_y);
        output->roads[i].maxy = max(output->roads[i].maxy, al_y);
        output->roads[i].members.push_back(it->first);
      }
    }

    // if not found add a new road
    if (!found) {
      if (verbosity > 5)
        printf("\t\tadding new road\n");
      RPRecoElasticEvent::road_type road;
      road.sumx = al_x;
      road.sumy = al_y;
      road.minx = road.maxx = al_x;
      road.miny = road.maxy = al_y;
      road.members.push_back(it->first);
      output->roads.push_back(road);
    }
  }

  // analyze roads
  unsigned maxSize = 0;
  for (unsigned int i = 0; i < output->roads.size(); i++) {
    RPRecoElasticEvent::road_type& road = output->roads[i];
    vector<unsigned short> &members = road.members;

    // sum number of RP left and right
    unsigned int left = 0, right = 0;
    for (unsigned int j = 0; j < members.size(); j++) {
      if (IsRPLeft(members[j]))
        left++;
      else
        right++;
    }

    // print out
    if (verbosity > 5) {
      printf("ROAD #%i\n", i);
      printf("\tcenter x = %.1f urad\n", road.centerX()*1E6);
      printf("\tcenter y = %.1f urad\n", road.centerY()*1E6);
      printf("\tleft = %i, right = %i\n", left, right);
      printf("\tmembers = ");
      for (unsigned int j = 0; j < members.size(); j++)
        printf("%i, ", members[j]);
      printf("\n");
    }

    // continue only if left >= 1 and right >=1
    if (left < 1 || right < 1) continue;

    // check the maximum size
    if (output->preferredRoad < 0 || maxSize < members.size()) {
      output->preferredRoad = i;
      maxSize = members.size();
    }
  }

  // set status and print out
  if (output->roads.size() == 0) {
    output->status = output->sNoRoad;
    if (verbosity > 1) printf("no road\n");
    return 1;
  }

  if (output->preferredRoad < 0) {
    output->status = output->sNoGoodRoad;
    if (verbosity > 1) printf("preferred road: %i\n", output->preferredRoad);
    return 2;
  }

  return 0;
}

//----------------------------------------------------------------------------------------------------

RPRecoElasticEvent::fit_type RPElasticReconstruction::DoFit(std::vector<recoHit_type> &hits)
{
  RPRecoElasticEvent::fit_type output;

  double sLL_x, svv_x, svL_x, sav_x, saL_x, saa_x;
  double sLL_y, svv_y, svL_y, sav_y, saL_y, saa_y;
  sLL_x = svv_x = svL_x = sav_x = saL_x = saa_x = 0.;
  sLL_y = svv_y = svL_y = sav_y = saL_y = saa_y = 0.;

  // loop over hits
  for (unsigned int i = 0; i < hits.size(); i++) {
    if (verbosity > 5)
      printf("\tx = %.3f, y = %.3f, si_x = %.3f, si_y = %.3f | Lx = %.2f, Ly = %.2f, vx = %.2E, vy = %.2E\n",
          hits[i].x*1E3, hits[i].y*1E3, hits[i].si_x*1E3, hits[i].si_y*1E3, hits[i].of->Lx, hits[i].of->Ly, hits[i].of->vx, hits[i].of->vy);

    recoHit_type &h = hits[i];

    sLL_x += h.of->Lx * h.of->Lx / h.si_x / h.si_x;
    svv_x += h.of->vx * h.of->vx / h.si_x / h.si_x;
    svL_x += h.of->vx * h.of->Lx / h.si_x / h.si_x;
    sav_x += h.x * h.of->vx / h.si_x / h.si_x;
    saL_x += h.x * h.of->Lx / h.si_x / h.si_x;
    saa_x += h.x * h.x / h.si_x / h.si_x;

    sLL_y += h.of->Ly * h.of->Ly / h.si_y / h.si_y;
    svv_y += h.of->vy * h.of->vy / h.si_y / h.si_y;
    svL_y += h.of->vy * h.of->Ly / h.si_y / h.si_y;
    sav_y += h.y * h.of->vy / h.si_y / h.si_y;
    saL_y += h.y * h.of->Ly / h.si_y / h.si_y;
    saa_y += h.y * h.y / h.si_y / h.si_y;
  }

  double det_x = svv_x * sLL_x - svL_x * svL_x;
  double det_y = svv_y * sLL_y - svL_y * svL_y;

  if (det_x != 0.) {
    output.th_x = (svv_x * saL_x - svL_x * sav_x) / det_x;
    output.x = (sLL_x * sav_x - svL_x * saL_x) / det_x;
    output.si_th_x = sqrt( svv_x / det_x );
    output.si_x = sqrt( sLL_x / det_x );
    output.ndf_x = hits.size() - 2;
    output.s2min_x = saa_x - output.th_x * saL_x - output.x * sav_x;
  } else
    if (verbosity > 0) printf("Zero determinant: det_x = %f\n", det_x);

  if (det_y != 0.) {
    output.th_y = (svv_y * saL_y - svL_y * sav_y) / det_y;
    output.y = (sLL_y * sav_y - svL_y * saL_y) / det_y;
    output.si_th_y = sqrt( svv_y / det_y );
    output.si_y = sqrt( sLL_y / det_y );
    output.ndf_y = hits.size() - 2;
    output.s2min_y = saa_y - output.th_y * saL_y - output.y * sav_y;
  } else
    if (verbosity > 0) printf("Zero determinant: det_y = %f\n", det_y);

  return output;
}

//----------------------------------------------------------------------------------------------------

unsigned char RPElasticReconstruction::MakeFits(const RPFittedTrackCollection *input, RPRecoElasticEvent *output, const BeamOpticsParams &optPar)
{
  if (verbosity > 5) cout << "*** MakeFits" << endl;

  if (output->preferredRoad < 0) return 1;

  // build hits collection
  std::vector<recoHit_type> leftHits, rightHits, allHits;
  for (unsigned int i = 0; i <output->roads[output->preferredRoad].members.size(); i++) {
    unsigned short ID = output->roads[output->preferredRoad].members[i];

    recoHit_type recoHit(input->find(ID)->second, &optics[ID]);

    if (IsRPLeft(ID))  leftHits.push_back(recoHit);
    else        rightHits.push_back(recoHit);
    allHits.push_back(recoHit);
  }

  // add hit in origin
  recoHit_type hitOrigin(0., optPar.GetPrimVertSizeX(), 0., optPar.GetPrimVertSizeY(), &opticsIP);
  if (addIPHit > 0 || leftHits.size() == 1) leftHits.push_back(hitOrigin);
  if (addIPHit > 0 || rightHits.size() == 1) rightHits.push_back(hitOrigin);
  if (addIPHit > 1) allHits.push_back(hitOrigin);

  // perform fits
  if (verbosity > 5) cout << "* hit info in mm, eff. length in m" << endl;
  if (verbosity > 5) cout << "* left fit" << endl;
  output->leftFit = DoFit(leftHits);
  if (verbosity > 5) printf("\t\ttheta_x = %.1f urad, theta_y = %.1f urad\n", output->leftFit.th_x*1E6, output->leftFit.th_y*1E6);
  if (verbosity > 5) cout << "* right fit" << endl;
  output->rightFit = DoFit(rightHits);
  if (verbosity > 5) printf("\t\ttheta_x = %.1f urad, theta_y = %.1f urad\n", output->rightFit.th_x*1E6, output->rightFit.th_y*1E6);
  if (verbosity > 5) cout << "* global fit" << endl;
  output->globalFit = DoFit(allHits);
  if (verbosity > 5) printf("\t\ttheta_x = %.1f urad, theta_y = %.1f urad\n", output->globalFit.th_x*1E6, output->globalFit.th_y*1E6);

  return 0;
}

//----------------------------------------------------------------------------------------------------

unsigned char RPElasticReconstruction::MakeResult(const RPFittedTrackCollection *input, RPRecoElasticEvent *output, const BeamOpticsParams &optPar)
{
  if (verbosity > 5) cout << "*** MakeResult" << endl;

  // apply cuts
  if ( fabs(output->leftFit.th_x - output->rightFit.th_x) > angleTolerance_x ) {
    if (verbosity > 0) printf("Inconsistent theta_x: left = %.2E, right = %.2E\n", output->leftFit.th_x, output->rightFit.th_x);
    output->status = output->sRejected;
    output->rejectReason |= output->rAngleX;
  }

  if ( fabs(output->leftFit.th_y - output->rightFit.th_y) > angleTolerance_y ) {
    if (verbosity > 0) printf("Inconsistent theta_y: left = %.2E, right = %.2E\n", output->leftFit.th_y, output->rightFit.th_y);
    output->status = output->sRejected;
    output->rejectReason |= output->rAngleY;
  }

  if ( fabs(output->leftFit.x - output->rightFit.x) >  vertexTolerance_x ) {
    if (verbosity > 0) printf("Inconsistent x*: left = %.2E, right = %.2E\n", output->leftFit.x, output->rightFit.x);
    output->status = output->sRejected;
    output->rejectReason |= output->rVertexX;  
  }

  if ( fabs(output->leftFit.y - output->rightFit.y) >  vertexTolerance_y ) {
    if (verbosity > 0) printf("Inconsistent y*: left = %.2E, right = %.2E\n", output->leftFit.y, output->rightFit.y);
    output->status = output->sRejected;
    output->rejectReason |= output->rVertexY;  
  }

  // make result
  output->result = output->globalFit;

  // apply correction
  if (eec_enabled) ApplyEdgeEffectCorrection(output->result);

  return output->rejectReason;
}

//----------------------------------------------------------------------------------------------------

void RPElasticReconstruction::InitializeEdgeEffectCorrection()
{
  if (verbosity > 2)
    printf(">> RPElasticReconstruction::InitializeEdgeEffectCorrection > si_bd = %.2f urad, th_L = %.2f urad\n", eec_si_bd*1E6, eec_th_L*1E6);

  if (!eec_graph) eec_graph = new TGraph();
  else eec_graph->Set(0);

  const double range = 500E-6;
  for (double x = -range; x < range; x += 1E-6) {
    double z1 = (eec_th_L - x) / eec_si_bd;
    double z2 = (eec_th_L + x) / eec_si_bd;
    double si_eff = eec_si_bd * sqrt(2.);
    double y = si_eff / sqrt(TMath::Pi()) * (exp(-z1*z1/2.) - exp(-z2*z2/2.)) / (TMath::Erfc(z1/sqrt(2.)) + TMath::Erfc(z2/sqrt(2.)));
    eec_graph->SetPoint(eec_graph->GetN(), y + x, x);
  }
}

//----------------------------------------------------------------------------------------------------

void RPElasticReconstruction::ApplyEdgeEffectCorrection(RPRecoElasticEvent::fit_type& fit)
{
  fit.th_y = eec_graph->Eval(fit.th_y);  
}

DEFINE_FWK_MODULE(RPElasticReconstruction);


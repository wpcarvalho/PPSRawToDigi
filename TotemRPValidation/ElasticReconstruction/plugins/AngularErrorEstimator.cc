/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*	Jan Kašpar (jan.kaspar@gmail.com)
*    
* $$RCSfile: AngularErrorEstimator.cc,v $: $
* $Revision: 9977 $
* $Date: 2015-01-12 15:00:26 +0100 (Mon, 12 Jan 2015) $
*
****************************************************************************/


#include <stdio.h>
#include <iostream>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include "TotemCondFormats/ProtonTransportFunctions/interface/ProtonTransportFunctions.h"
#include "SimG4CMS/TotemRPProtTranspPar/interface/LHCOpticsApproximator.h"

//----------------------------------------------------------------------------------------------------
/**
 *\brief To estimate angular resolution of the elastic reconstruction procedure.
**/
class AngularErrorEstimator : public edm::EDAnalyzer
{
 public:
  AngularErrorEstimator(const edm::ParameterSet&);
  ~AngularErrorEstimator();

  struct sums_type {
    unsigned int N;
    double vv, LL, vL;
    sums_type() : N(0), vv(0.), LL(0.), vL(0.) {};
    void add(double v, double L) { N++; vv+= v*v; LL += L*L; vL += v*L; }
    double den() { return vv*LL - vL*vL; }
    void Print() { printf("N = %i, vv = %.2E, LL = %.2E, vL = %.2E, den = %.2E\n", N, vv, LL, vL, den()); }
  };

  sums_type slx, sly, srx, sry, sgx, sgy;

  enum station_activity { saNone, saVer, saHor, saBoth } st12, st10, st02, st00;

 private:
  double si_p;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&) {}
  virtual void endJob() {}

  void ParseStationActivity(station_activity &sa, const std::string &str);
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

using namespace edm;
using namespace std;

AngularErrorEstimator::AngularErrorEstimator(const edm::ParameterSet& conf)
{
  si_p = conf.getParameter<double>("si_p");
  
  string s;
  s = conf.getParameter<string>("station_12"); ParseStationActivity(st12, s);
  s = conf.getParameter<string>("station_10"); ParseStationActivity(st10, s);
  s = conf.getParameter<string>("station_00"); ParseStationActivity(st00, s);
  s = conf.getParameter<string>("station_02"); ParseStationActivity(st02, s);
}

//----------------------------------------------------------------------------------------------------

void AngularErrorEstimator::ParseStationActivity(station_activity &sa, const std::string &s)
{
  sa = saNone;
  if (s == "v") sa = saVer;
  if (s == "h") sa = saHor;
  if (s == "vh") sa = saBoth;
  if (s == "hv") sa = saBoth;
  printf("%i\n", sa);
}

//----------------------------------------------------------------------------------------------------

AngularErrorEstimator::~AngularErrorEstimator()
{
}

//----------------------------------------------------------------------------------------------------

void Round(double &x)
{
  x = floor(x*100)/100;
}


void AngularErrorEstimator::beginRun(edm::Run const&, edm::EventSetup const& es)
{
  // get optics and beam parameters
  edm::ESHandle<BeamOpticsParams> optPar;
  es.get<BeamOpticsParamsRcd>().get(optPar);
  if(!optPar.isValid())
    throw cms::Exception("ElasticRecoValLibrary::analyze") << " edm::ESHandle<BeamOpticsParams> is invalid";

  edm::ESHandle<ProtonTransportFunctions> optFun;
  es.get<BeamOpticsParamsRcd>().get(optFun);
  if(!optFun.isValid())
	throw cms::Exception("RPElasticReconstruction::beginJob") << " edm::ESHandle<ProtonTransportFunctions> is invalid";

  printf("\n>> selected detectors\n");
  printf("ID      Cx          Lx        vx      Cy          Ly        vy\n");
  const ProtonTransportFunctions::MapType& fncts = optFun->GetFunctionMap();
  for (ProtonTransportFunctions::MapType::const_iterator it = fncts.begin(); it != fncts.end(); ++it) {
    unsigned int ID = it->first;
    bool left = (ID < 100);
    unsigned int station = ID / 10;
    unsigned int RP = ID % 10;

    // select geometry
    if (station == 12) {
      if ((RP == 0 || RP == 1) && (st12 == saNone || st12 == saHor)) continue;
      if ((RP == 2 || RP == 3) && (st12 == saNone || st12 == saVer)) continue;
      if ((RP == 4 || RP == 5) && (st12 == saNone || st12 == saHor)) continue;
    }
    if (station == 10) {
      if ((RP == 0 || RP == 1) && (st10 == saNone || st10 == saHor)) continue;
      if ((RP == 2 || RP == 3) && (st10 == saNone || st10 == saVer)) continue;
      if ((RP == 4 || RP == 5) && (st10 == saNone || st10 == saHor)) continue;
    }
    if (station == 0) {
      if ((RP == 0 || RP == 1) && (st00 == saNone || st00 == saHor)) continue;
      if ((RP == 2 || RP == 3) && (st00 == saNone || st00 == saVer)) continue;
      if ((RP == 4 || RP == 5) && (st00 == saNone || st00 == saHor)) continue;
    }
    if (station == 2) {
      if ((RP == 0 || RP == 1) && (st02 == saNone || st02 == saHor)) continue;
      if ((RP == 2 || RP == 3) && (st02 == saNone || st02 == saVer)) continue;
      if ((RP == 4 || RP == 5) && (st02 == saNone || st02 == saHor)) continue;
    }

    // not to count twice the vertical detectors - they never work at the same time and they have the same opt. functions
    if (RP == 1 || RP == 5) continue; 

    // extract optical functions
    double Cx = 0, Cy = 0, Lx = 0, Ly = 0, vx = 0, vy = 0, D = 0.;
    double in_null[5] = {optPar->GetBeamDisplacementX(), optPar->GetCrossingAngleX(), optPar->GetBeamDisplacementY(), optPar->GetCrossingAngleY(), 0.};
    it->second.real->GetLinearApproximation(in_null, Cx, Lx, vx, Cy, Ly, vy, D);
    if (left) {
		Lx = -Lx;
		Ly = -Ly;
	}

    Round(Lx);
    Round(Ly);
    Round(vx);
    Round(vy);
    printf("%3i:%12.2E%9.2f%9.2f%12.2E%9.2f%9.2f\n", ID, Cx, Lx, vx, Cy, Ly, vy);

    // add to sums
    if (left) { slx.add(vx, Lx); sly.add(vy, Ly); }
    else      { srx.add(vx, Lx); sry.add(vy, Ly); }
    sgx.add(vx, Lx); sgy.add(vy, Ly);
  }

  /*
  printf("\n>> detectors left: %i, right: %i, global: %i\n", slx.N, srx.N, sgx.N);
  printf(">> LL: slx = %.2E, sly = %.2E, srx = %.2E, sry = %.2E, sgx = %.2E, sgy = %.2E\n", slx.LL, sly.LL, srx.LL, sry.LL, sgx.LL, sgy.LL);
  printf(">> denominators: slx = %.2E, sly = %.2E, srx = %.2E, sry = %.2E, sgx = %.2E, sgy = %.2E\n", slx.den(), sly.den(), srx.den(), sry.den(), sgx.den(), sgy.den());
  */

  printf("\n>> summary of collected data:\n");
  printf("slx: "); slx.Print();
  printf("sly: "); sly.Print();
  printf("srx: "); srx.Print();
  printf("sry: "); sry.Print();
  printf("sgx: "); sgx.Print();
  printf("sgy: "); sgy.Print();

  double si_bd = optPar->GetBeamDivergenceX();
  printf("\n>> beam divergence in projection: %.2f urad\n", si_bd*1E6);
  printf(">> pitch error std. deviation: %.2f um\n", si_p*1E6);

  printf("\n\\multispan3\\bvrule\\hfil\\hbox{\\strut    left arm    }\\hfil\\cr\\bln\n");
  printf("\\hbox{beam divergence} & %6.2f & %6.2f \\cr\\ln\n", si_bd * 1E6, si_bd*1E6);
  double p_comp_x = sqrt(slx.vv / slx.den()) * si_p;
  double p_comp_y = sqrt(sly.vv / sly.den()) * si_p;
  printf("\\hbox{pitch}           & %6.2f & %6.2f \\cr\\ln\n", p_comp_x * 1E6, p_comp_y * 1E6);
  double full_x = sqrt(si_bd*si_bd + p_comp_x*p_comp_x);
  double full_y = sqrt(si_bd*si_bd + p_comp_y*p_comp_y);
  printf("\\hbox{full}            & %6.2f & %6.2f \\cr\\bln\n", full_x * 1E6, full_y * 1E6);

  printf("\\multispan3\\bvrule\\hfil\\hbox{\\strut    right arm    }\\hfil\\cr\\bln\n");
  printf("\\hbox{beam divergence} & %6.2f & %6.2f \\cr\\ln\n", si_bd * 1E6, si_bd*1E6);
  p_comp_x = sqrt(srx.vv / srx.den()) * si_p;
  p_comp_y = sqrt(sry.vv / sry.den()) * si_p;
  printf("\\hbox{pitch}           & %6.2f & %6.2f \\cr\\ln\n", p_comp_x * 1E6, p_comp_y * 1E6);
  full_x = sqrt(si_bd*si_bd + p_comp_x*p_comp_x);
  full_y = sqrt(si_bd*si_bd + p_comp_y*p_comp_y);
  printf("\\hbox{full}            & %6.2f & %6.2f \\cr\\bln\n", full_x * 1E6, full_y * 1E6);

  printf("\\multispan3\\bvrule\\hfil\\hbox{\\strut    double arm    }\\hfil\\cr\\bln\n");
  double fR_x = (sgx.vv*srx.LL - sgx.vL*srx.vL) / sgx.den(), fR_y = (sgy.vv*sry.LL - sgy.vL*sry.vL) / sgy.den();
  double fL_x = (sgx.vv*slx.LL - sgx.vL*slx.vL) / sgx.den(), fL_y = (sgy.vv*sly.LL - sgy.vL*sly.vL) / sgy.den();
  printf("%%x ratios: fR = %.3f, fL = %.3f, sum = %.3f\n", fR_x, fL_x, fR_x + fL_x);
  printf("%%y ratios: fR = %.3f, fL = %.3f, sum = %.3f\n", fR_y, fL_y, fR_y + fL_y);
  double bd_comp_x = si_bd * sqrt(fR_x*fR_x + fL_x*fL_x);
  double bd_comp_y = si_bd * sqrt(fR_y*fR_y + fL_y*fL_y);
  printf("\\hbox{beam divergence} & %6.2f & %6.2f \\cr\\ln\n", bd_comp_x * 1E6, bd_comp_y*1E6);
  p_comp_x = sqrt(sgx.vv / sgx.den()) * si_p;
  p_comp_y = sqrt(sgy.vv / sgy.den()) * si_p;
  printf("\\hbox{pitch}           & %6.2f & %6.2f \\cr\\ln\n", p_comp_x * 1E6, p_comp_y * 1E6);
  full_x = sqrt(bd_comp_x*bd_comp_x + p_comp_x*p_comp_x);
  full_y = sqrt(bd_comp_x*bd_comp_x + p_comp_y*p_comp_y);
  printf("\\hbox{full}            & %6.2f & %6.2f \\cr\\bln\n", full_x * 1E6, full_y * 1E6);
  
  printf("\\multispan3\\bvrule\\hfil\\hbox{\\strut    right-left differences    }\\hfil\\cr\\bln\n");
  double bd_comp = si_bd * sqrt(2.);
  p_comp_x = sqrt(srx.vv / srx.den() + slx.vv / slx.den()) * si_p;
  p_comp_y = sqrt(sry.vv / sry.den() + sly.vv / sly.den()) * si_p;
  full_x = sqrt(bd_comp*bd_comp + p_comp_x*p_comp_x);
  full_y = sqrt(bd_comp*bd_comp + p_comp_y*p_comp_y);
  printf("\\hbox{beam divergence} & %6.2f & %6.2f \\cr\\ln\n", bd_comp * 1E6, bd_comp*1E6);
  printf("\\hbox{pitch}           & %6.2f & %6.2f \\cr\\ln\n", p_comp_x * 1E6, p_comp_y * 1E6);
  printf("\\hbox{full}            & %6.2f & %6.2f \\cr\\bln\n", full_x * 1E6, full_y * 1E6);
}

DEFINE_FWK_MODULE(AngularErrorEstimator);


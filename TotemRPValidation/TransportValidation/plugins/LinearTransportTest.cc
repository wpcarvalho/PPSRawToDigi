/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $$RCSfile: LinearTransportTest.cc,v $: $
* $Revision: 1.1 $
* $Date: 2009/02/04 16:19:57 $
*
****************************************************************************/

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include "TotemCondFormats/ProtonTransportFunctions/interface/ProtonTransportFunctions.h"
#include "TotemCondFormats/DataRecord/interface/ProtonTransportRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "SimG4CMS/TotemRPProtTranspPar/interface/LHCOpticsApproximator.h"

#include "TFile.h"
#include "TGraph.h"

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

/**
 *\brief A template for profile data structures
 **/
template <class bin_type>
class Profile
{
  public:
    std::string name;

    Profile(std::string _name = "") : name(_name) {}
    Profile(std::string _name, unsigned int _bins, double _from, double _to) : name(_name), from(_from), to(_to) {
      for (unsigned int i = 0; i < _bins; i++)
        bins.push_back(bin_type());
    }

    void Fill(double x, double y, double w = 1.) {
      if (bins.size() < 1)
        return;
      double idd = floor((x - from)/(to - from) * bins.size());
      if (idd < 0 || idd > (bins.size() - 1))
        return;
      bins[int(idd)].Fill(y, w);
    }

    double GetBinCenter(unsigned int i) const { return from + (to - from)/bins.size() * (double(i) + 0.5); }
    const std::vector<bin_type>& GetBins() const { return bins; }

  protected:
    double from, to;
    std::vector<bin_type> bins;
};

//----------------------------------------------------------------------------------------------------

/**
 *\brief Data set for one bin
 **/
struct bin_type
{
  double S1, Sy, Syy;
  double max_abs_y;

  bin_type() : S1(0.), Sy(0.), Syy(0.), max_abs_y(0.) {}
  void Fill(double y, double ) {
    S1 += 1.; Sy += y; Syy += y*y;
    if (fabs(y) > max_abs_y)
      max_abs_y = fabs(y);
  }
  double Mean() const { return Sy / S1; }
  double Sigma() const { return sqrt((Syy - Sy*Sy/S1) / (S1 - 1.)); }
};

//----------------------------------------------------------------------------------------------------

/**
 *\brief To check the validity of the linear approximation to the proton transport parameterization.
**/
class LinearTransportTest : public edm::EDAnalyzer
{
  public:
    LinearTransportTest(const edm::ParameterSet &ps);
    ~LinearTransportTest() {}

  private:
    unsigned int verbosity;

    double th_min, th_max;
    double vertex_n_sigma;
    double xi_min, xi_max;
    unsigned int N;
    std::vector<unsigned int> RPs;
    std::string outputFile;

    edm::ESHandle<BeamOpticsParams> optPar;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void analyze(const edm::Event &e, const edm::EventSetup &es);
    virtual void endJob();

    void DoScan(unsigned int RP, LHCOpticsApproximator *);
    void WriteProfiles(const Profile<bin_type>&);
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

LinearTransportTest::LinearTransportTest(const edm::ParameterSet &ps)
{
  verbosity = ps.getUntrackedParameter<unsigned int>("verbosity", 0);

  // conversion to theta is done in beginJob
  th_min = ps.getParameter<double>("t_min");
  th_max = ps.getParameter<double>("t_max");

  vertex_n_sigma = ps.getParameter<double>("vertex_n_sigma");
  xi_min = ps.getParameter<double>("xi_min");
  xi_max = ps.getParameter<double>("xi_max");

  N = ps.getParameter<unsigned int>("samples");
  if (N < 2) throw cms::Exception("LinearTransportTest::LinearTransportTest") << " The number of samples (N) must be greater than 1" << endl;

  RPs = ps.getParameter< vector<unsigned int> >("RPs");

  outputFile = ps.getParameter<string>("outputFile");
}

//----------------------------------------------------------------------------------------------------

void LinearTransportTest::beginRun(edm::Run const&, edm::EventSetup const& es)
{
  es.get<BeamOpticsParamsRcd>().get(optPar);
  if(!optPar.isValid())
    throw cms::Exception("LinearTransportTest::beginJob") << " edm::ESHandle<BeamOpticsParams> is invalid";

  edm::ESHandle<ProtonTransportFunctions> optFun;
  es.get<ProtonTransportRcd>().get(optFun);
  if(!optFun.isValid())
	throw cms::Exception("LinearTransportTest::beginJob") << " edm::ESHandle<ProtonTransportFunctions> is invalid";

  // convert t to theta
  th_min = sqrt(th_min) / optPar->GetBeamMomentum(); 
  th_max = sqrt(th_max) / optPar->GetBeamMomentum();

  TFile *of = new TFile(outputFile.c_str(), "recreate");

  for (unsigned int i = 0; i < RPs.size(); i++) {
    LHCOpticsApproximator *f = optFun->GetFunction(RPs[i]);

    char buf[20];
    sprintf(buf, "%i", RPs[i]);
    gDirectory = gDirectory->mkdir(buf);

    DoScan(RPs[i], f);   

    gDirectory->cd("..");
  }

  delete of;
}

//----------------------------------------------------------------------------------------------------

void LinearTransportTest::DoScan(unsigned int RP, LHCOpticsApproximator *o)
{
  printf("scanning %i\n", RP);

  // get the linear approximation
  double Cx = 0., Lx = 0., vx = 0., Cy = 0., Ly = 0., vy = 0., D = 0.;
  double in_null[5] = {optPar->GetBeamDisplacementX(), optPar->GetCrossingAngleX(), optPar->GetBeamDisplacementY(), 
    optPar->GetCrossingAngleY(), (xi_min + xi_max) / 2.};
  o->GetLinearApproximation(in_null, Cx, Lx, vx, Cy, Ly, vy, D);

  // prepare division
  double dth = (th_max - th_min) / (N - 1);
  double xs_range = optPar->GetPrimVertSizeX() * vertex_n_sigma;
  double dxs = 2*xs_range / (N - 1);
  double ys_range = optPar->GetPrimVertSizeY() * vertex_n_sigma;
  double dys = 2*ys_range / (N - 1);
  double dxi = (xi_max - xi_min) / (N - 1);
  //
  // book histograms
  TH1D *diff_x = new TH1D("diff_x", ";x_{lin} - x_{par}   (#mum)", 100, -25., 25.);
  TH1D *diff_y = new TH1D("diff_y", ";y_{lin} - y_{par}   (#mum)", 100, -25., 25.);

  Profile<bin_type> px_thx("px_thx", N*N, -th_max - dth/2, th_max + dth/2); // ... this is a dirty trick
  Profile<bin_type> px_thy("px_thy", N*N, -th_max - dth/2, th_max + dth/2);
  Profile<bin_type> px_xs("px_xs", N, -xs_range - dxs/2, +xs_range + dxs/2);
  Profile<bin_type> px_ys("px_ys", N, -ys_range - dys/2, +ys_range + dys/2);
  Profile<bin_type> px_xi("px_xi", N, xi_min - dxi/2, xi_max + dxi/2);

  Profile<bin_type> py_thx("py_thx", N*N, -th_max - dth/2, th_max + dth/2); // ... this is a dirty trick
  Profile<bin_type> py_thy("py_thy", N*N, -th_max - dth/2, th_max + dth/2);
  Profile<bin_type> py_xs("py_xs", N, -xs_range - dxs/2, +xs_range + dxs/2);
  Profile<bin_type> py_ys("py_ys", N, -ys_range - dys/2, +ys_range + dys/2);
  Profile<bin_type> py_xi("py_xi", N, xi_min - dxi/2, xi_max + dxi/2);
  
  // scan phase space
  for (double th_x = -th_max; th_x <= th_max; th_x += dth) {
    if (th_x > -th_min && th_x < 0.) th_x = th_min;
    for (double th_y = th_min; th_y <= th_max; th_y += dth) {
      if (th_y > -th_min && th_y < 0.) th_y = th_min;
      for (double xs = -xs_range; xs <= xs_range; xs += dxs)
        for (double ys = -ys_range; ys <= ys_range; ys += dys)
          for (double xi = xi_min; xi <= xi_max; xi += dxi) {
            // evaluate linear approximation
            double x_lin = Cx + Lx * th_x + vx * xs + D*xi;
            double y_lin = Cy + Ly * th_y + vy * ys;

            // evaluate exact position
            double in[5] = {optPar->GetBeamDisplacementX() + xs, optPar->GetCrossingAngleX() + th_x, 
                            optPar->GetBeamDisplacementY() + ys, optPar->GetCrossingAngleY() + th_y, xi};

            double out[2] = {0., 0.};
            o->Transport2D(in, out, false);
            double x_par = out[0];
            double y_par = out[1];

            // print out
            if (verbosity > 5)
              printf("%.1f, %.1f, %.1f, %.1f, %.1E --> %E, %E | %E, %E\n", th_x*1E6, th_y*1E6, xs*1E6, ys*1E6, xi, x_lin, x_par, y_lin, y_par);

            // comparison
            diff_x->Fill((x_lin - x_par) * 1E6);
            diff_y->Fill((y_lin - y_par) * 1E6);
            px_thx.Fill(th_x, (x_lin - x_par) * 1E6);
            px_thy.Fill(th_y, (x_lin - x_par) * 1E6);
            px_xs.Fill(xs, (x_lin - x_par) * 1E6);
            px_ys.Fill(ys, (x_lin - x_par) * 1E6);
            px_xi.Fill(xi, (x_lin - x_par) * 1E6);

            py_thx.Fill(th_x, (y_lin - y_par) * 1E6);
            py_thy.Fill(th_y, (y_lin - y_par) * 1E6);
            py_xs.Fill(xs, (y_lin - y_par) * 1E6);
            py_ys.Fill(ys, (y_lin - y_par) * 1E6);
            py_xi.Fill(xi, (y_lin - y_par) * 1E6);
          }
    }
  }

  // write histograms
  diff_x->Write();
  diff_y->Write();
  WriteProfiles(px_thx);
  WriteProfiles(px_thy);
  WriteProfiles(px_xs);
  WriteProfiles(px_ys);
  WriteProfiles(px_xi);
  WriteProfiles(py_thx);
  WriteProfiles(py_thy);
  WriteProfiles(py_xs);
  WriteProfiles(py_ys);
  WriteProfiles(py_xi);
}

//----------------------------------------------------------------------------------------------------

void LinearTransportTest::WriteProfiles(const Profile<bin_type> &pr)
{
  TGraph *g_max = new TGraph();
  TGraph *g_mean = new TGraph();
  TGraph *g_sigma = new TGraph();

  const vector<bin_type> &bins = pr.GetBins();
  for (unsigned int i = 0; i < bins.size(); i++) {
    double x = pr.GetBinCenter(i);
    if (bins[i].S1 <= 0.)
      continue;
    int point = g_max->GetN();
    g_max->SetPoint(point, x, bins[i].max_abs_y);
    g_mean->SetPoint(point, x, bins[i].Mean());
    g_sigma->SetPoint(point, x, bins[i].Sigma());
  }

  gDirectory = gDirectory->mkdir(pr.name.c_str());
  g_max->Write("max");
  g_mean->Write("mean");
  g_sigma->Write("sigma");
  gDirectory->cd("..");
}

//----------------------------------------------------------------------------------------------------

void LinearTransportTest::analyze(const edm::Event &e, const edm::EventSetup &es)
{
}

//----------------------------------------------------------------------------------------------------

void LinearTransportTest::endJob()
{
}



DEFINE_FWK_MODULE(LinearTransportTest);


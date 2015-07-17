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

#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecoElasticEvent.h"
#include "TotemRawDataLibrary/DataFormats/interface/RawEvent.h"

#include <map>
#include <vector>

#include "TH1D.h"
#include "TFile.h"

#include "parameters.h"


/**
 * \ingroup ElasticReconstruction
 * \brief TODO
 **/
class RPElasticReconstructionNew : public edm::EDProducer
{
  public:
    explicit RPElasticReconstructionNew(const edm::ParameterSet& conf);
    virtual ~RPElasticReconstructionNew();

    virtual void beginRun(edm::Run&, edm::EventSetup const&) {}
    virtual void produce(edm::Event& e, const edm::EventSetup& c);
    virtual void endJob();
  
  protected:
    bool SkipBunch(unsigned int run, unsigned bunch);

    std::map<unsigned int, double> cqa, cqb, cca, ccb, ccc, csi;
    std::vector<unsigned int> cuts;	// list of active cuts

    void BuildCuts();

    void ProcessDiagonal(const RPFittedTrack &L_F, const RPFittedTrack &L_N,
      const RPFittedTrack &R_N, const RPFittedTrack &R_F, RPRecoElasticEvent &output);

  private:
      edm::InputTag fittedTrackCollectionInputLabel;
      edm::InputTag rawEventDataInputLabel;
#ifdef DEBUG
    TH1D *h_th_x, *h_th_y, *h_vtx_x;
#endif
};

//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

RPElasticReconstructionNew::RPElasticReconstructionNew(const edm::ParameterSet& conf)
#ifdef DEBUG
  : h_th_x(new TH1D("h_th_x", "", 100, 0., 0.)),
  h_th_y(new TH1D("h_th_y", "", 100, 0., 0.)),
  h_vtx_x(new TH1D("h_vtx_x", "", 100, 0., 0.))
#endif
{

  rawEventDataInputLabel = conf.getParameter<edm::InputTag>("RawEventLabel");
  fittedTrackCollectionInputLabel = conf.getParameter<edm::InputTag>("RPFittedTrackCollLabel");

  produces<RPRecoElasticEvent>();
}

//----------------------------------------------------------------------------------------------------

RPElasticReconstructionNew::~RPElasticReconstructionNew()
{
}

//----------------------------------------------------------------------------------------------------

void RPElasticReconstructionNew::endJob()
{
#ifdef DEBUG
  TFile *outF = new TFile("RPElasticReconstructionNew_debug.root", "recreate");
  h_th_x->Write();
  h_th_y->Write();
  h_vtx_x->Write();
  delete outF;
#endif
}

//----------------------------------------------------------------------------------------------------

void RPElasticReconstructionNew::produce(edm::Event& e, const edm::EventSetup& c)
{
  // get intput
  Handle< RPFittedTrackCollection > input; 
  e.getByLabel(fittedTrackCollectionInputLabel, input);

  Handle< Totem::RawEvent > rawData;
  e.getByLabel(rawEventDataInputLabel, rawData);

  // prepare output
  std::auto_ptr<RPRecoElasticEvent> output(new RPRecoElasticEvent());

  bool skip = false;

  // check trigger bits
  unsigned int input_status_bits = rawData->triggerData.input_status_bits;
  if ((input_status_bits & 3) == 0) {
    skip = true;
    output->status = RPRecoElasticEvent::sRejected;
  }
  
  // check bunch number
  unsigned int run_num = rawData->triggerData.run_num / 10000;
  unsigned int bunch_num = rawData->triggerData.bunch_num;

  if (SkipBunch(run_num, bunch_num)) {
    skip = true;
    output->status = RPRecoElasticEvent::sRejected;
  }

  if (!skip) {
    // check counter
    unsigned int count = 0;
  
    // try diagonal 45 bottom - 56 top
    RPFittedTrackCollection::const_iterator end = input->end();
    RPFittedTrackCollection::const_iterator L_F = input->find(25);
    RPFittedTrackCollection::const_iterator L_N = input->find(21);
    RPFittedTrackCollection::const_iterator R_N = input->find(120);
    RPFittedTrackCollection::const_iterator R_F = input->find(124);
  
    if (L_F != end && L_N != end && R_N != end && R_F != end) {
      Init_45b_56t();
      ProcessDiagonal(L_F->second, L_N->second, R_N->second, R_F->second, *output);
      count++;
    }
  
    // try diagonal 45 top - 56 bottom
    L_F = input->find(24);
    L_N = input->find(20);
    R_N = input->find(121);
    R_F = input->find(125);
  
    if (L_F != end && L_N != end && R_N != end && R_F != end) {
      Init_45t_56b();
      ProcessDiagonal(L_F->second, L_N->second, R_N->second, R_F->second, *output);
      count++;
    }
  
    if (count > 1)
      printf(">> RPElasticReconstructionNew::produce > WARNING: double diagonal signature.\n");
  }

  // save output
  e.put(output);
}

//----------------------------------------------------------------------------------------------------

void RPElasticReconstructionNew::ProcessDiagonal(const RPFittedTrack &L_F, const RPFittedTrack &L_N,
  const RPFittedTrack &R_N, const RPFittedTrack &R_F, RPRecoElasticEvent &output)
{
  // get hit coordinates
  double x_L_F = L_F.X0(), y_L_F = L_F.Y0();
  double x_L_N = L_N.X0(), y_L_N = L_N.Y0();
  double x_R_N = R_N.X0(), y_R_N = R_N.Y0();
  double x_R_F = R_F.X0(), y_R_F = R_F.Y0();

  // global constants
  double d_RP = 5372.;	// in mm
  double n_si = 3.;

  // vertex reconstruction
  double la_L_F = -271., la_L_N = +1059.;    // experimental L_x / v_x in mm
  double la_R_F = 221., la_R_N = -1109.;    

//  double nu_L_F = x_L_F / v_x_F;  // in mm
//  double nu_R_F = x_R_F / v_x_F;
//  double nu_L_N = x_L_N / v_x_N;
//  double nu_R_N = x_R_N / v_x_N;

  double vtx_x_L = ( x_L_F / v_x_F * la_L_N - x_L_N / v_x_N * la_L_F ) / (la_L_N - la_L_F); // in mm
  double vtx_x_R = ( x_R_F / v_x_F * la_R_N - x_R_N / v_x_N * la_R_F ) / (la_R_N - la_R_F);
  double vtx_x = (vtx_x_L + vtx_x_R) / 2.;

  // theta_y reconstruction
  double th_y_L_F = y_L_F / L_y_F / 1.00;  // in rad
  double th_y_L_N = y_L_N / L_y_N / 1.00;
  double th_y_L = - (th_y_L_N + th_y_L_F) / 2.;

  double th_y_R_F = y_R_F / L_y_F / 1.00;
  double th_y_R_N = y_R_N / L_y_N / 1.00;
  double th_y_R = + (th_y_R_N + th_y_R_F) / 2.;

  double th_y = (th_y_R + th_y_L) / 2.;

  // theta_x reconstruction
  //double th_x_R = + (x_R_N - v_x_rat * x_R_F) / L_x_N + th_x_al_R + th_x_be_R * (th_y_R - th_y_shift);

  double th_x_R = + ( (x_R_F - x_R_N) / d_RP - dvds_x * vtx_x_R ) / dLds_x * 0.999 - (th_x_al_R + th_x_be_R * th_y_R); // in rad
  double th_x_L = - ( (x_L_F - x_L_N) / d_RP - dvds_x * vtx_x_L ) / dLds_x * 1.000 - (th_x_al_L + th_x_be_L * th_y_L);
  
//  double th_x_R_noYCorr = + ( (x_R_F - x_R_N) / d_RP - dvds_x * vtx_x_R ) / dLds_x;
//  double th_x_L_noYCorr = - ( (x_L_F - x_L_N) / d_RP - dvds_x * vtx_x_L ) / dLds_x;

  double th_x = (th_x_R + th_x_L) / 2.;

  // theta reconstruction
//  double th_sq = th_x*th_x + th_y*th_y;
//  double th = sqrt(th_sq);
//  double phi = atan2(th_y, th_x);
  
  // t reconstruction
//  double t_x = p*p * th_x * th_x;
//  double t_y = p*p * th_y * th_y;
//  double t = t_x + t_y;

  // cut evaluation
  BuildCuts();
  
  cqa[1] = th_x_R;  cqb[1] = th_x_L;
  cqa[2] = th_y_R;  cqb[2] = th_y_L;
  cqa[3] = th_x_R;  cqb[3] = vtx_x_R;
  cqa[4] = th_x_L;  cqb[4] = vtx_x_L;
  cqa[5] = y_R_N;   cqb[5] = y_R_F - y_R_N;
  cqa[6] = y_L_N;   cqb[6] = y_L_F - y_L_N;
  cqa[7] = th_x;    cqb[7] = vtx_x_R - vtx_x_L;
  cqa[8] = vtx_x_R;  cqb[8] = vtx_x_L;

  bool pass = true;
  map<unsigned int, double> cv; // cut value
  map<unsigned int, bool> ct; // cut true or false
  for (map<unsigned int, double>::iterator cit = csi.begin(); cit != csi.end(); ++cit) {
    unsigned ci = cit->first;
    cv[ci] = cca[ci]*cqa[ci] + ccb[ci]*cqb[ci] + ccc[ci];
    ct[ci] = (fabs(cv[ci]) <= n_si * csi[ci]);
    pass &= ct[ci];
  }

  // fill output
  if (pass) {
    output.leftFit.th_x = th_x_L; output.leftFit.th_y = th_y_L; output.leftFit.x = vtx_x_L;
    output.rightFit.th_x = th_x_R; output.leftFit.th_y = th_y_R; output.leftFit.x = vtx_x_R;
    output.globalFit.th_x = th_x; output.leftFit.th_y = th_y; output.leftFit.x = vtx_x;
    output.result = output.globalFit;
    output.status = RPRecoElasticEvent::sOK;
  } else {
    output.status = RPRecoElasticEvent::sRejected;
  }

  // debug plots
#ifdef DEBUG
  if (pass) {
    h_th_x->Fill(th_x);
    h_th_y->Fill(th_y);
    h_vtx_x->Fill(vtx_x);
  }
#endif
}

//----------------------------------------------------------------------------------------------------

void RPElasticReconstructionNew::BuildCuts()
{
  cuts.clear();

  // cut structure:
  //  | a*qa + b*qb + c| < n_si * si

  // a: th_x_R, b: th_x_L
  cca[1] = 1.;
  ccb[1] = -1.;
  ccc[1] = cut1_c;
  csi[1] = cut1_si;
  cuts.push_back(1);
  
  // a: th_y_R, b: th_y_L
  cca[2] = 1.;
  ccb[2] = -1.;
  ccc[2] = 0.;
  csi[2] = cut2_si;
  cuts.push_back(2);
  
  // a: th_x_R, b: vtx_x_R
  cca[3] = 0.;
  ccb[3] = 1.;
  ccc[3] = 0.;
  csi[3] = cut34_si;
  cuts.push_back(3);
  
  // a: th_x_L, b: vtx_x_L
  cca[4] = 0.;
  ccb[4] = 1.;
  ccc[4] = 0.;
  csi[4] = cut34_si;
  cuts.push_back(4);
  
  // a: y_R_N, b: y_R_F - y_R_N
  cca[5] = -cut5_a;
  ccb[5] = 1.;
  ccc[5] = cut5_b;
  csi[5] = cut5_si;
  cuts.push_back(5);
  
  // a: y_L_N, b: y_L_F - y_L_N
  cca[6] = -cut6_a;
  ccb[6] = 1.;
  ccc[6] = cut6_b;
  csi[6] = cut6_si;
  cuts.push_back(6);
  
  // a: th_x, b: vtx_x_R - vtx_x_L
  cca[7] = -cut7_al;
  ccb[7] = 1.;
  ccc[7] = cut7_c;
  csi[7] = cut7_si;
  cuts.push_back(7);
}

//----------------------------------------------------------------------------------------------------

bool RPElasticReconstructionNew::SkipBunch(unsigned int run, unsigned bunch)
{
  if (run == 8318)
    if (bunch != 1885)
      return true;

  if (run >= 8333 && run <= 8341)
    if (bunch != 100)
      return true;

  if (run >= 8367 && run <= 8368)
    if (bunch != 648)
      return true;
  if (run >= 8369 && run <= 8371)
    if (bunch != 648 && bunch != 2990)
      return true;
  if (run == 8372)
    if (bunch != 648 && bunch != 2990 && bunch != 26)
      return true;

  return false;
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(RPElasticReconstructionNew);


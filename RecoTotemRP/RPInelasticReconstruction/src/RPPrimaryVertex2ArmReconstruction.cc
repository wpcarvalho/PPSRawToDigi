/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
* 	Hubert Niewiadomski
*
****************************************************************************/

#include "TotemProtonTransport/TotemRPProtonTransportParametrization/interface/LHCOpticsApproximator.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "RecoTotemRP/RPRomanPotResolutionService/interface/RPFitResolution.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProtonPair.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProtonPairCollection.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "HepMC/GenEvent.h"
#include "CLHEP/Vector/LorentzVector.h"

#include <string>
#include <map>
#include <cassert>
#include <vector>
#include <iostream>
#include <memory>

#include "TFile.h"
#include "TVector3.h"
#include "TRandom2.h"

#include "RecoTotemRP/RPInverseParameterization/interface/RPInverse2SidedParameterization.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RP2DHit.h"
#include "RecoTotemRP/RPRomanPotResolutionService/interface/RPFitResolution.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"


//----------------------------------------------------------------------------------------------------

class RPPrimaryVertex2ArmReconstruction : public edm::EDProducer
{
  public:
    typedef std::map<unsigned int, RP2DHit> rec_tracks_collection;
    
    explicit RPPrimaryVertex2ArmReconstruction(const edm::ParameterSet& conf);
    virtual ~RPPrimaryVertex2ArmReconstruction();
    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void produce(edm::Event& e, const edm::EventSetup& c);
  
  private:
    typedef std::vector<int> station_rp_ids_type;

    edm::InputTag rpFittedTrackCollectionLabel;
    edm::EDGetTokenT<edm::DetSetVector<TotemRPLocalTrack>> rpFittedTrackCollectionToken;

    edm::InputTag HepMCProductLabel;
    edm::EDGetTokenT<edm::HepMCProduct> HepMCProductToken;

    const edm::ParameterSet conf_;
    int verbosity_;
    std::auto_ptr<RPInverse2SidedParameterization> inv_param_;

    std::string param_file_name_220_right_;
    std::string param_file_name_220_left_;
    std::string param_file_name_210_right_;
    std::string param_file_name_210_left_;

    std::string param_prefix_220_right_;
    std::string param_prefix_220_left_;
    std::string param_prefix_210_right_;
    std::string param_prefix_210_left_;

    std::string right_beam_postfix_;
    std::string left_beam_postfix_;

    bool external_primary_vertex_;
    bool set_primary_vertex_to_zero_;
    bool elastic_scattering_reconstruction_;

    TVector3 primary_vertex_; 
    TVector3 primary_vertex_error_; 
    TVector3 primary_vertex_nom_pos_;

    double rp_multiple_scattering_sigma_;
    
    TRandom2 rand_;
    RPFitResolution resol_degrad_service_;

    BeamOpticsParams BOPar_;


    void InitInverseParametrizationFitter();

    int SelectHits(const edm::DetSetVector<TotemRPLocalTrack> &tracks, unsigned int armId, rec_tracks_collection &coll);

    void AddMultipleScatteringVarianceContribution(rec_tracks_collection &coll, 
      double relative_effective_length_x, double relative_effective_length_y, 
      double rp_multiple_scattering_sigma, int no_of_pots);

    void AddInStationMultipleScatteringContribution(rec_tracks_collection &col, 
          double rp_multiple_scattering_sigma);

    bool FindPrimaryVertex(edm::Event& e);

    bool CollectionContainsBotRP(const rec_tracks_collection& track_coll);
    bool CollectionContainsTopRP(const rec_tracks_collection& track_coll);   

    void Reconstruct(const rec_tracks_collection &rec_col_1, 
      const rec_tracks_collection &rec_col_2, RPInverse2SidedParameterization &inv_par, 
      RPReconstructedProtonPairCollection & rec_prot_pair_col, bool eternal_prim_vert);
    
    // TODO: remove ??
    /*

    bool station_210_mandatory_in_reconstruction_;
    bool station_220_mandatory_in_reconstruction_;
    bool any_station_in_reconstruction_;
    
    double relative_effective_length_210_220_x_right_;
    double relative_effective_length_210_220_y_right_;
    double relative_effective_length_210_220_x_left_;
    double relative_effective_length_210_220_y_left_;
    
    station_rp_ids_type station_210_right_ids_, station_210_left_ids_, 
          station_220_right_ids_, station_220_left_ids_;

    station_rp_ids_type left_arm_top_pot_ids_, left_arm_bot_pot_ids_, 
          right_arm_top_pot_ids_, right_arm_bot_pot_ids_;
    */
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

RPPrimaryVertex2ArmReconstruction::RPPrimaryVertex2ArmReconstruction(const edm::ParameterSet& conf)
 : conf_(conf), resol_degrad_service_(conf)
{
  rpFittedTrackCollectionLabel = conf.getParameter<edm::InputTag>("RPFittedTrackCollectionLabel");
  rpFittedTrackCollectionToken = consumes<DetSetVector<TotemRPLocalTrack>>(rpFittedTrackCollectionLabel);

  produces< RPReconstructedProtonPairCollection > ();
}

//----------------------------------------------------------------------------------------------------

void RPPrimaryVertex2ArmReconstruction::InitInverseParametrizationFitter()
{
  typedef std::map<unsigned int, LHCOpticsApproximator> rp_param_map_type;
  rp_param_map_type rp_param_map_right;
  rp_param_map_type rp_param_map_left;
  std::string par_name;
  
  // 220 station right
  char *cmsswPath_220 = getenv("CMSSW_BASE");
  std::string fileName_220_r = std::string(cmsswPath_220) + std::string("/src/") + param_file_name_220_right_;
  TFile *f_220_r = TFile::Open(fileName_220_r.c_str(), "read");
  par_name = param_prefix_220_right_+"_v_1_"+right_beam_postfix_;
  rp_param_map_right[120] = *((LHCOpticsApproximator*) f_220_r->Get(par_name.c_str()));
  
  rp_param_map_right[121] = *((LHCOpticsApproximator*) f_220_r->Get(par_name.c_str()));
  
  par_name = param_prefix_220_right_+"_h_1_"+right_beam_postfix_;
  rp_param_map_right[122] = *((LHCOpticsApproximator*) f_220_r->Get(par_name.c_str()));
  
  par_name = param_prefix_220_right_+"_h_2_"+right_beam_postfix_;
  rp_param_map_right[123] = *((LHCOpticsApproximator*) f_220_r->Get(par_name.c_str()));
  
  par_name = param_prefix_220_right_+"_v_2_"+right_beam_postfix_;
  rp_param_map_right[124] = *((LHCOpticsApproximator*) f_220_r->Get(par_name.c_str()));
  rp_param_map_right[125] = *((LHCOpticsApproximator*) f_220_r->Get(par_name.c_str()));
  f_220_r->Close();
  
  // 220 station left
  std::string fileName_220_l = std::string(cmsswPath_220) + std::string("/src/") + param_file_name_220_left_;
  TFile *f_220_l = TFile::Open(fileName_220_l.c_str(), "read");
  par_name = param_prefix_220_left_+"_v_1_"+left_beam_postfix_;
  rp_param_map_left[20] = *((LHCOpticsApproximator*) f_220_l->Get(par_name.c_str()));
  rp_param_map_left[21] = *((LHCOpticsApproximator*) f_220_l->Get(par_name.c_str()));
  
  par_name = param_prefix_220_left_+"_h_1_"+left_beam_postfix_;
  rp_param_map_left[22] = *((LHCOpticsApproximator*) f_220_l->Get(par_name.c_str()));
  
  par_name = param_prefix_220_left_+"_h_2_"+left_beam_postfix_;
  rp_param_map_left[23] = *((LHCOpticsApproximator*) f_220_l->Get(par_name.c_str()));
  
  par_name = param_prefix_220_left_+"_v_2_"+left_beam_postfix_;
  rp_param_map_left[24] = *((LHCOpticsApproximator*) f_220_l->Get(par_name.c_str()));
  rp_param_map_left[25] = *((LHCOpticsApproximator*) f_220_l->Get(par_name.c_str()));
  f_220_l->Close();
  
  // 210 station right
  char *cmsswPath_210 = getenv("CMSSW_BASE");
  std::string fileName_210_r = std::string(cmsswPath_210) + std::string("/src/") + param_file_name_210_right_;
  TFile *f_210_r = TFile::Open(fileName_210_r.c_str(), "read");
  par_name = param_prefix_210_right_+"_v_1_"+right_beam_postfix_;
  rp_param_map_right[100] = *((LHCOpticsApproximator*) f_210_r->Get(par_name.c_str()));
  rp_param_map_right[101] = *((LHCOpticsApproximator*) f_210_r->Get(par_name.c_str()));
  
  par_name = param_prefix_210_right_+"_h_1_"+right_beam_postfix_;
  rp_param_map_right[102] = *((LHCOpticsApproximator*) f_210_r->Get(par_name.c_str()));
  
  par_name = param_prefix_210_right_+"_h_2_"+right_beam_postfix_;
  rp_param_map_right[103] = *((LHCOpticsApproximator*) f_210_r->Get(par_name.c_str()));
  
  par_name = param_prefix_210_right_+"_v_2_"+right_beam_postfix_;
  rp_param_map_right[104] = *((LHCOpticsApproximator*) f_210_r->Get(par_name.c_str()));
  rp_param_map_right[105] = *((LHCOpticsApproximator*) f_210_r->Get(par_name.c_str()));
  f_210_r->Close();
  
  // 210 station left
  std::string fileName_210_l = std::string(cmsswPath_210) + std::string("/src/") + param_file_name_210_left_;
  TFile *f_210_l = TFile::Open(fileName_210_l.c_str(), "read");
  par_name = param_prefix_210_left_+"_v_1_"+left_beam_postfix_;
  rp_param_map_left[0] = *((LHCOpticsApproximator*) f_210_l->Get(par_name.c_str()));
  rp_param_map_left[1] = *((LHCOpticsApproximator*) f_210_l->Get(par_name.c_str()));
  
  par_name = param_prefix_210_left_+"_h_1_"+left_beam_postfix_;
  rp_param_map_left[2] = *((LHCOpticsApproximator*) f_210_l->Get(par_name.c_str()));
  
  par_name = param_prefix_210_left_+"_h_2_"+left_beam_postfix_;
  rp_param_map_left[3] = *((LHCOpticsApproximator*) f_210_l->Get(par_name.c_str()));
  
  par_name = param_prefix_210_left_+"_v_2_"+left_beam_postfix_;
  rp_param_map_left[4] = *((LHCOpticsApproximator*) f_210_l->Get(par_name.c_str()));
  rp_param_map_left[5] = *((LHCOpticsApproximator*) f_210_l->Get(par_name.c_str()));
  f_210_l->Close();
  
  inv_param_->AddParameterizationsRight(rp_param_map_right);
  inv_param_->AddParameterizationsLeft(rp_param_map_left);
  
  if(verbosity_)
  {
    for(rp_param_map_type::iterator it = rp_param_map_right.begin(); 
      it!=rp_param_map_right.end(); ++it)
    {
      std::cout<<"RPId : "<<it->first;
      it->second.PrintOpticalFunctions();
      std::cout<<std::endl;
    }
    
    for(rp_param_map_type::iterator it = rp_param_map_left.begin(); 
      it!=rp_param_map_left.end(); ++it)
    {
      std::cout<<"RPId : "<<it->first;
      it->second.PrintOpticalFunctions();
      std::cout<<std::endl;
    }
  }
}

//----------------------------------------------------------------------------------------------------

RPPrimaryVertex2ArmReconstruction::~RPPrimaryVertex2ArmReconstruction()
{
}

//----------------------------------------------------------------------------------------------------

void RPPrimaryVertex2ArmReconstruction::beginRun(edm::Run const& r, edm::EventSetup const& es)
{
  edm::ESHandle<BeamOpticsParams> BOParH;
  es.get<BeamOpticsParamsRcd>().get(BOParH);
  if(!BOParH.isValid())
    throw cms::Exception("RPPrimaryVertex2ArmReconstruction::beginRun") << " edm::ESHandle<BeamOpticsParams> is invalid";
  BOPar_ = *BOParH;
  
  inv_param_ = std::auto_ptr<RPInverse2SidedParameterization>(new RPInverse2SidedParameterization(conf_, BOPar_));

  verbosity_ = conf_.getParameter<int>("Verbosity");
  inv_param_->Verbosity(verbosity_);
  
  external_primary_vertex_ = conf_.getParameter<bool>("ExternalPrimaryVertex");
  set_primary_vertex_to_zero_ = conf_.getParameter<bool>("ConstrainPrimaryVertex");
  elastic_scattering_reconstruction_ = conf_.getParameter<bool>("ElasticScatteringReconstruction"); 
  
  external_primary_vertex_ = external_primary_vertex_ && !set_primary_vertex_to_zero_;
  
  param_file_name_210_right_ = conf_.getParameter<std::string>("ParameterizationFileName210Right");
  param_file_name_210_left_ = conf_.getParameter<std::string>("ParameterizationFileName210Left");
  param_prefix_210_right_ = conf_.getParameter<std::string>("ParameterizationNamePrefix210Right");
  param_prefix_210_left_ = conf_.getParameter<std::string>("ParameterizationNamePrefix210Left");
  
  param_file_name_220_right_ = conf_.getParameter<std::string>("ParameterizationFileName220Right");
  param_file_name_220_left_ = conf_.getParameter<std::string>("ParameterizationFileName220Left");
  param_prefix_220_right_ = conf_.getParameter<std::string>("ParameterizationNamePrefix220Right");
  param_prefix_220_left_ = conf_.getParameter<std::string>("ParameterizationNamePrefix220Left");
  
  if(external_primary_vertex_)
  {
    HepMCProductLabel = conf_.getParameter<edm::InputTag>("HepMCProductLabel");
    HepMCProductToken = consumes<edm::HepMCProduct>(HepMCProductLabel);
    primary_vertex_error_.SetX(conf_.getParameter<double>("PrimaryVertexXSigma"));
    primary_vertex_error_.SetY(conf_.getParameter<double>("PrimaryVertexYSigma"));
    primary_vertex_error_.SetZ(conf_.getParameter<double>("PrimaryVertexZSigma"));
  }
  
  if(set_primary_vertex_to_zero_)
  {
    // conversion from meters to mm
    primary_vertex_error_.SetX(BOPar_.GetPrimVertSizeX()*1000.0);
    primary_vertex_error_.SetY(BOPar_.GetPrimVertSizeY()*1000.0);
    primary_vertex_error_.SetZ(BOPar_.GetPrimVertSizeZ()*1000.0);
    primary_vertex_nom_pos_.SetX(BOPar_.GetBeamDisplacementX()*1000.0);
    primary_vertex_nom_pos_.SetY(BOPar_.GetBeamDisplacementY()*1000.0);
    primary_vertex_nom_pos_.SetZ(BOPar_.GetBeamDisplacementZ()*1000.0);
  }
  
  right_beam_postfix_ = conf_.getParameter<std::string>("RightBeamPostfix");
  left_beam_postfix_ = conf_.getParameter<std::string>("LeftBeamPostfix");
  
  rp_multiple_scattering_sigma_ = conf_.getParameter<double>("RPMultipleScatteringSigma");
  
  InitInverseParametrizationFitter();
}

//----------------------------------------------------------------------------------------------------

void RPPrimaryVertex2ArmReconstruction::produce(edm::Event& e, const edm::EventSetup& c)
{
  // get input
  edm::Handle< DetSetVector<TotemRPLocalTrack> > input; 
  RPReconstructedProtonPairCollection reconstructed_proton_pair_collection;
  e.getByToken(rpFittedTrackCollectionToken, input);

  // process primary vertex
  if (external_primary_vertex_)
  {
    if(!FindPrimaryVertex(e))
      throw cms::Exception("RPPrimaryVertex2ArmReconstruction::produce") << "Primary vertex not found." << endl;
  }
  
  if (set_primary_vertex_to_zero_)
  {
    primary_vertex_ = primary_vertex_nom_pos_;
  }

  // split input hits per station
  rec_tracks_collection hits_l;
  SelectHits(*input, 0, hits_l);

  rec_tracks_collection hits_r;
  SelectHits(*input, 1, hits_r);
  
  // can proton in an arm be reconstructed?
  bool left_reconstructable = (hits_l.size() >= 2);
  bool right_reconstructable = (hits_r.size() >= 2);
  
  // check for symmetric location of vertical roman pots in elastic reconstruction
  bool allow_elastic_recon = false;
  if(elastic_scattering_reconstruction_)
  {
    allow_elastic_recon = allow_elastic_recon || (CollectionContainsBotRP(hits_r) && CollectionContainsTopRP(hits_l));
    allow_elastic_recon = allow_elastic_recon || (CollectionContainsBotRP(hits_l) && CollectionContainsTopRP(hits_r));
  } else {
    allow_elastic_recon = true;
  }
  
  if (right_reconstructable && left_reconstructable && allow_elastic_recon)
  {
    AddInStationMultipleScatteringContribution(hits_l, rp_multiple_scattering_sigma_);
    AddInStationMultipleScatteringContribution(hits_r, rp_multiple_scattering_sigma_);
    
    Reconstruct(hits_l, hits_r, *inv_param_, 
        reconstructed_proton_pair_collection, 
        external_primary_vertex_ || set_primary_vertex_to_zero_);
  }
  
  e.put(make_unique<RPReconstructedProtonPairCollection>(reconstructed_proton_pair_collection));
}

//----------------------------------------------------------------------------------------------------

bool RPPrimaryVertex2ArmReconstruction::CollectionContainsTopRP(const rec_tracks_collection& track_coll)
{
  rec_tracks_collection::const_iterator it;
  rec_tracks_collection::const_iterator beg = track_coll.begin();
  rec_tracks_collection::const_iterator end = track_coll.end();
  
  for(it = beg; it!=end; ++it)
  {
    int rp_pos = it->first % 10;
    if(rp_pos == 0 || rp_pos == 4)
      return true;
  }
  return false;
}

//----------------------------------------------------------------------------------------------------

bool RPPrimaryVertex2ArmReconstruction::CollectionContainsBotRP(const rec_tracks_collection& track_coll)
{
  rec_tracks_collection::const_iterator it;
  rec_tracks_collection::const_iterator beg = track_coll.begin();
  rec_tracks_collection::const_iterator end = track_coll.end();
  
  for(it = beg; it!=end; ++it)
  {
    int rp_pos = it->first % 10;
    if(rp_pos == 1 || rp_pos == 5)
      return true;
  }
  return false;
}

//----------------------------------------------------------------------------------------------------

void RPPrimaryVertex2ArmReconstruction::AddInStationMultipleScatteringContribution(rec_tracks_collection &col, 
    double rp_multiple_scattering_sigma)
{
  if (verbosity_)
  {
    rec_tracks_collection::iterator it = col.begin();
    std::cout<<"RPPrimaryVertexInelasticReconstruction::AddInStationMultipleScatteringContribution"<<std::endl;
    for(; it!=col.end(); ++it)
    {
      std::cout<<"rp:"<<it->first<<" sigma=("<<it->second.Sx()<<","<<it->second.Sy()<<")"<<std::endl;
      std::cout<<"\tpos=("<<it->second.X()<<","<<it->second.Y()<<","<<it->second.Z()<<")"<<std::endl;
    }
  }

  double sig2 = rp_multiple_scattering_sigma*rp_multiple_scattering_sigma;

  // traverse all hit pairs
  for (auto &h_dest : col)
  {
    //printf("RP %u: multiple scattering contributions from ", h_dest.first);

    double scat_var_contrib = 0.0;

    for (const auto &h_src : col)
    {
      // skip combinations of twice the same hit
      if (h_src.first == h_dest.first)
        continue;

      // calculate distance from src to dest, with Z increasing from the IP
      double dist = fabs(h_dest.second.Z()) - fabs(h_src.second.Z());

      // ensure causality: dist > 0
      if (dist < 0.)
        continue;

      //printf("%u, ", h_src.first);

      scat_var_contrib += sig2 * dist*dist;
    }

    double vx = h_dest.second.Vx() + scat_var_contrib;
    double vy = h_dest.second.Vy() + scat_var_contrib;

    //printf("\n    multiple_scat_position_sigma=%.3E, total sigma x=%.3E, y=%.3E\n", sqrt(scat_var_contrib), sqrt(vx), sqrt(vy));

    h_dest.second.Vx(vx);
    h_dest.second.Vy(vy);
  }
}

//----------------------------------------------------------------------------------------------------

bool RPPrimaryVertex2ArmReconstruction::FindPrimaryVertex(edm::Event& e)
{
  edm::Handle<edm::HepMCProduct> HepMCEvt;
  e.getByToken(HepMCProductToken, HepMCEvt);

  if(!HepMCEvt.isValid())
  {
    throw cms::Exception("RPPrimaryVertex2ArmReconstruction::FindPrimaryVertex") <<
      "Unable to find HepMCProduct(HepMC::GenEvent) in edm::Event" << endl;
  }
  
  const HepMC::GenEvent *evt = HepMCEvt->GetEvent();
  
  int vertex_number = 0;
  for(HepMC::GenEvent::vertex_const_iterator vitr = evt->vertices_begin(); 
        vitr != evt->vertices_end(); ++vitr )
  {
    ++vertex_number;
    const HepMC::FourVector &xvtx = (*vitr)->position();
    primary_vertex_.SetX(xvtx.x() + rand_.Gaus(0, primary_vertex_error_.X()) );
    primary_vertex_.SetY(xvtx.y() + rand_.Gaus(0, primary_vertex_error_.Y()) );
    primary_vertex_.SetZ(xvtx.z() + rand_.Gaus(0, primary_vertex_error_.Z()) );
    if(verbosity_)
    {
      std::cout<<"Vertex found: [" << xvtx.x() << ", " << xvtx.y() << ", " << xvtx.z() << ", " << xvtx.t()
    		<< "]  smeared:("<<primary_vertex_.X()<<", " <<primary_vertex_.Y()<<", "<<
    		primary_vertex_.Z()<<")"<<std::endl;
    }
  }

  return vertex_number==1;
}

//----------------------------------------------------------------------------------------------------

// return the number of pots through which the proton goes
int RPPrimaryVertex2ArmReconstruction::SelectHits(const DetSetVector<TotemRPLocalTrack> &tracks, 
    unsigned int armId_request, rec_tracks_collection &coll)
{
  // filter hits with selected arm id
  bool top = false;
  bool bottom = false;
  set<unsigned int> units;
  for (const auto &ds : tracks)
  {
    const unsigned int &rpId = ds.detId();
    unsigned int armId = rpId / 100;
    if (armId != armId_request)
      continue;

    // currently support only 1 track per RP
    if (ds.size() != 1)
      throw cms::Exception("RPPrimaryVertex2ArmReconstruction::SelectHits") << ds.size() << "tracks found in RP " << rpId << endl;

    const TotemRPLocalTrack &tr = ds[0];

    if (!tr.isValid())
      continue;

    coll[rpId] = resol_degrad_service_.Create2DHit(rpId, tr);

    // for quality checks
    unsigned int rpNum = rpId % 10;
    if (rpNum == 0 || rpNum == 4)
      top = true;
    if (rpNum == 1 || rpNum == 5)
      bottom = true;

    unsigned int unitId = (rpId / 10) * 10;
    if (rpNum > 2)
      unitId++;
    units.insert(unitId);
  }

  //printf("\n");

  // quality check
  bool collection_accepted = (units.size() >= 2) && ( !(top && bottom) );
  
  if (!collection_accepted)
    coll.clear();
  
  return coll.size();
}

//----------------------------------------------------------------------------------------------------

void RPPrimaryVertex2ArmReconstruction::Reconstruct(const rec_tracks_collection &rec_col_1, 
    const rec_tracks_collection &rec_col_2, RPInverse2SidedParameterization &inv_par, 
    RPReconstructedProtonPairCollection & rec_prot_pair_col, bool external_prim_vertex)
{
  inv_par.ClearEvent();
  inv_par.AddProtonAtRPCollection(rec_col_1);
  inv_par.AddProtonAtRPCollection(rec_col_2);
  
  if(external_prim_vertex)
    inv_par.SetPrimaryVertex(primary_vertex_, primary_vertex_error_);
  
  RPReconstructedProtonPair rec_prot_pair;
  if(verbosity_)
    inv_par.PrintFittedHitsInfo(std::cout);
    
  inv_par.Fit(rec_prot_pair);
  rec_prot_pair_col.push_back(rec_prot_pair);
}

DEFINE_FWK_MODULE(RPPrimaryVertex2ArmReconstruction);

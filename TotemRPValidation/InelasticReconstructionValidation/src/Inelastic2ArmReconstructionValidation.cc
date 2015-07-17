#include "TotemRPValidation/InelasticReconstructionValidation/interface/Inelastic2ArmReconstructionValidation.h"
#include "TFile.h"
#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
//#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "SimG4Core/Notification/interface/SimG4Exception.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <boost/shared_ptr.hpp>
#include <string>
#include "TRandom2.h"
#include "SimG4CMS/TotemRPProtTranspPar/interface/LHCOpticsApproximator.h"
#include "TNamed.h"
#include "TObject.h"
#include "TROOT.h"
#include "TClassTable.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProtonPair.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProtonPairCollection.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include <memory>

Inelastic2ArmReconstructionValidation::Inelastic2ArmReconstructionValidation(const edm::ParameterSet& conf)
:  conf_(conf)
{
  //test_profile_.SetDirectory(0);
  hist_file_name_ = conf.getParameter<std::string>("HistogramFileName");
  HepMCProductLabelVertex_ = conf.getParameter<std::string>("SmearedVertexHepMCProductLabel");
  HepMCProductLabelProton_ = conf.getParameter<std::string>("PrimaryProtonHepMCProductLabel");
  verbosity_ = conf.getParameter<int>("Verbosity");
  rand_gen_.SetSeed(conf.getParameter<int>("RandomGenSeed"));
  HepMCVertexModuleName_ = conf.getParameter<std::string>("SmearedVertexHepMCModuleName");
  HepMCProtonModuleName_ = conf.getParameter<std::string>("PrimaryProtonHepMCModuleName");
  SDValidation_ = conf.getParameter<bool>("SDValidation");
  rpReconstructedProtonPairCollectionLabel = conf.getParameter<edm::InputTag>("RPReconstructedProtonPairCollectionLabel");

}


Inelastic2ArmReconstructionValidation::~Inelastic2ArmReconstructionValidation()
{
}


void Inelastic2ArmReconstructionValidation::beginRun(edm::Run const&, edm::EventSetup const& es)
{
  edm::ESHandle<BeamOpticsParams> BOParH;
  es.get<BeamOpticsParamsRcd>().get(BOParH);
  if(!BOParH.isValid())
    throw cms::Exception("Inelastic2ArmReconstructionValidation") << " edm::ESHandle<BeamOpticsParams> is invalid";
  
  BOPar_ = *BOParH;
  mad_param_transport_ = std::auto_ptr<ParamMADRefTransport>(new ParamMADRefTransport(conf_, es));
  prot_pair_histograms_ = std::auto_ptr<RPProtonPairReconstructionInfo>(new RPProtonPairReconstructionInfo("/prot_pair_rec_val/", conf_));
  residual_histograms_ = std::auto_ptr<BaseCollectionManager<PoolsReconstructionInfo, int, edm::ParameterSet> >(new BaseCollectionManager<PoolsReconstructionInfo, int, edm::ParameterSet>("/residuals/RP_", conf_));
}


void Inelastic2ArmReconstructionValidation::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  if(!FindPrimaryProtons(e) )
  {
    if(verbosity_)
      std::cout<<"Inelastic2ArmReconstructionValidation : 2 Primary protons not found!! Skipping the event."<<std::endl;
    return;
  }
  if(!FindPrimaryVertex(e) )
  {
    if(verbosity_)
      std::cout<<"Inelastic2ArmReconstructionValidation : Primary vertex not found!! Skipping the event."<<std::endl;
    return;
  }
  bool rec_pair_found;
  if(!(rec_pair_found=FindReconstrucedProtons(e)) )
  {
    if(verbosity_)
      std::cout<<"Inelastic2ArmReconstructionValidation : Reconstructed protons not found!! Skipping the event."<<std::endl;
  }
  
  if(verbosity_)
    std::cout<<"Filling the reference histograms."<<std::endl;
  
  FillReferenceHistograms();
  if(rec_pair_found)
  {
    FillResidualHistograms();
  }
}


bool Inelastic2ArmReconstructionValidation::FindPrimaryProtons(const edm::Event& e)
{
  edm::Handle<edm::HepMCProduct> HepMCEvt;
  e.getByLabel(HepMCProtonModuleName_, HepMCProductLabelProton_, HepMCEvt ) ;
  
  primary_protons_.clear();
  if(!HepMCEvt.isValid())
  {
    throw SimG4Exception("Inelastic2ArmReconstructionValidation::FindPrimaryProtons : Unable to find HepMCProduct(HepMC::GenEvent) in edm::Event  ");
  }
  
  const HepMC::GenEvent *evt = HepMCEvt->GetEvent();

  int right_count = 0;
  int left_count = 0;
  
  for(HepMC::GenEvent::particle_const_iterator it = evt->particles_begin(); 
      it != evt->particles_end(); ++it )
  {

    HepMC::GenParticle * g = (*it);
    int g_status = g->status();
    
    int pdg_id = g->pdg_id();
  
    if(verbosity_)
    {
      g->print(std::cout);
      if (g->production_vertex()){
	    const HepMC::FourVector &v = g->production_vertex()->position();
        std::cout<< "[" << v.x() << ", " << v.y() << ", " << v.z() << ", " << v.t() << "]"<<std::endl;
      }
    
    }
    
    // scanning only for particles with status == 1 
    if (g_status == 1 && pdg_id == 2212)
    {
      if(verbosity_)
        std::cout<<"Setting the primary protons."<<std::endl;
      
      if(g->momentum().z()>0 && (!right_count || g->momentum().rho() > primary_protons_[1]->momentum().rho()))
      {
        ++right_count;
        primary_protons_[1]=g;
      }
      if(g->momentum().z()<0 && (!left_count ||  g->momentum().rho() > primary_protons_[0]->momentum().rho()))
      {
        ++left_count;
        primary_protons_[0]=g;
      }
    }
  } // end loop on HepMC particles
  
  if(SDValidation_ && right_count && left_count)
  {
    if(primary_protons_[1]->momentum().rho() > primary_protons_[0]->momentum().rho())
    {
      left_count = 0;
    }
    else
    {
      right_count = 0;
    }
  }
  bool result = right_count>0 || left_count>0;
  
//  bool result = right_count<=1 && left_count<=1 && (right_count>0 || left_count>0);

  if(verbosity_)
    std::cout<<"right_count="<<right_count<<" left_count="<<left_count
        <<" result="<<result<<std::endl;
  return result;
}

bool Inelastic2ArmReconstructionValidation::FindPrimaryVertex(const edm::Event& e)
{
  edm::Handle<edm::HepMCProduct> HepMCEvt;
  e.getByLabel(HepMCVertexModuleName_, HepMCProductLabelVertex_, HepMCEvt ) ;
  
  primary_vertex_.clear();
  if(!HepMCEvt.isValid())
  {
    throw SimG4Exception("Inelastic2ArmReconstructionValidation::FindPrimaryVertex : Unable to find HepMCProduct(HepMC::GenEvent) in edm::Event  ");
  }
  
  const HepMC::GenEvent *evt = HepMCEvt->GetEvent();
  
  int right_count = 0;
  int left_count = 0;
  
  for(HepMC::GenEvent::particle_const_iterator it = evt->particles_begin(); 
      it != evt->particles_end(); ++it )
  {
    HepMC::GenParticle * g = (*it);
    int g_status = g->status();
    int pdg_id = g->pdg_id();
    
    if(verbosity_)
    {
      g->print(std::cout);
      if (g->production_vertex()){
  	    const HepMC::FourVector &v = g->production_vertex()->position();
        std::cout<< "[" << v.x() << ", " << v.y() << ", " << v.z() << ", " << v.t() << "]"<<std::endl;
      }
    }
    
    // scanning only for particles with status == 1 
    if (g_status == 1 && pdg_id == 2212)
    {
      if(verbosity_)
        std::cout<<"Setting the vertex."<<std::endl;
      if(g->momentum().z()>0 && (!right_count || g->momentum().rho() > primary_protons_[1]->momentum().rho()))
      {
        ++right_count;
        primary_vertex_[1]=g;
      }
      if(g->momentum().z()<0 && (!left_count ||  g->momentum().rho() > primary_protons_[0]->momentum().rho()))
      {
        ++left_count;
        primary_vertex_[0]=g;
      }
    }
  } // end loop on HepMC particles
  
  if(SDValidation_ && right_count && left_count)
  {
    if(primary_protons_[1]->momentum().rho() > primary_protons_[0]->momentum().rho())
    {
      left_count = 0;
    }
    else
    {
      right_count = 0;
    }
  }
  
  bool result = right_count>0 || left_count>0;
  
  if(verbosity_)
    std::cout<<"right_count="<<right_count<<" left_count="<<left_count
        <<" result="<<result<<std::endl;
  return result;
}


bool Inelastic2ArmReconstructionValidation::FindReconstrucedProtons(const edm::Event& e)
{
  if(verbosity_)
    std::cout<<"Finding the reconstructed protons"<<std::endl;
  edm::Handle< RPReconstructedProtonPairCollection > input;
  e.getByLabel(rpReconstructedProtonPairCollectionLabel, input);
  
  if(!input.isValid() || input->size()!=1)
  {
    if(verbosity_)
      std::cout<<"reconstructed protons not found"<<std::endl;
    return false;
  }
  else
  {
    if(verbosity_)
      std::cout<<"reconstructed protons found"<<std::endl;
    reconstructed_protons_ = (*input)[0];
    return true;
  }
}


void Inelastic2ArmReconstructionValidation::FillReferenceHistograms()
{
  if(primary_protons_.find(0)==primary_protons_.end() 
      || primary_protons_.find(1)==primary_protons_.end()
      || primary_vertex_.find(0)==primary_vertex_.end()
      || primary_vertex_.find(1)==primary_vertex_.end())
  {
    if(verbosity_)
      std::cout<<"Primary vertex or primary protons not found"<<std::endl;
    return;
  }
  
  prot_pair_histograms_->FillReferenceHistograms(BOPar_, primary_protons_[0], 
      primary_protons_[1], primary_vertex_[0], primary_vertex_[1], verbosity_);
}


void Inelastic2ArmReconstructionValidation::FillResidualHistograms()
{
  if(primary_protons_.find(0)==primary_protons_.end() 
      || primary_protons_.find(1)==primary_protons_.end()
      || primary_vertex_.find(0)==primary_vertex_.end()
      || primary_vertex_.find(1)==primary_vertex_.end())
    return;
  
  prot_pair_histograms_->FillResidualHistograms(BOPar_, primary_protons_[0], 
      primary_protons_[1], primary_vertex_[0], primary_vertex_[1], 
      &reconstructed_protons_, verbosity_);
  FillRPPools();
}


void Inelastic2ArmReconstructionValidation::FillRPPools()
{
  if(primary_protons_.find(0)==primary_protons_.end() 
      || primary_protons_.find(1)==primary_protons_.end() 
      || primary_vertex_.find(0)==primary_vertex_.end()
      || primary_vertex_.find(1)==primary_vertex_.end()
      || !reconstructed_protons_.Valid()
      || reconstructed_protons_.KsiLeft()>0.01 || reconstructed_protons_.KsiRight()>0.01)
    return;
  
  //left
  double pr_px0 = primary_protons_[0]->momentum().x();
  double pr_py0 = primary_protons_[0]->momentum().y();
  double pr_pz0 = primary_protons_[0]->momentum().z();
  
  double pr_vx0 = primary_vertex_[0]->production_vertex()->position().x();//+BOPar_.GetBeamDisplacementX()*1000.0;
  double pr_vy0 = primary_vertex_[0]->production_vertex()->position().y();//+BOPar_.GetBeamDisplacementY()*1000.0;
  double pr_vz0 = primary_vertex_[0]->production_vertex()->position().z();//+BOPar_.GetBeamDisplacementZ()*1000.0;
  
  //right
  double pr_px1 = primary_protons_[1]->momentum().x();
  double pr_py1 = primary_protons_[1]->momentum().y();
  double pr_pz1 = primary_protons_[1]->momentum().z();
  
  double pr_vx1 = primary_vertex_[1]->production_vertex()->position().x();//+BOPar_.GetBeamDisplacementX()*1000.0;
  double pr_vy1 = primary_vertex_[1]->production_vertex()->position().y();//+BOPar_.GetBeamDisplacementY()*1000.0;
  double pr_vz1 = primary_vertex_[1]->production_vertex()->position().z();//+BOPar_.GetBeamDisplacementZ()*1000.0;
  
  //left
  RPRecoProtMADXVariables madx_left_var = reconstructed_protons_.GetMADXVariablesLeft();
  HepMC::FourVector rec_p0 = BOPar_.LeftProtonMADXCanonicalVariablesToP(madx_left_var); 
  //  double rec_px0 = reconstructed_protons_.PxLeft();
  //  double rec_py0 = reconstructed_protons_.PyLeft();
  //  double rec_pz0 = reconstructed_protons_.PzLeft();
  //right
  RPRecoProtMADXVariables madx_right_var = reconstructed_protons_.GetMADXVariablesRight();
  HepMC::FourVector rec_p1 = BOPar_.RightProtonMADXCanonicalVariablesToP(madx_right_var);
  //  double rec_px1 = reconstructed_protons_.PxRight();
  //  double rec_py1 = reconstructed_protons_.PyRight();
  //  double rec_pz1 = reconstructed_protons_.PzRight();
  
  double rec_vx = reconstructed_protons_.X3D();
  double rec_vy = reconstructed_protons_.Y3D();
  double rec_vz = reconstructed_protons_.Z3D();
  
  TVector2 prim_pos;
  TVector2 rec_pos;
  
  for(RPReconstructedProton::debug_hits_map_type::const_iterator 
      it = reconstructed_protons_.DebugHits().begin();
  it!=reconstructed_protons_.DebugHits().end(); ++it)
  {
    bool prim_res;
    bool rec_res;
    
    if(LeftArm(it->first))
    {
      prim_res = mad_param_transport_->MADXTeoreticalRPTransversePosition(it->first, 
          pr_vx0, pr_vy0, pr_vz0, pr_px0, pr_py0, pr_pz0, prim_pos);
      rec_res = mad_param_transport_->MADXTeoreticalRPTransversePosition(it->first, 
          rec_vx, rec_vy, rec_vz, rec_p0.px(), rec_p0.py(), rec_p0.pz(), rec_pos);
    }
    else
    {
      prim_res = mad_param_transport_->MADXTeoreticalRPTransversePosition(it->first, 
          pr_vx1, pr_vy1, pr_vz1, pr_px1, pr_py1, pr_pz1, prim_pos);
      rec_res = mad_param_transport_->MADXTeoreticalRPTransversePosition(it->first, 
          rec_vx, rec_vy, rec_vz, rec_p1.px(), rec_p1.py(), rec_p1.pz(), rec_pos);
    }
    
    if(prim_res && rec_res)
    {
      residual_histograms_->GetObj(it->first)->Fill(prim_pos, it->second, rec_pos);
    }
  }
}


void Inelastic2ArmReconstructionValidation::endJob()
{
  WriteHistograms(hist_file_name_);
}


void Inelastic2ArmReconstructionValidation::WriteHistograms(const std::string &root_file_name)
{
  TFile *f = TFile::Open(root_file_name.c_str(), "recreate");
  if(!f || !f->IsWritable())
  {
    std::cout<<"Output file not opened correctly!!"<<std::endl;
  }
  
  if(verbosity_)
    std::cout<<"Writting histograms..."<<std::endl;
  prot_pair_histograms_->WriteRootFile(f);
  residual_histograms_->Write(f);
  if(verbosity_)
    std::cout<<"Writting histograms finnished."<<std::endl;
  f->Close();
}

DEFINE_FWK_MODULE(Inelastic2ArmReconstructionValidation);

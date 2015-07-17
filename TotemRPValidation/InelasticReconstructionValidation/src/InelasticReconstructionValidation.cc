#include "TotemRPValidation/InelasticReconstructionValidation/interface/InelasticReconstructionValidation.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProtonCollection.h"
#include "TFile.h"
#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
//#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "SimG4Core/Notification/interface/SimG4Exception.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProton.h"
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
#include "FWCore/Framework/interface/MakerMacros.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include <memory>


InelasticReconstructionValidation::InelasticReconstructionValidation(const edm::ParameterSet& conf)
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
  rpReconstructedProtonCollectionLabel = conf.getParameter<edm::InputTag>("RPReconstructedProtonCollectionLabel");

}


InelasticReconstructionValidation::~InelasticReconstructionValidation()
{
}


void InelasticReconstructionValidation::beginRun(edm::Run const&, edm::EventSetup const& es)
{
  edm::ESHandle<BeamOpticsParams> BOParH;
  es.get<BeamOpticsParamsRcd>().get(BOParH);
  if(!BOParH.isValid())
    throw cms::Exception("InelasticReconstructionValidation") << " edm::ESHandle<BeamOpticsParams> is invalid";
  
  BOPar_ = *BOParH;
  mad_param_transport_ = std::auto_ptr<ParamMADRefTransport>(new ParamMADRefTransport(conf_, es));
  prot_histograms_ = std::auto_ptr<BaseCollectionManager<RPProtonReconstructionInfo, int, edm::ParameterSet> >(new BaseCollectionManager<RPProtonReconstructionInfo, int, edm::ParameterSet>("/prot_rec_val/Arm_", conf_));
  residual_histograms_ = std::auto_ptr<BaseCollectionManager<PoolsReconstructionInfo, int, edm::ParameterSet> >(new BaseCollectionManager<PoolsReconstructionInfo, int, edm::ParameterSet>("/residuals/RP_", conf_));
}


void InelasticReconstructionValidation::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  if(!FindPrimaryProtons(e) )
  {
    if(verbosity_)
      std::cout<<"InelasticReconstructionValidation : Primary protons not found!! Skipping the event."<<std::endl;
    return;
  }
  if(!FindPrimaryVertex(e) )
  {
    if(verbosity_)
      std::cout<<"InelasticReconstructionValidation : Primary vertex not found!! Skipping the event."<<std::endl;
    return;
  }
  if(!FindReconstrucedProtons(e) )
  {
    if(verbosity_)
      std::cout<<"InelasticReconstructionValidation : Reconstructed protons not found!! Skipping the event."<<std::endl;
  }
  
  for(primary_prot_map_type::const_iterator pr_it = primary_protons_.begin(); 
        pr_it!=primary_protons_.end(); ++pr_it )
  {
    primary_prot_map_type::const_iterator pr_vert_it = primary_vertex_.find(pr_it->first);
    if(pr_vert_it==primary_vertex_.end())
    {
      std::cout<<"Primary vertex arm_id="<<pr_it->first<<" not found! Skipping..."<<std::endl;
      continue;
    }
      
    FillReferenceHistograms(pr_it->first, pr_it->second, pr_vert_it->second);
    reconstructed_prot_map_type::const_iterator rec_it = reconstructed_protons_.find(pr_it->first);
    if(rec_it!=reconstructed_protons_.end())
    {
      FillResidualHistograms(pr_it->first, pr_it->second, pr_vert_it->second, rec_it->second);
    }
  }
}


bool InelasticReconstructionValidation::FindPrimaryProtons(const edm::Event& e)
{
  edm::Handle<edm::HepMCProduct> HepMCEvt;
  e.getByLabel(HepMCProtonModuleName_, HepMCProductLabelProton_, HepMCEvt ) ; //change the way you get it!!!!!! 

  primary_protons_.clear();

  if(!HepMCEvt.isValid())
  {
     throw SimG4Exception("InelasticReconstructionValidation : Unable to find HepMCProduct(HepMC::GenEvent) in edm::Event  ");
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
    
//    const HepMC::FourVector &vtx = g->production_vertex()->position();
//    const HepMC::FourVector &mom  = g->momentum();
    
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
      primary_protons_.erase(0);
    }
    else
    {
      right_count = 0;
      primary_protons_.erase(1);
    }
  }
  bool result = right_count>0 || left_count>0;
  
//  bool result = right_count<=1 && left_count<=1 && (right_count>0 || left_count>0);

  if(verbosity_)
    std::cout<<"right_count="<<right_count<<" left_count="<<left_count
        <<" result="<<result<<std::endl;
  return result;
}


bool InelasticReconstructionValidation::FindPrimaryVertex(const edm::Event& e)
{
  edm::Handle<edm::HepMCProduct> HepMCEvt;
  e.getByLabel(HepMCVertexModuleName_, HepMCProductLabelVertex_, HepMCEvt ) ;
  
  primary_vertex_.clear();
  if(!HepMCEvt.isValid())
  {
     throw SimG4Exception("InelasticReconstructionValidation : Unable to find HepMCProduct(HepMC::GenEvent) in edm::Event  ");
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
      if(g->momentum().z()>0 && (!right_count || g->momentum().rho() > primary_vertex_[1]->momentum().rho()))
      {
        ++right_count;
        primary_vertex_[1]=g;
      }
      if(g->momentum().z()<0 && (!left_count ||  g->momentum().rho() > primary_vertex_[0]->momentum().rho()))
      {
        ++left_count;
        primary_vertex_[0]=g;
      }
    }
  } // end loop on HepMC particles
  
  if(SDValidation_ && right_count && left_count)
  {
    if(primary_vertex_[1]->momentum().rho() > primary_vertex_[0]->momentum().rho())
    {
      left_count = 0;
      primary_vertex_.erase(0);
    }
    else
    {
      right_count = 0;
      primary_vertex_.erase(1);
    }
  }
  
  bool result = right_count>0 || left_count>0;
  
  if(verbosity_)
    std::cout<<"right_count="<<right_count<<" left_count="<<left_count
        <<" result="<<result<<std::endl;
  return result;
}


bool InelasticReconstructionValidation::FindReconstrucedProtons(const edm::Event& e)
{
  reconstructed_protons_.clear();
  
  if(verbosity_)
    std::cout<<"Finding the reconstructed protons"<<std::endl;
  edm::Handle< RPReconstructedProtonCollection > input;
  e.getByLabel(rpReconstructedProtonCollectionLabel, input);
  
  int right = 0;
  int left = 0;
  for(RPReconstructedProtonCollection::const_iterator it = input->begin(); 
        it!=input->end(); ++it)
  {
    if(verbosity_)
      std::cout<<"Reconstructed proton ksi"<< it->Ksi()
            <<"zdirection:"<<it->ZDirection()<<std::endl;
    if(it->ZDirection() > 0)
    {
      ++right;
      reconstructed_protons_[1] = &(*it);
      if(verbosity_)
        std::cout<<"Setting the reconstructed proton right"<<std::endl;
    }
    if(it->ZDirection() < 0)
    {
      ++left;
      reconstructed_protons_[0] = &(*it);
      if(verbosity_)
        std::cout<<"Setting the reconstructed proton left"<<std::endl;
    }
  }
  bool result = right<=1 && left<=1 && (right>0 || left>0);
  return result;
}


void InelasticReconstructionValidation::FillReferenceHistograms(int arm_id, 
      const HepMC::GenParticle* pr_prot, const HepMC::GenParticle* pr_vert)
{
  prot_histograms_->GetObj(arm_id)->FillReferenceHistograms(BOPar_, pr_prot, pr_vert, verbosity_);
}


void InelasticReconstructionValidation::FillResidualHistograms(int arm_id, 
      const HepMC::GenParticle* pr_prot, const HepMC::GenParticle* pr_vert,
      const RPReconstructedProton* rec_prot)
{
  prot_histograms_->GetObj(arm_id)->FillResidualHistograms(BOPar_, pr_prot, pr_vert, rec_prot, verbosity_);
  FillRPPools(pr_prot, pr_vert, rec_prot);
}


void InelasticReconstructionValidation::FillRPPools(
      const HepMC::GenParticle* pr_prot, const HepMC::GenParticle* pr_vert, 
      const RPReconstructedProton* rec_prot)
{
  if(!pr_prot || !rec_prot || !rec_prot->Valid() || !(rec_prot->Ksi()<0.01))
    return;
  
  double pr_px = pr_prot->momentum().x();
  double pr_py = pr_prot->momentum().y();
  double pr_pz = pr_prot->momentum().z();
  
  double pr_vx = pr_vert->production_vertex()->position().x();//+BOPar_.GetBeamDisplacementX()*1000.0;
  double pr_vy = pr_vert->production_vertex()->position().y();//+BOPar_.GetBeamDisplacementY()*1000.0;
  double pr_vz = pr_vert->production_vertex()->position().z();//+BOPar_.GetBeamDisplacementZ()*1000.0;
  
  //project the primary vertex on z=0 along the primary proton momentum direction
  pr_vx -= pr_px/pr_pz*pr_vz;
  pr_vy -= pr_py/pr_pz*pr_vz;
  
  RPRecoProtMADXVariables madx_var = rec_prot->GetMADXVariables();
  HepMC::FourVector rec_p = BOPar_.ProtonMADXCanonicalVariablesToP(madx_var, rec_prot->ZDirection()); 
//  double rec_px = rec_prot->Px();
//  double rec_py = rec_prot->Py();
//  double rec_pz = rec_prot->Pz();
  
  double rec_vx = rec_prot->X();
  double rec_vy = rec_prot->Y();
  
  TVector2 prim_pos;
  TVector2 rec_pos;
  
  for(RPReconstructedProton::debug_hits_map_type::const_iterator 
        it = rec_prot->DebugHits().begin();
        it!=rec_prot->DebugHits().end(); ++it)
  {
    bool prim_res = mad_param_transport_->MADXTeoreticalRPTransversePosition(it->first, 
        pr_vx, pr_vy, 0.0, pr_px, pr_py, pr_pz, prim_pos);
    bool rec_res = mad_param_transport_->MADXTeoreticalRPTransversePosition(it->first, 
        rec_vx, rec_vy, 0.0, rec_p.px(), rec_p.py(), rec_p.pz(), rec_pos);
    
    if(prim_res && rec_res)
    {
      residual_histograms_->GetObj(it->first)->Fill(prim_pos, it->second, rec_pos);
    }
  }
}


void InelasticReconstructionValidation::endJob()
{
  WriteHistograms(hist_file_name_);
}


void InelasticReconstructionValidation::WriteHistograms(const std::string &root_file_name)
{
  TFile *f = TFile::Open(root_file_name.c_str(), "recreate");
  if(!f || !f->IsWritable())
  {
    std::cout<<"Output file not opened correctly!!"<<std::endl;
  }
//  std::cout<<"name:"<<test_profile_.Class()->GetName()<<std::endl;
//  std::cout<<"check sum:"<<test_profile_.Class()->GetCheckSum()<<std::endl;
//  std::cout<<test_profile_.Class()->GetClassInfo()<<std::endl;
//  std::cout<<"version:"<<test_profile_.Class()->GetClassVersion()<<std::endl<<std::endl;
//  
//  std::cout<<"name:"<<TNamed::Class()->GetName()<<std::endl;
//  std::cout<<"check sum:"<<TNamed::Class()->GetCheckSum()<<std::endl;
//  std::cout<<TNamed::Class()->GetClassInfo()<<std::endl;
//  std::cout<<"version:"<<TNamed::Class()->GetClassVersion()<<std::endl<<std::endl;
//  
//  std::cout<<"name:"<<TObject::Class()->GetName()<<std::endl;
//  std::cout<<"check sum:"<<TObject::Class()->GetCheckSum()<<std::endl;
//  std::cout<<TObject::Class()->GetClassInfo()<<std::endl;
//  std::cout<<"version:"<<TObject::Class()->GetClassVersion()<<std::endl<<std::endl;
//  
//  std::cout<<"name:"<<LHCOpticsApproximator::Class()->GetName()<<std::endl;
//  std::cout<<"check sum:"<<LHCOpticsApproximator::Class()->GetCheckSum()<<std::endl;
//  std::cout<<LHCOpticsApproximator::Class()->GetClassInfo()<<std::endl;
//  std::cout<<"version:"<<LHCOpticsApproximator::Class()->GetClassVersion()<<std::endl<<std::endl;
//  
//  std::cout<<"name:"<<ReconstructionProfile::Class()->GetName()<<std::endl;
//  std::cout<<"check sum:"<<ReconstructionProfile::Class()->GetCheckSum()<<std::endl;
//  std::cout<<ReconstructionProfile::Class()->GetClassInfo()<<std::endl;
//  std::cout<<"version:"<<ReconstructionProfile::Class()->GetClassVersion()<<std::endl<<std::endl;
  
  if(verbosity_)
    std::cout<<"Writting histograms..."<<std::endl;
  prot_histograms_->Write(f);
  residual_histograms_->Write(f);
  if(verbosity_)
    std::cout<<"Writting histograms finnished."<<std::endl;
  f->Close();
}

DEFINE_FWK_MODULE(InelasticReconstructionValidation);

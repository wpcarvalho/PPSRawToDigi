#include "TotemRPValidation/ParamMADRefTransport/interface/ParamMADRefTransport.h"
#include "TVector3.h"
#include "TMath.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"


ParamMADRefTransport::ParamMADRefTransport(const edm::ParameterSet& conf, const edm::EventSetup&es)
{
  edm::ESHandle<BeamOpticsParams> BOParH;
  es.get<BeamOpticsParamsRcd>().get(BOParH);
  if(!BOParH.isValid())
    throw cms::Exception("ParamMADRefTransport") << " edm::ESHandle<BeamOpticsParams> is invalid";
  
  beam_energy_ = BOParH->GetBeamEnergy();
  mp0_ = BOParH->GetProtonMass();
  beam_momentum_ = BOParH->GetBeamMomentum();
     
  edm::LogInfo("ParamMADRefTransport") << "[ParamMADRefTransport::ParamMADRefTransport] Constructing object...";
  verbosity_ = conf.getParameter<int>("Verbosity");
  
  param_file_name_220_right_ = conf.getParameter<std::string>("ParameterizationFileName220Right");
  param_file_name_220_left_ = conf.getParameter<std::string>("ParameterizationFileName220Left");
  param_file_name_150_right_ = conf.getParameter<std::string>("ParameterizationFileName150Right");
  param_file_name_150_left_ = conf.getParameter<std::string>("ParameterizationFileName150Left");
  
  param_prefix_220_right_ = conf.getParameter<std::string>("ParameterizationNamePrefix220Right");
  param_prefix_220_left_ = conf.getParameter<std::string>("ParameterizationNamePrefix220Left");
  param_prefix_150_right_ = conf.getParameter<std::string>("ParameterizationNamePrefix150Right");
  param_prefix_150_left_ = conf.getParameter<std::string>("ParameterizationNamePrefix150Left");
  
  right_beam_postfix_ = conf.getParameter<std::string>("RightBeamPostfix");
  left_beam_postfix_ = conf.getParameter<std::string>("LeftBeamPostfix");
  
  InitInverseParametrizationFitter();
}


void ParamMADRefTransport::InitInverseParametrizationFitter()
{
  typedef std::map<RPId, LHCOpticsApproximator> rp_param_map_type;
  rp_param_map_type rp_param_map_right;
  rp_param_map_type rp_param_map_left;
  std::string par_name;
  
  //220 station right
  if(verbosity_)
    std::cout<<"ParamMADRefTransport::InitInverseParametrizationFitter"<<std::cout;


  char *cmsswPath = getenv("CMSSW_BASE");
  std::string fileName_220_r = std::string(cmsswPath) + std::string("/src/") + param_file_name_220_right_;
  TFile *f_220_r = TFile::Open(fileName_220_r.c_str(),"read");
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
//  f_220_r->Close();

  //220 station left
  std::string fileName_220_l = std::string(cmsswPath) + std::string("/src/") + param_file_name_220_left_;
  TFile *f_220_l = TFile::Open(fileName_220_l.c_str(),"read");
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
//  f_220_l->Close();

  //150 station right
  std::string fileName_150_r = std::string(cmsswPath) + std::string("/src/") + param_file_name_150_right_;
  TFile *f_150_r = TFile::Open(fileName_150_r.c_str(),"read");
  par_name = param_prefix_150_right_+"_v_1_"+right_beam_postfix_;
  rp_param_map_right[100] = *((LHCOpticsApproximator*) f_150_r->Get(par_name.c_str()));
  rp_param_map_right[101] = *((LHCOpticsApproximator*) f_150_r->Get(par_name.c_str()));

  par_name = param_prefix_150_right_+"_h_1_"+right_beam_postfix_;
  rp_param_map_right[102] = *((LHCOpticsApproximator*) f_150_r->Get(par_name.c_str()));

  par_name = param_prefix_150_right_+"_h_2_"+right_beam_postfix_;
  rp_param_map_right[103] = *((LHCOpticsApproximator*) f_150_r->Get(par_name.c_str()));

  par_name = param_prefix_150_right_+"_v_2_"+right_beam_postfix_;
  rp_param_map_right[104] = *((LHCOpticsApproximator*) f_150_r->Get(par_name.c_str()));
  rp_param_map_right[105] = *((LHCOpticsApproximator*) f_150_r->Get(par_name.c_str()));
//  f_150_r->Close();

  //150 station left
  std::string fileName_150_l = std::string(cmsswPath) + std::string("/src/") + param_file_name_150_left_;
  TFile *f_150_l = TFile::Open(fileName_150_l.c_str(),"read");
  par_name = param_prefix_150_left_+"_v_1_"+left_beam_postfix_;
  rp_param_map_left[0] = *((LHCOpticsApproximator*) f_150_l->Get(par_name.c_str()));
  rp_param_map_left[1] = *((LHCOpticsApproximator*) f_150_l->Get(par_name.c_str()));

  par_name = param_prefix_150_left_+"_h_1_"+left_beam_postfix_;
  rp_param_map_left[2] = *((LHCOpticsApproximator*) f_150_l->Get(par_name.c_str()));

  par_name = param_prefix_150_left_+"_h_2_"+left_beam_postfix_;
  rp_param_map_left[3] = *((LHCOpticsApproximator*) f_150_l->Get(par_name.c_str()));

  par_name = param_prefix_150_left_+"_v_2_"+left_beam_postfix_;
  rp_param_map_left[4] = *((LHCOpticsApproximator*) f_150_l->Get(par_name.c_str()));
  rp_param_map_left[5] = *((LHCOpticsApproximator*) f_150_l->Get(par_name.c_str()));
//  f_150_l->Close();

  SetParameterizations(rp_param_map_right);
  SetParameterizations(rp_param_map_left);
  
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

void ParamMADRefTransport::SetParameterizations(const transport_to_rp_type& param_map)
{
  transport_to_rp_.insert(param_map.begin(), param_map.end() );
}


//position in [mm], momentum in [GeV]
bool ParamMADRefTransport::MADXTeoreticalRPTransversePosition(RPId id, double vx, double vy, double vz, 
          double px, double py, double pz, 
          TVector2& rp_position, bool include_apertures)
{
  if((pz>0 && id<100) || (pz<0 && id>=100))
    return false;
    
  TVector3 prot_dir_z(px/pz, py/pz, 1.0);
  double tr_x = vx - prot_dir_z.X()*vz;
  double tr_y = vy - prot_dir_z.Y()*vz;
  double theta_x = px/beam_momentum_;
  double theta_y = py/beam_momentum_;
  double ksi = (TMath::Sqrt(px*px+py*py+pz*pz)-beam_momentum_)/beam_momentum_;
  
  double in_par[5];
  in_par[0] = tr_x/1000.;  //convert to [m]
  in_par[1] = theta_x;
  in_par[2] = tr_y/1000.;
  in_par[3] = theta_y;
  in_par[4] = ksi;
  
  double out[2];
  
  transport_to_rp_type::iterator it = transport_to_rp_.find(id);
  if(it==transport_to_rp_.end())
    return false;
  
  bool valid = it->second.Transport2D(in_par, out, include_apertures);
  rp_position.Set(out[0]*1000., out[1]*1000.);
  return valid;
}

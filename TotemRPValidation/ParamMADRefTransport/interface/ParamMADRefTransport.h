#ifndef TotemRPValidation_ParamMADRefTransport_ParamMADRefTransport_h
#define TotemRPValidation_ParamMADRefTransport_ParamMADRefTransport_h
#include "TVector2.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimG4Core/TotemRPProtonTransportParametrization/interface/LHCOpticsApproximator.h"
#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

//Data Formats
#include <iostream>
#include <string>

class ParamMADRefTransport
{
  public:
    typedef std::map<RPId, LHCOpticsApproximator> transport_to_rp_type;
    ParamMADRefTransport(const edm::ParameterSet& conf, const edm::EventSetup&es);
    bool MADXTeoreticalRPTransversePosition(RPId id, double vx, double vy, double vz, 
          double px, double py, double pz, 
          TVector2& rp_position, bool include_apertures=false);
  private:
    void InitInverseParametrizationFitter();
    void SetParameterizations(const transport_to_rp_type& param_map);
    
    int verbosity_;
    std::string param_file_name_220_right_;
    std::string param_file_name_220_left_;
    std::string param_file_name_150_right_;
    std::string param_file_name_150_left_;
    std::string param_prefix_220_right_;
    std::string param_prefix_220_left_;
    std::string param_prefix_150_right_;
    std::string param_prefix_150_left_;
    std::string right_beam_postfix_;
    std::string left_beam_postfix_; 
    
    transport_to_rp_type transport_to_rp_;
    double beam_energy_;  //GeV
    double mp0_;  //GeV
    double beam_momentum_;  //GeV
};


#endif

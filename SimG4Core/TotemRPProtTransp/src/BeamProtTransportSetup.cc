/****************************************************************************
*
* This is a part of the TOTEM offline software.
* Authors:
*    Hubert Niewiadomski
*    Jan Kašpar (jan.kaspar@gmail.com)
*
* $$RCSfile: BeamProtTransportSetup.cc,v $: $
* $Revision: 1.5.2.3 $
* $Date: 2009/12/07 12:50:56 $
*
****************************************************************************/

//#define G4V7
#define DEBUG 0

#include "SimG4Core/TotemRPProtTransp/interface/ProtTranspFastSimModel.h"
#include "SimG4Core/TotemRPProtTransp/interface/BeamProtTransportSetup.h"
#include "TotemCondFormats/DataRecord/interface/ProtonTransportRcd.h"
#include "SimG4CMS/TotemRPProtTranspPar/interface/LHCOpticsApproximator.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include "TotemCondFormats/DataRecord/interface/ProtonTransportRcd.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "Geometry/TotemRecords/interface/RealGeometryRecord.h"

#include "SimG4Core/Geometry/interface/G4LogicalVolumeToDDLogicalPartMap.h" 
#include "DetectorDescription/Core/interface/DDName.h"
#include "DetectorDescription/Core/interface/DDMaterial.h"
#include "DetectorDescription/Core/interface/DDCurrentNamespace.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Electron.hh"
#include "G4FastSimulationManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"

#include <string>

#include "TFile.h"

#include "G4NistManager.hh"
#include "G4RegionStore.hh"
#include <stdlib.h>
#include <thread>

using namespace edm;
using namespace std;


BeamProtTransportSetup* BeamProtTransportSetup::instance = NULL;


BeamProtTransportSetup::BeamProtTransportSetup(const edm::ParameterSet & p) :
  m_pBeamProtTransportSetup(p)
{
  if (instance)
    throw cms::Exception("BeamProtTransportSetup") << "BeamProtTransportSetup is a singleton class and has already been initialized.\n";
  instance = this;

  edm::ParameterSet m_BeamProtTran = p.getParameter<edm::ParameterSet>("BeamProtTransportSetup");
  verbosity_ = m_BeamProtTran.getParameter<bool>("Verbosity");

  std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!! BeamProtTransportSetup::BeamProtTransportSetup !!!!!!!!!!!!!!!"<<std::endl;;

  Beam_IP_150_R_LV = NULL;
  Beam_IP_150_L_LV = NULL;

  model_ip_150_r = NULL;
  model_ip_150_l = NULL;

  Beam_IP_150_R_LV_Name = "Beam_IP_150_R";
  Beam_IP_150_L_LV_Name = "Beam_IP_150_L";

  FindLogicalVolumes();
  BuildTransportModels(p);
}


void BeamProtTransportSetup::FindLogicalVolumes()
{
  //Finding correct G4LogicalVolume for parameterisation
//  ConcreteG4LogicalVolumeToDDLogicalPartMapper::Vector vec =
//  G4LogicalVolumeToDDLogicalPartMapper::instance()->all("volumes");

  edm::LogInfo("TotemRP") << "TotemRP::Proton parameterisation initialization begin !!";

  G4LogicalVolumeStore * theStore = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolumeStore::const_iterator it;
  for (it = theStore->begin(); it != theStore->end(); it++)
  {
    G4LogicalVolume * v = *it;

    if (v->GetName()==Beam_IP_150_R_LV_Name)
    {
      edm::LogInfo("TotemRP") << "TotemRP::Parameterization being initialized for "<<
          v->GetName();
      Beam_IP_150_R_LV = v;
    }
    else if (v->GetName()==Beam_IP_150_L_LV_Name)
    {
      edm::LogInfo("TotemRP") << "TotemRP::Parameterization being initialized for "<<
          v->GetName();
      Beam_IP_150_L_LV = v;
    }
  }
}


void BeamProtTransportSetup::BuildTransportModels(const edm::ParameterSet & p)
{
  edm::ParameterSet m_BeamProtTran = p.getParameter<edm::ParameterSet>("BeamProtTransportSetup");
  std::string param_root_file = m_BeamProtTran.getParameter<std::string>("ModelRootFile");
  std::string model_ip_150_r_name = m_BeamProtTran.getParameter<std::string>("Model_IP_150_R_Name");
  std::string model_ip_150_l_name = m_BeamProtTran.getParameter<std::string>("Model_IP_150_L_Name");

  model_ip_150_r_zmin = m_BeamProtTran.getParameter<double>("Model_IP_150_R_Zmin");
  model_ip_150_r_zmax = m_BeamProtTran.getParameter<double>("Model_IP_150_R_Zmax");
  model_ip_150_l_zmin = m_BeamProtTran.getParameter<double>("Model_IP_150_L_Zmin");
  model_ip_150_l_zmax = m_BeamProtTran.getParameter<double>("Model_IP_150_L_Zmax");

  char *cmsswPath = getenv("CMSSW_BASE");
  std::string fileName = std::string(cmsswPath) + std::string("/src/") + param_root_file;
  
  TFile *f = TFile::Open(fileName.c_str(),"read");
  edm::LogInfo("TotemRP")<<"Root file opened, pointer:"<<f<<std::endl;
  if(!f)
  {
    std::cout<<"BeamProtTransportSetup: File "<<fileName<<" not found. Exiting."<<std::endl;
    exit(0);
  }
  
  LHCOpticsApproximator *aprox_ip_150_r = (LHCOpticsApproximator *) f->Get(model_ip_150_r_name.c_str());
  LHCOpticsApproximator *aprox_ip_150_l = (LHCOpticsApproximator *) f->Get(model_ip_150_l_name.c_str());
  
  if(aprox_ip_150_r == NULL)
  {
    std::cout<<"BeamProtTransportSetup: Parameterisation "<<model_ip_150_r_name<<" missing in file "<<fileName<<std::endl;
    std::cout<<"Job stopped!!"<<std::endl;
    exit(0);
  }
  if(aprox_ip_150_l == NULL)
  {
    std::cout<<"BeamProtTransportSetup: Parameterisation "<<model_ip_150_l_name<<" missing in file "<<fileName<<std::endl;
    std::cout<<"Job stopped!!"<<std::endl;
    exit(0);
  }
    
  edm::LogInfo("TotemRP")<<"Parameterizations read from file, pointers:"<<aprox_ip_150_r<<" "<<aprox_ip_150_l<<" "<<std::endl;

  if(aprox_ip_150_r && aprox_ip_150_l && Beam_IP_150_R_LV && Beam_IP_150_L_LV)
  {
#ifdef G4V7
    model_ip_150_r = new ProtTranspFastSimModel(Beam_IP_150_R_LV_Name,
        Beam_IP_150_R_LV, *aprox_ip_150_r, model_ip_150_r_zmin, model_ip_150_r_zmax, verbosity_);

    model_ip_150_l = new ProtTranspFastSimModel(Beam_IP_150_L_LV_Name,
        Beam_IP_150_L_LV, *aprox_ip_150_l, model_ip_150_l_zmin, model_ip_150_l_zmax, verbosity_);
#else
    G4ProductionCuts *dummyPC = new G4ProductionCuts();
    G4Region *region_ip_150_r = new G4Region(Beam_IP_150_R_LV_Name);
    region_ip_150_r->SetProductionCuts(dummyPC);
    Beam_IP_150_R_LV->SetRegion(region_ip_150_r);
  edm::LogInfo("TotemRP")<<"k7";
//  edm::LogInfo("TotemRP")<<"k8 " << Beam_IP_150_R_LV->GetNoDaughters() << std::endl;
//  edm::LogInfo("TotemRP")<<"k9 " << Beam_IP_150_R_LV->GetMaterial() << std::endl;
//  edm::LogInfo("TotemRP")<<"k10 " << Beam_IP_150_R_LV->GetDaughter(0) << std::endl;
  region_ip_150_r->AddRootLogicalVolume(Beam_IP_150_R_LV);
  edm::LogInfo("TotemRP")<<"k8" <<std::endl;
    region_ip_150_r->SetProductionCuts(dummyPC);
    model_ip_150_r = new ProtTranspFastSimModel(Beam_IP_150_R_LV_Name,
        region_ip_150_r, *aprox_ip_150_r, model_ip_150_r_zmin, model_ip_150_r_zmax, verbosity_);

    G4Region *region_ip_150_l = new G4Region(Beam_IP_150_L_LV_Name);
    region_ip_150_l->SetProductionCuts(dummyPC);
    Beam_IP_150_L_LV->SetRegion(region_ip_150_l);
    region_ip_150_l->AddRootLogicalVolume(Beam_IP_150_L_LV);
    model_ip_150_l = new ProtTranspFastSimModel(Beam_IP_150_L_LV_Name,
        region_ip_150_l, *aprox_ip_150_l, model_ip_150_l_zmin, model_ip_150_l_zmax, verbosity_);
#endif

    edm::LogInfo("TotemRP") << "TotemRP::Fast transport models have been initialized. ";
  }
  else
  {
    edm::LogError("TotemRP") << "TotemRP::Fast transport models failed to be initialized\nLHCApproximators: " <<
		aprox_ip_150_r << ", " << aprox_ip_150_l << ", " << 
		"\nBeam volumes: " << Beam_IP_150_R_LV << ", " << Beam_IP_150_L_LV << ", ";
	// TODO this causes segmentation fault
	throw cms::Exception("TotemRP") << "TotemRP::Fast transport models failed to be initialized.";
//	exit(1);
  }
  f->Close();
}


BeamProtTransportSetup::~BeamProtTransportSetup()
{
#if DEBUG > 0
  printf(">> ~BeamProtTransportSetup\n");
#endif

  if(model_ip_150_r) delete model_ip_150_r;
  if(model_ip_150_l) delete model_ip_150_l;
}


void BeamProtTransportSetup::UpdateSetup(const edm::EventSetup &es)
{

}


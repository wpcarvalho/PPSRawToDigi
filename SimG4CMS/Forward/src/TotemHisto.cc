//#include "Utilities/Configuration/interface/Architecture.h"
#include "SimG4CMS/Forward/interface/TotemHisto.h"
//#include "Utilities/GenUtil/interface/CMSexception.h"

#include <iostream>
#include <cmath>


//#define IT_DEBUG

TotemHisto::TotemHisto(std::string name) {

  std::cout << std::endl << "===>>>  Start booking user histograms with Root" << std::endl;

  nome_file=name;
ntuple = new TNtuple("ntuple", "Ntuple","Event:UnitID:Ptype:TrackID:ParentID:ELoss:PABS:vx:vy:vz:x:y:z:Px:Py:Pz:VPx:VPy:VPz");

  std::cout << std::endl << "===>>> Done booking user histograms and ntuples " << std::endl;

  set_EVT(0);
  set_X(0.);
  set_Y(0.);
  set_Z(0.);
  set_UID( 0);
  set_Ptype( 0) ;
  set_TID( 0) ;
  set_PID(0);
  set_ELoss( 0.) ;
  set_PABS( 0.) ;
  set_VX( 0.);
  set_VY( 0.) ;
  set_VZ( 0.);
  set_PX(0.);
  set_PY(0.);
  set_PZ(0.);
  set_VPX(0.);
  set_VPY(0.);
  set_VPZ(0.);

}
TotemHisto::TotemHisto() {

  std::cout << std::endl << "===>>>  Start booking user histograms with Root" << std::endl;

  nome_file="TotemHits.root";
ntuple = new TNtuple("ntuple", "Ntuple","Event:UnitID:Ptype:TrackID:ParentID:ELoss:PABS:vx:vy:vz:x:y:z:Px:Py:Pz:VPx:VPy:VPz");

  std::cout << std::endl << "===>>> Done booking user histograms and ntuples " << std::endl;

  set_EVT(0);
  set_X(0.);
  set_Y(0.);
  set_Z(0.);
  set_UID( 0);
  set_Ptype( 0) ;
  set_TID( 0) ;
  set_PID(0);
  set_ELoss( 0.) ;
  set_PABS( 0.) ;
  set_VX( 0.);
  set_VY( 0.) ;
  set_VZ( 0.);
  set_PX(0.);
  set_PY(0.);
  set_PZ(0.);
  set_VPX(0.);
  set_VPY(0.);
  set_VPZ(0.);

}


TotemHisto::~TotemHisto() {

  std::cout << "========================================================" << std::endl;  
  std::cout << "=== TotemHisto: Start writing user histograms ===" << std::endl;
  const char* c_nome_file;
  c_nome_file=nome_file.c_str();
 
  TFile rt_hf(c_nome_file,"RECREATE");
  rt_hf.SetCompressionLevel(2);

  std::cout << "I have already created root file" << std::endl;

 

  TotemHisto::ntuple->Write();

  rt_hf.Close();  

  std::cout << std::endl << "TotemHisto: End writing user histograms " << std::endl;

}


void TotemHisto::fillNtuple(){


#ifdef IT_DEBUG 
 std::cout << "RIEMPO L'NTUPLA DI ROOT " <<std::endl;
#endif
  rootvec[0]=evt;
  rootvec[1]=UID;
  rootvec[2]=Ptype;
  rootvec[3]=TID;
  rootvec[4]=PID;
  rootvec[5]=ELoss;
  rootvec[6]=PABS;
  rootvec[7]=vx;
  rootvec[8]=vy;
  rootvec[9]=vz;
  rootvec[10]=x;
  rootvec[11]=y;
  rootvec[12]=z;
  rootvec[13]=Px;
  rootvec[14]=Py;
  rootvec[15]=Pz;
  rootvec[16]=VPx;
  rootvec[17]=VPy;
  rootvec[18]=VPz;

  if(ntuple){
    ntuple->Fill(&(rootvec[0]));
  }
}

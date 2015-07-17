#include <memory>
/*
  Created by Fabrizio Ferro - INFN Genova for TOTEM
*/

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "DataFormats/T1Road/interface/T1Road.h"
#include "DataFormats/T1Road/interface/T1RecHitGlobal.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/Provenance/interface/Provenance.h"
#include "DataFormats/Provenance/interface/BranchDescription.h"
#include "DataFormats/T1Road/interface/T1Road.h"
#include "DataFormats/T1T2Track/interface/T1T2TrackCollection.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"
#include "TotemAlignment/T1Alignment/interface/func3.h"

//
// class decleration
//


#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TMatrixD.h"


class T1Sex2 : public edm::EDAnalyzer 
{
 

 public:
  explicit T1Sex2(const edm::ParameterSet&);
  ~T1Sex2();
  float Eta(float ,float ,float );
   float Phi(float,float);
   int sextant(float ,float ,float );
   int plane(float);
 private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  vector<TVector3> _plus;
  vector<TVector2> _plus_prev;
  vector<TVector2> _plus_error;

  vector<TVector3> _minus;
  vector<TVector2> _minus_prev;
  vector<TVector2> _minus_error;




  TFile* theFile; 
  int _Verbosity;

  TH1D * hDDDx_0_0_01 ; 
  TH1D * hDDDy_0_0_01 ; 
  
  TH1D * hDDDx_0_0_12 ; 
  TH1D * hDDDy_0_0_12 ; 
  
  TH1D * hDDDx_0_0_23 ; 
  TH1D * hDDDy_0_0_23 ; 
  
  TH1D * hDDDx_0_0_34 ; 
  TH1D * hDDDy_0_0_34 ; 
  
  TH1D * hDDDx_0_0_45 ; 
  TH1D * hDDDy_0_0_45 ; 
  
  TH1D * hDDDx_0_0_50 ; 
  TH1D * hDDDy_0_0_50 ; 
  
    TH1D * hDDDx_0_1_01 ; 
  TH1D * hDDDy_0_1_01 ; 
  
  TH1D * hDDDx_0_1_12 ; 
  TH1D * hDDDy_0_1_12 ; 
  
  TH1D * hDDDx_0_1_23 ; 
  TH1D * hDDDy_0_1_23 ; 
  
  TH1D * hDDDx_0_1_34 ; 
  TH1D * hDDDy_0_1_34 ; 
  
  TH1D * hDDDx_0_1_45 ; 
  TH1D * hDDDy_0_1_45 ; 
  
  TH1D * hDDDx_0_1_50 ; 
  TH1D * hDDDy_0_1_50 ; 
  
  TH1D * hDDDx_0_2_01 ; 
  TH1D * hDDDy_0_2_01 ; 
  
  TH1D * hDDDx_0_2_12 ; 
  TH1D * hDDDy_0_2_12 ; 
  
  TH1D * hDDDx_0_2_23 ; 
  TH1D * hDDDy_0_2_23 ; 
  
  TH1D * hDDDx_0_2_34 ; 
  TH1D * hDDDy_0_2_34 ; 
  
  TH1D * hDDDx_0_2_45 ; 
  TH1D * hDDDy_0_2_45 ; 
  
  TH1D * hDDDx_0_2_50 ; 
  TH1D * hDDDy_0_2_50 ; 
  

  TH1D * hDDDx_0_3_01 ; 
  TH1D * hDDDy_0_3_01 ; 
  
  TH1D * hDDDx_0_3_12 ; 
  TH1D * hDDDy_0_3_12 ; 
  
  TH1D * hDDDx_0_3_23 ; 
  TH1D * hDDDy_0_3_23 ; 
  
  TH1D * hDDDx_0_3_34 ; 
  TH1D * hDDDy_0_3_34 ; 
  
  TH1D * hDDDx_0_3_45 ; 
  TH1D * hDDDy_0_3_45 ; 
  
  TH1D * hDDDx_0_3_50 ; 
  TH1D * hDDDy_0_3_50 ; 
  
  TH1D * hDDDx_0_4_01 ; 
  TH1D * hDDDy_0_4_01 ; 
  
  TH1D * hDDDx_0_4_12 ; 
  TH1D * hDDDy_0_4_12 ; 
  
  TH1D * hDDDx_0_4_23 ; 
  TH1D * hDDDy_0_4_23 ; 
  
  TH1D * hDDDx_0_4_34 ; 
  TH1D * hDDDy_0_4_34 ; 
  
  TH1D * hDDDx_0_4_45 ; 
  TH1D * hDDDy_0_4_45 ; 
  
  TH1D * hDDDx_0_4_50 ; 
  TH1D * hDDDy_0_4_50 ; 
  





  TH1D * hDDDx_1_0_01 ; 
  TH1D * hDDDy_1_0_01 ; 
  
  TH1D * hDDDx_1_0_12 ; 
  TH1D * hDDDy_1_0_12 ; 
  
  TH1D * hDDDx_1_0_23 ; 
  TH1D * hDDDy_1_0_23 ; 
  
  TH1D * hDDDx_1_0_34 ; 
  TH1D * hDDDy_1_0_34 ; 
  
  TH1D * hDDDx_1_0_45 ; 
  TH1D * hDDDy_1_0_45 ; 
  
  TH1D * hDDDx_1_0_50 ; 
  TH1D * hDDDy_1_0_50 ; 
  
    TH1D * hDDDx_1_1_01 ; 
  TH1D * hDDDy_1_1_01 ; 
  
  TH1D * hDDDx_1_1_12 ; 
  TH1D * hDDDy_1_1_12 ; 
  
  TH1D * hDDDx_1_1_23 ; 
  TH1D * hDDDy_1_1_23 ; 
  
  TH1D * hDDDx_1_1_34 ; 
  TH1D * hDDDy_1_1_34 ; 
  
  TH1D * hDDDx_1_1_45 ; 
  TH1D * hDDDy_1_1_45 ; 
  
  TH1D * hDDDx_1_1_50 ; 
  TH1D * hDDDy_1_1_50 ; 
  
  TH1D * hDDDx_1_2_01 ; 
  TH1D * hDDDy_1_2_01 ; 
  
  TH1D * hDDDx_1_2_12 ; 
  TH1D * hDDDy_1_2_12 ; 
  
  TH1D * hDDDx_1_2_23 ; 
  TH1D * hDDDy_1_2_23 ; 
  
  TH1D * hDDDx_1_2_34 ; 
  TH1D * hDDDy_1_2_34 ; 
  
  TH1D * hDDDx_1_2_45 ; 
  TH1D * hDDDy_1_2_45 ; 
  
  TH1D * hDDDx_1_2_50 ; 
  TH1D * hDDDy_1_2_50 ; 
  

  TH1D * hDDDx_1_3_01 ; 
  TH1D * hDDDy_1_3_01 ; 
  
  TH1D * hDDDx_1_3_12 ; 
  TH1D * hDDDy_1_3_12 ; 
  
  TH1D * hDDDx_1_3_23 ; 
  TH1D * hDDDy_1_3_23 ; 
  
  TH1D * hDDDx_1_3_34 ; 
  TH1D * hDDDy_1_3_34 ; 
  
  TH1D * hDDDx_1_3_45 ; 
  TH1D * hDDDy_1_3_45 ; 
  
  TH1D * hDDDx_1_3_50 ; 
  TH1D * hDDDy_1_3_50 ; 
  
  TH1D * hDDDx_1_4_01 ; 
  TH1D * hDDDy_1_4_01 ; 
  
  TH1D * hDDDx_1_4_12 ; 
  TH1D * hDDDy_1_4_12 ; 
  
  TH1D * hDDDx_1_4_23 ; 
  TH1D * hDDDy_1_4_23 ; 
  
  TH1D * hDDDx_1_4_34 ; 
  TH1D * hDDDy_1_4_34 ; 
  
  TH1D * hDDDx_1_4_45 ; 
  TH1D * hDDDy_1_4_45 ; 
  
  TH1D * hDDDx_1_4_50 ; 
  TH1D * hDDDy_1_4_50 ; 



};

// -*- C++ -*-
//
// Package:    T1Analysis2
// Class:      T1Analysis2
// 
/**\class T1Analysis2 T1Analysis2.cc TotemAnalysis/T1Analysis2/src/T1Analysis2.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  fabrizio ferro
//         Created:  Mon Nov 28 09:33:34 CET 2011
// $Id$
//
//


// system include files
#include <memory>
#include <map>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "TotemRawDataLibrary/DataFormats/interface/RawEvent.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/T1Road/interface/T1RecHitGlobal.h"
#include "DataFormats/T1T2Track/interface/T1T2TrackCollection.h"

#include <boost/shared_ptr.hpp>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

#define PI 3.141592654

//
// class declaration
//

class T1Analysis3 : public edm::EDAnalyzer {
public:
  explicit T1Analysis3(const edm::ParameterSet&);
  ~T1Analysis3();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

// ----------member data ---------------------------

  float EFF[2][120][4];
  float EFFplus[120][4];
  float EFFminus[120][4];

  std::auto_ptr<TH1D> hDeta;
  std::auto_ptr<TH1D> hDphi;

  std::auto_ptr<TH1D> hPhiPlus1;
  std::auto_ptr<TH1D> hPhiPlus2;
  std::auto_ptr<TH1D> hPhiPlus3;
  std::auto_ptr<TH1D> hPhiPlus4;

  std::auto_ptr<TH1D> hPhiPlus1near;
  std::auto_ptr<TH1D> hPhiPlus2near;
  std::auto_ptr<TH1D> hPhiPlus3near;
  std::auto_ptr<TH1D> hPhiPlus4near;

  std::auto_ptr<TH1D> hPhiMinus1;
  std::auto_ptr<TH1D> hPhiMinus2;
  std::auto_ptr<TH1D> hPhiMinus3;
  std::auto_ptr<TH1D> hPhiMinus4;

  std::auto_ptr<TH1D> hPhiMinus1near;
  std::auto_ptr<TH1D> hPhiMinus2near;
  std::auto_ptr<TH1D> hPhiMinus3near;
  std::auto_ptr<TH1D> hPhiMinus4near;

  std::auto_ptr<TH1D> hPhiDivPlus1;
  std::auto_ptr<TH1D> hPhiDivPlus2;
  std::auto_ptr<TH1D> hPhiDivPlus3;
  std::auto_ptr<TH1D> hPhiDivPlus4;

  std::auto_ptr<TH1D> hPhiDivPlus1near;
  std::auto_ptr<TH1D> hPhiDivPlus2near;
  std::auto_ptr<TH1D> hPhiDivPlus3near;
  std::auto_ptr<TH1D> hPhiDivPlus4near;

  std::auto_ptr<TH1D> hPhiDivMinus1;
  std::auto_ptr<TH1D> hPhiDivMinus2;
  std::auto_ptr<TH1D> hPhiDivMinus3;
  std::auto_ptr<TH1D> hPhiDivMinus4;

  std::auto_ptr<TH1D> hPhiDivMinus1near;
  std::auto_ptr<TH1D> hPhiDivMinus2near;
  std::auto_ptr<TH1D> hPhiDivMinus3near;
  std::auto_ptr<TH1D> hPhiDivMinus4near;



  std::auto_ptr<TH2D> 	hDetaRmin;
  std::auto_ptr<TH2D>	hDetaChi;
  std::auto_ptr<TH2D>	hDetaZ;

  std::auto_ptr<TH2D>	hDphiRmin;
  std::auto_ptr<TH2D>	hDphiChi;
  std::auto_ptr<TH2D>	hDphiZ;
  std::auto_ptr<TH2D>	hDphiPhi;

  TH1D * effplus1 ; 
  TH1D * effminus1 ;
  TH1D * effplus2 ; 
  TH1D * effminus2 ;
  TH1D * effplus3 ; 
  TH1D * effminus3 ;
  TH1D * effplus4 ; 
  TH1D * effminus4 ;

  TFile* theFile;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
T1Analysis3::T1Analysis3(const edm::ParameterSet& iConfig)

{
//now do what ever initialization is needed

}


T1Analysis3::~T1Analysis3()
{
 
// do anything here that needs to be done at desctruction time
// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
T1Analysis3::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;



//

// T1 data

  edm::Handle<T1T2TrackCollection> T1trackCollection;
  iEvent.getByLabel("t1tracks2","T1TrackColl",T1trackCollection);

  T1T2TrackCollection::const_iterator T1TC_it;

// tagli su numero di tracce e chiquadro? (50)

//  if(T1trackCollection->size() < 40)
  for(T1TC_it=T1trackCollection->begin(); T1TC_it!=T1trackCollection->end(); T1TC_it++){
    float Deta_max=0;
    float Dphi_max=0;
//    if( (*T1TC_it).ChiSquared() < 50)
      for(unsigned int i=0; i<(*T1TC_it).GetHitEntries(); i++){

	float Deta = (*T1TC_it).Eta() - (*T1TC_it).GetHitT1(i).eta();
	float Dphi = (*T1TC_it).Phi() - (*T1TC_it).GetHitT1(i).phi();
	if(fabs(Dphi)>PI && Dphi>0)Dphi=-2*PI + Dphi;
	if(fabs(Dphi)>PI && Dphi<0)Dphi=2*PI + Dphi;
	if(fabs(Deta)>Deta_max)Deta_max=fabs(Deta);
	if(fabs(Dphi)>Dphi_max)Dphi_max=fabs(Dphi);

	hDeta->Fill(Deta);
	hDphi->Fill(Dphi);

	hDetaRmin->Fill(Deta,(*T1TC_it).Rmin());
	hDetaChi->Fill(Deta,(*T1TC_it).ChiSquaredOverN());
	hDetaZ->Fill(Deta,(*T1TC_it).Z_at_Rmin());

	hDphiRmin->Fill(Dphi,(*T1TC_it).Rmin());
	hDphiChi->Fill(Dphi,(*T1TC_it).ChiSquaredOverN());
	hDphiZ->Fill(Dphi,(*T1TC_it).Z_at_Rmin());
	hDphiPhi->Fill(Dphi,(*T1TC_it).Phi());
      }
// multicounting delle tracce?
	if((*T1TC_it).GetHitT1(0).eta() > 3. && (*T1TC_it).GetHitT1(0).eta()<3.4){
	  hPhiPlus1->Fill( (*T1TC_it).GetHitT1(0).phi() );
	  hPhiDivPlus1->Fill( (*T1TC_it).GetHitT1(0).phi() );
	  if(Deta_max<0.4 && Dphi_max<0.05){hPhiPlus1near->Fill( (*T1TC_it).GetHitT1(0).phi() );
hPhiDivPlus1near->Fill( (*T1TC_it).GetHitT1(0).phi() );
	}
	}

	if((*T1TC_it).GetHitT1(0).eta() > 3.4 && (*T1TC_it).GetHitT1(0).eta()<3.8){
	  hPhiPlus2->Fill( (*T1TC_it).GetHitT1(0).phi() );
	  hPhiDivPlus2->Fill( (*T1TC_it).GetHitT1(0).phi() );
	  if(Deta_max<0.4 && Dphi_max<0.05){hPhiPlus2near->Fill( (*T1TC_it).GetHitT1(0).phi() );
hPhiDivPlus2near->Fill( (*T1TC_it).GetHitT1(0).phi() );
	}
	}
	if((*T1TC_it).GetHitT1(0).eta() > 3.8 && (*T1TC_it).GetHitT1(0).eta()<4.2){
	  hPhiPlus3->Fill( (*T1TC_it).GetHitT1(0).phi() );
	  hPhiDivPlus3->Fill( (*T1TC_it).GetHitT1(0).phi() );
	  if(Deta_max<0.4 && Dphi_max<0.05){hPhiPlus3near->Fill( (*T1TC_it).GetHitT1(0).phi() );
hPhiDivPlus3near->Fill( (*T1TC_it).GetHitT1(0).phi() );
	}
	}
	if((*T1TC_it).GetHitT1(0).eta() > 4.2 && (*T1TC_it).GetHitT1(0).eta()<4.6){
	  hPhiPlus4->Fill( (*T1TC_it).GetHitT1(0).phi() );
	  hPhiDivPlus4->Fill( (*T1TC_it).GetHitT1(0).phi() );
	  if(Deta_max<0.4 && Dphi_max<0.05){hPhiPlus4near->Fill( (*T1TC_it).GetHitT1(0).phi() );
hPhiDivPlus4near->Fill( (*T1TC_it).GetHitT1(0).phi() );
	}
	}

	if((*T1TC_it).GetHitT1(0).eta() < -3. && (*T1TC_it).GetHitT1(0).eta()>-3.4){
	  hPhiMinus1->Fill( (*T1TC_it).GetHitT1(0).phi() );
	  hPhiDivMinus1->Fill( (*T1TC_it).GetHitT1(0).phi() );
	  if(Deta_max<0.4 && Dphi_max<0.05){hPhiMinus1near->Fill( (*T1TC_it).GetHitT1(0).phi() );
hPhiDivMinus1near->Fill( (*T1TC_it).GetHitT1(0).phi() );
	}
	}

	if((*T1TC_it).GetHitT1(0).eta() < -3.4 && (*T1TC_it).GetHitT1(0).eta()>-3.8){
	  hPhiMinus2->Fill( (*T1TC_it).GetHitT1(0).phi() );
	  hPhiDivMinus2->Fill( (*T1TC_it).GetHitT1(0).phi() );
	  if(Deta_max<0.4 && Dphi_max<0.05){hPhiMinus2near->Fill( (*T1TC_it).GetHitT1(0).phi() );
hPhiDivMinus2near->Fill( (*T1TC_it).GetHitT1(0).phi() );
	}
	}
	if((*T1TC_it).GetHitT1(0).eta() < -3.8 && (*T1TC_it).GetHitT1(0).eta()>-4.2){
	  hPhiMinus3->Fill( (*T1TC_it).GetHitT1(0).phi() );
	  hPhiDivMinus3->Fill( (*T1TC_it).GetHitT1(0).phi() );
	  if(Deta_max<0.4 && Dphi_max<0.05){hPhiMinus3near->Fill( (*T1TC_it).GetHitT1(0).phi() );
hPhiDivMinus3near->Fill( (*T1TC_it).GetHitT1(0).phi() );
	}
	}
	if((*T1TC_it).GetHitT1(0).eta() < -4.2 && (*T1TC_it).GetHitT1(0).eta()>-4.6){
	  hPhiMinus4->Fill( (*T1TC_it).GetHitT1(0).phi() );
	  hPhiDivMinus4->Fill( (*T1TC_it).GetHitT1(0).phi() );
	  if(Deta_max<0.4 && Dphi_max<0.05){hPhiMinus4near->Fill( (*T1TC_it).GetHitT1(0).phi() );
hPhiDivMinus4near->Fill( (*T1TC_it).GetHitT1(0).phi() );
	}
	}


      
      


  }










 

}


// ------------ method called once each job just before starting event loop  ------------
void 
T1Analysis3::beginJob()
{


 
 

  hDeta = std::auto_ptr<TH1D>(new TH1D("hDeta","hDeta",100,-2*PI,2*PI));
  hDeta->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hDeta","hDeta",100,-0.5,99.5)

  hDphi = std::auto_ptr<TH1D>(new TH1D("hDphi","hDphi",100,-2*PI,2*PI));
  hDphi->SetDirectory(0); //auto_ptr<TH1D>(new TH1D("hDphi","hDphi",100,-0.5,99.5)

  hDetaRmin = std::auto_ptr<TH2D>(new TH2D("DEtaRmin","DEtaRmin",100,-2*PI,2*PI,100,0,1000));
  hDetaChi = std::auto_ptr<TH2D>(new TH2D("DEtaChi","DEtaChi",100,-2*PI,2*PI,100,0,30));
  hDetaZ = std::auto_ptr<TH2D>(new TH2D("DEtaZ","DEtaZ",100,-2*PI,2*PI,100,-10000,10000));

  hDphiRmin = std::auto_ptr<TH2D>(new TH2D("DPhiRmin","DPhiRmin",100,-2*PI,2*PI,100,0,1000));
  hDphiChi = std::auto_ptr<TH2D>(new TH2D("DPhiChi","DPhiChi",100,-2*PI,2*PI,100,0,30));
  hDphiZ = std::auto_ptr<TH2D>(new TH2D("DPhiZ","DPhiZ",100,-2*PI,2*PI,100,-10000,10000));
  hDphiPhi = std::auto_ptr<TH2D>(new TH2D("DPhiPhi","DPhiPhi",100,-2*PI,2*PI,100,0,2*PI));

  hPhiPlus1 =  std::auto_ptr<TH1D>(new TH1D("hPhiPlus1","hPhiPlus1",120,0.,2*PI));
  hPhiPlus2 =  std::auto_ptr<TH1D>(new TH1D("hPhiPlus2","hPhiPlus2",120,0.,2*PI));
  hPhiPlus3 =  std::auto_ptr<TH1D>(new TH1D("hPhiPlus3","hPhiPlus3",120,0.,2*PI));
  hPhiPlus4 =  std::auto_ptr<TH1D>(new TH1D("hPhiPlus4","hPhiPlus4",120,0.,2*PI));

  hPhiPlus1near =  std::auto_ptr<TH1D>(new TH1D("hPhiPlus1near","hPhiPlus1near",120,0.,2*PI));
  hPhiPlus2near =  std::auto_ptr<TH1D>(new TH1D("hPhiPlus2near","hPhiPlus2near",120,0.,2*PI));
  hPhiPlus3near =  std::auto_ptr<TH1D>(new TH1D("hPhiPlus3near","hPhiPlus3near",120,0.,2*PI));
  hPhiPlus4near =  std::auto_ptr<TH1D>(new TH1D("hPhiPlus4near","hPhiPlus4near",120,0.,2*PI));

  hPhiMinus1 =  std::auto_ptr<TH1D>(new TH1D("hPhiMinus1","hPhiMinus1",120,0.,2*PI));
  hPhiMinus2 =  std::auto_ptr<TH1D>(new TH1D("hPhiMinus2","hPhiMinus2",120,0.,2*PI));
  hPhiMinus3 =  std::auto_ptr<TH1D>(new TH1D("hPhiMinus3","hPhiMinus3",120,0.,2*PI));
  hPhiMinus4 =  std::auto_ptr<TH1D>(new TH1D("hPhiMinus4","hPhiMinus4",120,0.,2*PI));

  hPhiMinus1near =  std::auto_ptr<TH1D>(new TH1D("hPhiMinus1near","hPhiMinus1near",120,0.,2*PI));
  hPhiMinus2near =  std::auto_ptr<TH1D>(new TH1D("hPhiMinus2near","hPhiMinus2near",120,0.,2*PI));
  hPhiMinus3near =  std::auto_ptr<TH1D>(new TH1D("hPhiMinus3near","hPhiMinus3near",120,0.,2*PI));
  hPhiMinus4near =  std::auto_ptr<TH1D>(new TH1D("hPhiMinus4near","hPhiMinus4near",120,0.,2*PI));

  hPhiDivPlus1 =  std::auto_ptr<TH1D>(new TH1D("hPhiDivPlus1","hPhiDivPlus1",120,0.,2*PI));
  hPhiDivPlus2 =  std::auto_ptr<TH1D>(new TH1D("hPhiDivPlus2","hPhiDivPlus2",120,0.,2*PI));
  hPhiDivPlus3 =  std::auto_ptr<TH1D>(new TH1D("hPhiDivPlus3","hPhiDivPlus3",120,0.,2*PI));
  hPhiDivPlus4 =  std::auto_ptr<TH1D>(new TH1D("hPhiDivPlus4","hPhiDivPlus4",120,0.,2*PI));

  hPhiDivPlus1near =  std::auto_ptr<TH1D>(new TH1D("hPhiDivPlus1near","hPhiDivPlus1near",120,0.,2*PI));
  hPhiDivPlus2near =  std::auto_ptr<TH1D>(new TH1D("hPhiDivPlus2near","hPhiDivPlus2near",120,0.,2*PI));
  hPhiDivPlus3near =  std::auto_ptr<TH1D>(new TH1D("hPhiDivPlus3near","hPhiDivPlus3near",120,0.,2*PI));
  hPhiDivPlus4near =  std::auto_ptr<TH1D>(new TH1D("hPhiDivPlus4near","hPhiDivPlus4near",120,0.,2*PI));

  hPhiDivMinus1 =  std::auto_ptr<TH1D>(new TH1D("hPhiDivMinus1","hPhiDivMinus1",120,0.,2*PI));
  hPhiDivMinus2 =  std::auto_ptr<TH1D>(new TH1D("hPhiDivMinus2","hPhiDivMinus2",120,0.,2*PI));
  hPhiDivMinus3 =  std::auto_ptr<TH1D>(new TH1D("hPhiDivMinus3","hPhiDivMinus3",120,0.,2*PI));
  hPhiDivMinus4 =  std::auto_ptr<TH1D>(new TH1D("hPhiDivMinus4","hPhiDivMinus4",120,0.,2*PI));

  hPhiDivMinus1near =  std::auto_ptr<TH1D>(new TH1D("hPhiDivMinus1near","hPhiDivMinus1near",120,0.,2*PI));
  hPhiDivMinus2near =  std::auto_ptr<TH1D>(new TH1D("hPhiDivMinus2near","hPhiDivMinus2near",120,0.,2*PI));
  hPhiDivMinus3near =  std::auto_ptr<TH1D>(new TH1D("hPhiDivMinus3near","hPhiDivMinus3near",120,0.,2*PI));
  hPhiDivMinus4near =  std::auto_ptr<TH1D>(new TH1D("hPhiDivMinus4near","hPhiDivMinus4near",120,0.,2*PI));

  effplus1 = (new TH1D("effplus1","effplus1",120,0,2*PI));
  effminus1 = (new TH1D("effminus1","effminus1",120,0,2*PI));
  effplus2 = (new TH1D("effplus2","effplus2",120,0,2*PI));
  effminus2 = (new TH1D("effminus2","effminus2",120,0,2*PI));
  effplus3 = (new TH1D("effplus3","effplus3",120,0,2*PI));
  effminus3 = (new TH1D("effminus3","effminus3",120,0,2*PI));
  effplus4 = (new TH1D("effplus4","effplus4",120,0,2*PI));
  effminus4 = (new TH1D("effminus4","effminus4",120,0,2*PI));


  hDetaRmin->SetDirectory(0);
  hDetaChi->SetDirectory(0);
  hDetaZ->SetDirectory(0);

  hDphiRmin->SetDirectory(0);
  hDphiChi->SetDirectory(0);
  hDphiZ->SetDirectory(0);
  hDphiPhi->SetDirectory(0);

  hPhiPlus1 ->SetDirectory(0); 
  hPhiPlus2 ->SetDirectory(0); 
  hPhiPlus3 ->SetDirectory(0); 
  hPhiPlus4 ->SetDirectory(0); 

  hPhiPlus1near ->SetDirectory(0); 
  hPhiPlus2near ->SetDirectory(0); 
  hPhiPlus3near ->SetDirectory(0); 
  hPhiPlus4near ->SetDirectory(0); 

  hPhiMinus1 ->SetDirectory(0); 
  hPhiMinus2 ->SetDirectory(0); 
  hPhiMinus3 ->SetDirectory(0); 
  hPhiMinus4 ->SetDirectory(0); 

  hPhiMinus1near ->SetDirectory(0); 
  hPhiMinus2near ->SetDirectory(0); 
  hPhiMinus3near ->SetDirectory(0); 
  hPhiMinus4near ->SetDirectory(0); 

  hPhiDivPlus1 ->SetDirectory(0); 
  hPhiDivPlus2 ->SetDirectory(0); 
  hPhiDivPlus3 ->SetDirectory(0); 
  hPhiDivPlus4 ->SetDirectory(0); 

  hPhiDivPlus1near ->SetDirectory(0); 
  hPhiDivPlus2near ->SetDirectory(0); 
  hPhiDivPlus3near ->SetDirectory(0); 
  hPhiDivPlus4near ->SetDirectory(0); 

  hPhiDivMinus1 ->SetDirectory(0); 
  hPhiDivMinus2 ->SetDirectory(0); 
  hPhiDivMinus3 ->SetDirectory(0); 
  hPhiDivMinus4 ->SetDirectory(0); 

  hPhiDivMinus1near ->SetDirectory(0); 
  hPhiDivMinus2near ->SetDirectory(0); 
  hPhiDivMinus3near ->SetDirectory(0); 
  hPhiDivMinus4near ->SetDirectory(0); 

  effplus1 ->SetDirectory(0); 
  effminus1 ->SetDirectory(0);
  effplus2 ->SetDirectory(0); 
  effminus2 ->SetDirectory(0);
  effplus3 ->SetDirectory(0); 
  effminus3 ->SetDirectory(0);
  effplus4 ->SetDirectory(0); 
  effminus4 ->SetDirectory(0);





  ifstream file;
  file.open("EFF_bin.dat");
  for(unsigned int j=0;j<2;j++){
    
    for(unsigned int i = 0; i < 120; i++){
      for(unsigned int l=0; l<4;l++){
	
	file >> EFF[j][i][l];
	if(!file.good()) break;
      }
    }
  }

  file.close();

  for(unsigned int i=0;i<120;i++)
    for(unsigned int j=0;j<4;j++){
      
      EFFplus[i][j] = EFF[0][i][j];
      EFFminus[i][j] = EFF[1][i][j];

    }
    

  for(unsigned int i=0;i<120;i++){
    float fi = (float)i*0.0523333 + 0.052333/2.;
    
    effplus1->Fill(fi,EFFplus[i][0]);
    effminus1->Fill(fi,EFFminus[i][0]);
    effplus2->Fill(fi,EFFplus[i][1]);
    effminus2->Fill(fi,EFFminus[i][1]);
    effplus3->Fill(fi,EFFplus[i][2]);
    effminus3->Fill(fi,EFFminus[i][2]);
    effplus4->Fill(fi,EFFplus[i][3]);
    effminus4->Fill(fi,EFFminus[i][3]);


  }


}
// ------------ method called once each job just after ending the event loop  ------------
void 
T1Analysis3::endJob() 
{

  theFile = TFile::Open("T1Analysis3.root", "recreate");
  if(!theFile || !theFile->IsWritable())
    {
      std::cout<<"Output file not opened correctly!!"<<std::endl;
    }

  hDeta->Write(); 
  hDphi->Write(); 


  hDetaRmin->Write();
  hDetaChi->Write();
  hDetaZ->Write();

  hDphiRmin->Write();
  hDphiChi->Write();
  hDphiZ->Write();
  hDphiPhi->Write();

  hPhiPlus1 ->Write(); 
  hPhiPlus2 ->Write(); 
  hPhiPlus3 ->Write(); 
  hPhiPlus4 ->Write(); 

  hPhiPlus1near ->Write(); 
  hPhiPlus2near ->Write(); 
  hPhiPlus3near ->Write(); 
  hPhiPlus4near ->Write(); 

  hPhiMinus1 ->Write(); 
  hPhiMinus2 ->Write(); 
  hPhiMinus3 ->Write(); 
  hPhiMinus4 ->Write(); 

  hPhiMinus1near ->Write(); 
  hPhiMinus2near ->Write(); 
  hPhiMinus3near ->Write(); 
  hPhiMinus4near ->Write(); 

 

  hPhiDivPlus1 ->Divide(effplus1); 
  hPhiDivPlus2 ->Divide(effplus2); 
  hPhiDivPlus3 ->Divide(effplus3); 
  hPhiDivPlus4 ->Divide(effplus4); 

  hPhiDivPlus1near ->Divide(effplus1); 
  hPhiDivPlus2near ->Divide(effplus2); 
  hPhiDivPlus3near ->Divide(effplus3); 
  hPhiDivPlus4near ->Divide(effplus4); 

  hPhiDivMinus1 ->Divide(effminus1); 
  hPhiDivMinus2 ->Divide(effminus2); 
  hPhiDivMinus3 ->Divide(effminus3); 
  hPhiDivMinus4 ->Divide(effminus4); 

  hPhiDivMinus1near ->Divide(effminus1); 
  hPhiDivMinus2near ->Divide(effminus2); 
  hPhiDivMinus3near ->Divide(effminus3); 
  hPhiDivMinus4near ->Divide(effminus4); 



  hPhiDivPlus1 ->Write(); 
  hPhiDivPlus2 ->Write(); 
  hPhiDivPlus3 ->Write(); 
  hPhiDivPlus4 ->Write(); 

  hPhiDivPlus1near ->Write(); 
  hPhiDivPlus2near ->Write(); 
  hPhiDivPlus3near ->Write(); 
  hPhiDivPlus4near ->Write(); 

  hPhiDivMinus1 ->Write(); 
  hPhiDivMinus2 ->Write(); 
  hPhiDivMinus3 ->Write(); 
  hPhiDivMinus4 ->Write(); 

  hPhiDivMinus1near ->Write(); 
  hPhiDivMinus2near ->Write(); 
  hPhiDivMinus3near ->Write(); 
  hPhiDivMinus4near ->Write(); 

  effplus1 ->Write(); 
  effminus1 ->Write();
  effplus2 ->Write(); 
  effminus2 ->Write();
  effplus3 ->Write(); 
  effminus3 ->Write();
  effplus4 ->Write(); 
  effminus4 ->Write();






  theFile->Close();





}

// ------------ method called when starting to processes a run  ------------
void 
T1Analysis3::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
T1Analysis3::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
T1Analysis3::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
T1Analysis3::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
T1Analysis3::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
//The following says we do not know what parameters are allowed so do no validation
// Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(T1Analysis3);

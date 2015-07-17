/*
  Created by Fabrizio Ferro - INFN Genova for TOTEM
*/
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "TotemAlignment/T1Alignment/interface/T1InternalAlignment2.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TVector3.h"
#include "TVector2.h"
#include <fstream>
#include <sstream>


#ifndef PI
#define PI 3.14159
#endif
#define _30DEG 0.5236
#define _25DEG 0.4363
#define _60DEG 1.0472
#define _120DEG 2.0944
#define _240DEG 4.1888
#define _300DEG 5.2360
#define _3DEG 0.05239878



//#define _DEBUG_

T1InternalAlignment2::T1InternalAlignment2(const edm::ParameterSet& iConfig):_Verbosity(0)
{

  _Verbosity =  iConfig.getParameter<int>("Verbosity");

}


T1InternalAlignment2::~T1InternalAlignment2()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
T1InternalAlignment2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup )
{


  using namespace edm;



  edm::Handle<T1T2TrackCollection> trackCollection;
  iEvent.getByLabel("t1tracks2","T1TrackColl",trackCollection);


  int TracceRicostruite=0;


#ifdef _DEBUG_
  std::cout << " Taglia della track collection: " << trackCollection->size() <<std::endl;
#endif
  if(_Verbosity>0)
    if (trackCollection->size() == 0){ std::cout << " WARNING: No tracks in T1 in the event !      ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]] " <<std::endl;
    }
    else{
      std::cout << "                                                                " <<trackCollection->size() << " tracks in the event " << std::endl;
    }
  TracceRicostruite=trackCollection->size();



  T1T2TrackCollection::const_iterator TC_it;

  if(TracceRicostruite<10)
  for(TC_it=trackCollection->begin(); TC_it!=trackCollection->end(); TC_it++){

    float x00 = (*TC_it).GetHitT1(0).GlobalPosition().x();
    float y00 = (*TC_it).GetHitT1(0).GlobalPosition().y();
    float z00 = (*TC_it).GetHitT1(0).GlobalPosition().z();


    if(
/*
  fabs((*TC_it).Z_at_Rmin()) < 500 && (*TC_it).Rmin()<50 && (*TC_it).ChiSquared()<40

	  && (fabs(Phi(x00,y00 )-(*TC_it).Phi())/3.14159*180.0) < 7
	  && (*TC_it).GetHitEntries()>=4
*/
0==0
	  ){



// tracce sestante 5+ e 5 hits per traccia
      if( (*TC_it).Eta()>0) {

	if ( sextant(x00,y00,z00) == 0 )
	  allTracks00->push_back((*TC_it));

	if ( sextant(x00,y00,z00) == 1 )
	  allTracks01->push_back((*TC_it));

	if ( sextant(x00,y00,z00) == 2 )
	  allTracks02->push_back((*TC_it));

	if ( sextant(x00,y00,z00) == 3 )
	  allTracks03->push_back((*TC_it));

	if ( sextant(x00,y00,z00) == 4 )
	  allTracks04->push_back((*TC_it));

	if ( sextant(x00,y00,z00) == 5 )
	  allTracks05->push_back((*TC_it));

    }
//negative arm
    if((*TC_it).Eta()<0){
	if ( sextant(x00,y00,z00) == 0 )
	  allTracks10->push_back((*TC_it));

	if ( sextant(x00,y00,z00) == 1 )
	  allTracks11->push_back((*TC_it));

	if ( sextant(x00,y00,z00) == 2 )
	  allTracks12->push_back((*TC_it));

	if ( sextant(x00,y00,z00) == 3 )
	  allTracks13->push_back((*TC_it));

	if ( sextant(x00,y00,z00) == 4 )
	  allTracks14->push_back((*TC_it));

	if ( sextant(x00,y00,z00) == 5 )
	  allTracks15->push_back((*TC_it));



    }
    }
  }







}


// ------------ method called once each job just before starting event loop  ------------
void
T1InternalAlignment2::beginJob()
{
  if(_Verbosity>0)
    cout << "T1InternalAlignment2::beginJob" << endl;




  hDeltaX_CSC0_Pl0 = new TH1D("DeltaX_CSC0_Pl0","DeltaX_CSC0_Pl0",100,-20.,20.);
  hDeltaY_CSC0_Pl0 = new TH1D("DeltaY_CSC0_Pl0","DeltaY_CSC0_Pl0",100,-20.,20.);

  hDeltaX_CSC0_Pl1 = new TH1D("DeltaX_CSC0_Pl1","DeltaX_CSC0_Pl1",100,-20.,20.);
  hDeltaY_CSC0_Pl1 = new TH1D("DeltaY_CSC0_Pl1","DeltaY_CSC0_Pl1",100,-20.,20.);

  hDeltaX_CSC0_Pl2 = new TH1D("DeltaX_CSC0_Pl2","DeltaX_CSC0_Pl2",100,-20.,20.);
  hDeltaY_CSC0_Pl2 = new TH1D("DeltaY_CSC0_Pl2","DeltaY_CSC0_Pl2",100,-20.,20.);

  hDeltaX_CSC0_Pl3 = new TH1D("DeltaX_CSC0_Pl3","DeltaX_CSC0_Pl3",100,-20.,20.);
  hDeltaY_CSC0_Pl3 = new TH1D("DeltaY_CSC0_Pl3","DeltaY_CSC0_Pl3",100,-20.,20.);

  hDeltaX_CSC0_Pl4 = new TH1D("DeltaX_CSC0_Pl4","DeltaX_CSC0_Pl4",100,-20.,20.);
  hDeltaY_CSC0_Pl4 = new TH1D("DeltaY_CSC0_Pl4","DeltaY_CSC0_Pl4",100,-20.,20.);


  hDeltaX_CSC1_Pl0 = new TH1D("DeltaX_CSC1_Pl0","DeltaX_CSC1_Pl0",100,-20.,20.);
  hDeltaY_CSC1_Pl0 = new TH1D("DeltaY_CSC1_Pl0","DeltaY_CSC1_Pl0",100,-20.,20.);

  hDeltaX_CSC1_Pl1 = new TH1D("DeltaX_CSC1_Pl1","DeltaX_CSC1_Pl1",100,-20.,20.);
  hDeltaY_CSC1_Pl1 = new TH1D("DeltaY_CSC1_Pl1","DeltaY_CSC1_Pl1",100,-20.,20.);

  hDeltaX_CSC1_Pl2 = new TH1D("DeltaX_CSC1_Pl2","DeltaX_CSC1_Pl2",100,-20.,20.);
  hDeltaY_CSC1_Pl2 = new TH1D("DeltaY_CSC1_Pl2","DeltaY_CSC1_Pl2",100,-20.,20.);

  hDeltaX_CSC1_Pl3 = new TH1D("DeltaX_CSC1_Pl3","DeltaX_CSC1_Pl3",100,-20.,20.);
  hDeltaY_CSC1_Pl3 = new TH1D("DeltaY_CSC1_Pl3","DeltaY_CSC1_Pl3",100,-20.,20.);

  hDeltaX_CSC1_Pl4 = new TH1D("DeltaX_CSC1_Pl4","DeltaX_CSC1_Pl4",100,-20.,20.);
  hDeltaY_CSC1_Pl4 = new TH1D("DeltaY_CSC1_Pl4","DeltaY_CSC1_Pl4",100,-20.,20.);


  hDeltaX_CSC2_Pl0 = new TH1D("DeltaX_CSC2_Pl0","DeltaX_CSC2_Pl0",100,-20.,20.);
  hDeltaY_CSC2_Pl0 = new TH1D("DeltaY_CSC2_Pl0","DeltaY_CSC2_Pl0",100,-20.,20.);

  hDeltaX_CSC2_Pl1 = new TH1D("DeltaX_CSC2_Pl1","DeltaX_CSC2_Pl1",100,-20.,20.);
  hDeltaY_CSC2_Pl1 = new TH1D("DeltaY_CSC2_Pl1","DeltaY_CSC2_Pl1",100,-20.,20.);

  hDeltaX_CSC2_Pl2 = new TH1D("DeltaX_CSC2_Pl2","DeltaX_CSC2_Pl2",100,-20.,20.);
  hDeltaY_CSC2_Pl2 = new TH1D("DeltaY_CSC2_Pl2","DeltaY_CSC2_Pl2",100,-20.,20.);

  hDeltaX_CSC2_Pl3 = new TH1D("DeltaX_CSC2_Pl3","DeltaX_CSC2_Pl3",100,-20.,20.);
  hDeltaY_CSC2_Pl3 = new TH1D("DeltaY_CSC2_Pl3","DeltaY_CSC2_Pl3",100,-20.,20.);

  hDeltaX_CSC2_Pl4 = new TH1D("DeltaX_CSC2_Pl4","DeltaX_CSC2_Pl4",100,-20.,20.);
  hDeltaY_CSC2_Pl4 = new TH1D("DeltaY_CSC2_Pl4","DeltaY_CSC2_Pl4",100,-20.,20.);


  hDeltaX_CSC3_Pl0 = new TH1D("DeltaX_CSC3_Pl0","DeltaX_CSC3_Pl0",100,-20.,20.);
  hDeltaY_CSC3_Pl0 = new TH1D("DeltaY_CSC3_Pl0","DeltaY_CSC3_Pl0",100,-20.,20.);

  hDeltaX_CSC3_Pl1 = new TH1D("DeltaX_CSC3_Pl1","DeltaX_CSC3_Pl1",100,-20.,20.);
  hDeltaY_CSC3_Pl1 = new TH1D("DeltaY_CSC3_Pl1","DeltaY_CSC3_Pl1",100,-20.,20.);

  hDeltaX_CSC3_Pl2 = new TH1D("DeltaX_CSC3_Pl2","DeltaX_CSC3_Pl2",100,-20.,20.);
  hDeltaY_CSC3_Pl2 = new TH1D("DeltaY_CSC3_Pl2","DeltaY_CSC3_Pl2",100,-20.,20.);

  hDeltaX_CSC3_Pl3 = new TH1D("DeltaX_CSC3_Pl3","DeltaX_CSC3_Pl3",100,-20.,20.);
  hDeltaY_CSC3_Pl3 = new TH1D("DeltaY_CSC3_Pl3","DeltaY_CSC3_Pl3",100,-20.,20.);

  hDeltaX_CSC3_Pl4 = new TH1D("DeltaX_CSC3_Pl4","DeltaX_CSC3_Pl4",100,-20.,20.);
  hDeltaY_CSC3_Pl4 = new TH1D("DeltaY_CSC3_Pl4","DeltaY_CSC3_Pl4",100,-20.,20.);


  hDeltaX_CSC4_Pl0 = new TH1D("DeltaX_CSC4_Pl0","DeltaX_CSC4_Pl0",100,-20.,20.);
  hDeltaY_CSC4_Pl0 = new TH1D("DeltaY_CSC4_Pl0","DeltaY_CSC4_Pl0",100,-20.,20.);

  hDeltaX_CSC4_Pl1 = new TH1D("DeltaX_CSC4_Pl1","DeltaX_CSC4_Pl1",100,-20.,20.);
  hDeltaY_CSC4_Pl1 = new TH1D("DeltaY_CSC4_Pl1","DeltaY_CSC4_Pl1",100,-20.,20.);

  hDeltaX_CSC4_Pl2 = new TH1D("DeltaX_CSC4_Pl2","DeltaX_CSC4_Pl2",100,-20.,20.);
  hDeltaY_CSC4_Pl2 = new TH1D("DeltaY_CSC4_Pl2","DeltaY_CSC4_Pl2",100,-20.,20.);

  hDeltaX_CSC4_Pl3 = new TH1D("DeltaX_CSC4_Pl3","DeltaX_CSC4_Pl3",100,-20.,20.);
  hDeltaY_CSC4_Pl3 = new TH1D("DeltaY_CSC4_Pl3","DeltaY_CSC4_Pl3",100,-20.,20.);

  hDeltaX_CSC4_Pl4 = new TH1D("DeltaX_CSC4_Pl4","DeltaX_CSC4_Pl4",100,-20.,20.);
  hDeltaY_CSC4_Pl4 = new TH1D("DeltaY_CSC4_Pl4","DeltaY_CSC4_Pl4",100,-20.,20.);


  hDeltaX_CSC5_Pl0 = new TH1D("DeltaX_CSC5_Pl0","DeltaX_CSC5_Pl0",100,-20.,20.);
  hDeltaY_CSC5_Pl0 = new TH1D("DeltaY_CSC5_Pl0","DeltaY_CSC5_Pl0",100,-20.,20.);

  hDeltaX_CSC5_Pl1 = new TH1D("DeltaX_CSC5_Pl1","DeltaX_CSC5_Pl1",100,-20.,20.);
  hDeltaY_CSC5_Pl1 = new TH1D("DeltaY_CSC5_Pl1","DeltaY_CSC5_Pl1",100,-20.,20.);

  hDeltaX_CSC5_Pl2 = new TH1D("DeltaX_CSC5_Pl2","DeltaX_CSC5_Pl2",100,-20.,20.);
  hDeltaY_CSC5_Pl2 = new TH1D("DeltaY_CSC5_Pl2","DeltaY_CSC5_Pl2",100,-20.,20.);

  hDeltaX_CSC5_Pl3 = new TH1D("DeltaX_CSC5_Pl3","DeltaX_CSC5_Pl3",100,-20.,20.);
  hDeltaY_CSC5_Pl3 = new TH1D("DeltaY_CSC5_Pl3","DeltaY_CSC5_Pl3",100,-20.,20.);

  hDeltaX_CSC5_Pl4 = new TH1D("DeltaX_CSC5_Pl4","DeltaX_CSC5_Pl4",100,-20.,20.);
  hDeltaY_CSC5_Pl4 = new TH1D("DeltaY_CSC5_Pl4","DeltaY_CSC5_Pl4",100,-20.,20.);





//minus


  hDeltaMENOX_CSC0_Pl0 = new TH1D("DeltaMENOX_CSC0_Pl0","DeltaMENOX_CSC0_Pl0",100,-20.,20.);
  hDeltaMENOY_CSC0_Pl0 = new TH1D("DeltaMENOY_CSC0_Pl0","DeltaMENOY_CSC0_Pl0",100,-20.,20.);

  hDeltaMENOX_CSC0_Pl1 = new TH1D("DeltaMENOX_CSC0_Pl1","DeltaMENOX_CSC0_Pl1",100,-20.,20.);
  hDeltaMENOY_CSC0_Pl1 = new TH1D("DeltaMENOY_CSC0_Pl1","DeltaMENOY_CSC0_Pl1",100,-20.,20.);

  hDeltaMENOX_CSC0_Pl2 = new TH1D("DeltaMENOX_CSC0_Pl2","DeltaMENOX_CSC0_Pl2",100,-20.,20.);
  hDeltaMENOY_CSC0_Pl2 = new TH1D("DeltaMENOY_CSC0_Pl2","DeltaMENOY_CSC0_Pl2",100,-20.,20.);

  hDeltaMENOX_CSC0_Pl3 = new TH1D("DeltaMENOX_CSC0_Pl3","DeltaMENOX_CSC0_Pl3",100,-20.,20.);
  hDeltaMENOY_CSC0_Pl3 = new TH1D("DeltaMENOY_CSC0_Pl3","DeltaMENOY_CSC0_Pl3",100,-20.,20.);

  hDeltaMENOX_CSC0_Pl4 = new TH1D("DeltaMENOX_CSC0_Pl4","DeltaMENOX_CSC0_Pl4",100,-20.,20.);
  hDeltaMENOY_CSC0_Pl4 = new TH1D("DeltaMENOY_CSC0_Pl4","DeltaMENOY_CSC0_Pl4",100,-20.,20.);


  hDeltaMENOX_CSC1_Pl0 = new TH1D("DeltaMENOX_CSC1_Pl0","DeltaMENOX_CSC1_Pl0",100,-20.,20.);
  hDeltaMENOY_CSC1_Pl0 = new TH1D("DeltaMENOY_CSC1_Pl0","DeltaMENOY_CSC1_Pl0",100,-20.,20.);

  hDeltaMENOX_CSC1_Pl1 = new TH1D("DeltaMENOX_CSC1_Pl1","DeltaMENOX_CSC1_Pl1",100,-20.,20.);
  hDeltaMENOY_CSC1_Pl1 = new TH1D("DeltaMENOY_CSC1_Pl1","DeltaMENOY_CSC1_Pl1",100,-20.,20.);

  hDeltaMENOX_CSC1_Pl2 = new TH1D("DeltaMENOX_CSC1_Pl2","DeltaMENOX_CSC1_Pl2",100,-20.,20.);
  hDeltaMENOY_CSC1_Pl2 = new TH1D("DeltaMENOY_CSC1_Pl2","DeltaMENOY_CSC1_Pl2",100,-20.,20.);

  hDeltaMENOX_CSC1_Pl3 = new TH1D("DeltaMENOX_CSC1_Pl3","DeltaMENOX_CSC1_Pl3",100,-20.,20.);
  hDeltaMENOY_CSC1_Pl3 = new TH1D("DeltaMENOY_CSC1_Pl3","DeltaMENOY_CSC1_Pl3",100,-20.,20.);

  hDeltaMENOX_CSC1_Pl4 = new TH1D("DeltaMENOX_CSC1_Pl4","DeltaMENOX_CSC1_Pl4",100,-20.,20.);
  hDeltaMENOY_CSC1_Pl4 = new TH1D("DeltaMENOY_CSC1_Pl4","DeltaMENOY_CSC1_Pl4",100,-20.,20.);


  hDeltaMENOX_CSC2_Pl0 = new TH1D("DeltaMENOX_CSC2_Pl0","DeltaMENOX_CSC2_Pl0",100,-20.,20.);
  hDeltaMENOY_CSC2_Pl0 = new TH1D("DeltaMENOY_CSC2_Pl0","DeltaMENOY_CSC2_Pl0",100,-20.,20.);

  hDeltaMENOX_CSC2_Pl1 = new TH1D("DeltaMENOX_CSC2_Pl1","DeltaMENOX_CSC2_Pl1",100,-20.,20.);
  hDeltaMENOY_CSC2_Pl1 = new TH1D("DeltaMENOY_CSC2_Pl1","DeltaMENOY_CSC2_Pl1",100,-20.,20.);

  hDeltaMENOX_CSC2_Pl2 = new TH1D("DeltaMENOX_CSC2_Pl2","DeltaMENOX_CSC2_Pl2",100,-20.,20.);
  hDeltaMENOY_CSC2_Pl2 = new TH1D("DeltaMENOY_CSC2_Pl2","DeltaMENOY_CSC2_Pl2",100,-20.,20.);

  hDeltaMENOX_CSC2_Pl3 = new TH1D("DeltaMENOX_CSC2_Pl3","DeltaMENOX_CSC2_Pl3",100,-20.,20.);
  hDeltaMENOY_CSC2_Pl3 = new TH1D("DeltaMENOY_CSC2_Pl3","DeltaMENOY_CSC2_Pl3",100,-20.,20.);

  hDeltaMENOX_CSC2_Pl4 = new TH1D("DeltaMENOX_CSC2_Pl4","DeltaMENOX_CSC2_Pl4",100,-20.,20.);
  hDeltaMENOY_CSC2_Pl4 = new TH1D("DeltaMENOY_CSC2_Pl4","DeltaMENOY_CSC2_Pl4",100,-20.,20.);


  hDeltaMENOX_CSC3_Pl0 = new TH1D("DeltaMENOX_CSC3_Pl0","DeltaMENOX_CSC3_Pl0",100,-20.,20.);
  hDeltaMENOY_CSC3_Pl0 = new TH1D("DeltaMENOY_CSC3_Pl0","DeltaMENOY_CSC3_Pl0",100,-20.,20.);

  hDeltaMENOX_CSC3_Pl1 = new TH1D("DeltaMENOX_CSC3_Pl1","DeltaMENOX_CSC3_Pl1",100,-20.,20.);
  hDeltaMENOY_CSC3_Pl1 = new TH1D("DeltaMENOY_CSC3_Pl1","DeltaMENOY_CSC3_Pl1",100,-20.,20.);

  hDeltaMENOX_CSC3_Pl2 = new TH1D("DeltaMENOX_CSC3_Pl2","DeltaMENOX_CSC3_Pl2",100,-20.,20.);
  hDeltaMENOY_CSC3_Pl2 = new TH1D("DeltaMENOY_CSC3_Pl2","DeltaMENOY_CSC3_Pl2",100,-20.,20.);

  hDeltaMENOX_CSC3_Pl3 = new TH1D("DeltaMENOX_CSC3_Pl3","DeltaMENOX_CSC3_Pl3",100,-20.,20.);
  hDeltaMENOY_CSC3_Pl3 = new TH1D("DeltaMENOY_CSC3_Pl3","DeltaMENOY_CSC3_Pl3",100,-20.,20.);

  hDeltaMENOX_CSC3_Pl4 = new TH1D("DeltaMENOX_CSC3_Pl4","DeltaMENOX_CSC3_Pl4",100,-20.,20.);
  hDeltaMENOY_CSC3_Pl4 = new TH1D("DeltaMENOY_CSC3_Pl4","DeltaMENOY_CSC3_Pl4",100,-20.,20.);


  hDeltaMENOX_CSC4_Pl0 = new TH1D("DeltaMENOX_CSC4_Pl0","DeltaMENOX_CSC4_Pl0",100,-20.,20.);
  hDeltaMENOY_CSC4_Pl0 = new TH1D("DeltaMENOY_CSC4_Pl0","DeltaMENOY_CSC4_Pl0",100,-20.,20.);

  hDeltaMENOX_CSC4_Pl1 = new TH1D("DeltaMENOX_CSC4_Pl1","DeltaMENOX_CSC4_Pl1",100,-20.,20.);
  hDeltaMENOY_CSC4_Pl1 = new TH1D("DeltaMENOY_CSC4_Pl1","DeltaMENOY_CSC4_Pl1",100,-20.,20.);

  hDeltaMENOX_CSC4_Pl2 = new TH1D("DeltaMENOX_CSC4_Pl2","DeltaMENOX_CSC4_Pl2",100,-20.,20.);
  hDeltaMENOY_CSC4_Pl2 = new TH1D("DeltaMENOY_CSC4_Pl2","DeltaMENOY_CSC4_Pl2",100,-20.,20.);

  hDeltaMENOX_CSC4_Pl3 = new TH1D("DeltaMENOX_CSC4_Pl3","DeltaMENOX_CSC4_Pl3",100,-20.,20.);
  hDeltaMENOY_CSC4_Pl3 = new TH1D("DeltaMENOY_CSC4_Pl3","DeltaMENOY_CSC4_Pl3",100,-20.,20.);

  hDeltaMENOX_CSC4_Pl4 = new TH1D("DeltaMENOX_CSC4_Pl4","DeltaMENOX_CSC4_Pl4",100,-20.,20.);
  hDeltaMENOY_CSC4_Pl4 = new TH1D("DeltaMENOY_CSC4_Pl4","DeltaMENOY_CSC4_Pl4",100,-20.,20.);


  hDeltaMENOX_CSC5_Pl0 = new TH1D("DeltaMENOX_CSC5_Pl0","DeltaMENOX_CSC5_Pl0",100,-20.,20.);
  hDeltaMENOY_CSC5_Pl0 = new TH1D("DeltaMENOY_CSC5_Pl0","DeltaMENOY_CSC5_Pl0",100,-20.,20.);

  hDeltaMENOX_CSC5_Pl1 = new TH1D("DeltaMENOX_CSC5_Pl1","DeltaMENOX_CSC5_Pl1",100,-20.,20.);
  hDeltaMENOY_CSC5_Pl1 = new TH1D("DeltaMENOY_CSC5_Pl1","DeltaMENOY_CSC5_Pl1",100,-20.,20.);

  hDeltaMENOX_CSC5_Pl2 = new TH1D("DeltaMENOX_CSC5_Pl2","DeltaMENOX_CSC5_Pl2",100,-20.,20.);
  hDeltaMENOY_CSC5_Pl2 = new TH1D("DeltaMENOY_CSC5_Pl2","DeltaMENOY_CSC5_Pl2",100,-20.,20.);

  hDeltaMENOX_CSC5_Pl3 = new TH1D("DeltaMENOX_CSC5_Pl3","DeltaMENOX_CSC5_Pl3",100,-20.,20.);
  hDeltaMENOY_CSC5_Pl3 = new TH1D("DeltaMENOY_CSC5_Pl3","DeltaMENOY_CSC5_Pl3",100,-20.,20.);

  hDeltaMENOX_CSC5_Pl4 = new TH1D("DeltaMENOX_CSC5_Pl4","DeltaMENOX_CSC5_Pl4",100,-20.,20.);
  hDeltaMENOY_CSC5_Pl4 = new TH1D("DeltaMENOY_CSC5_Pl4","DeltaMENOY_CSC5_Pl4",100,-20.,20.);











  allTracks00 = new T1T2TrackCollection();
  allTracks01 = new T1T2TrackCollection();
  allTracks02 = new T1T2TrackCollection();
  allTracks03 = new T1T2TrackCollection();
  allTracks04 = new T1T2TrackCollection();
  allTracks05 = new T1T2TrackCollection();

  allTracks10 = new T1T2TrackCollection();
  allTracks11 = new T1T2TrackCollection();
  allTracks12 = new T1T2TrackCollection();
  allTracks13 = new T1T2TrackCollection();
  allTracks14 = new T1T2TrackCollection();
  allTracks15 = new T1T2TrackCollection();


}

// ------------ method called once each job just after ending the event loop  ------------
void T1InternalAlignment2::endJob() {
  if(_Verbosity>0)
    cout << "T1InternalAlignment2::endJob" << endl;


  theFile = new TFile("analysisFile.root","RECREATE");


// leggo precedente file dell'allineamento
  double Align_x[2][5][6];
  double Align_y[2][5][6];
  double Align_theta[2][5][6];
  for(int i=0; i<2; i++)
    for(int j=0; j<5; j++)
      for(int k=0; k<6; k++){
	Align_x[i][j][k]=0;
	Align_y[i][j][k]=0;
	Align_theta[i][j][k]=0;
      }
  unsigned int a,p,c;
  float x_,y_,th_;

  FILE *fileA = fopen("T1Align.txt","r");
  if(fileA!=NULL){
    while(!feof(fileA)){

      fscanf(fileA,"%u%u%u%f%f%f",&a,&p,&c,&x_,&y_,&th_);

      Align_x[a][p][c]=x_;
      Align_y[a][p][c]=y_;
      Align_theta[a][p][c]=th_;
    }
    fclose(fileA);
  }
//--------------------------------------------------


  double temp_Align_x[2][5][6];
  double temp_Align_y[2][5][6];
  double temp_Align_theta[2][5][6];

  int iret00 = fitParam(0,0,allTracks00, temp_Align_x, temp_Align_y, temp_Align_theta);
  int iret01 = fitParam(0,1,allTracks01, temp_Align_x, temp_Align_y, temp_Align_theta);
  int iret02 = fitParam(0,2,allTracks02, temp_Align_x, temp_Align_y, temp_Align_theta);
  int iret03 = fitParam(0,3,allTracks03, temp_Align_x, temp_Align_y, temp_Align_theta);
  int iret04 = fitParam(0,4,allTracks04, temp_Align_x, temp_Align_y, temp_Align_theta);
  int iret05 = fitParam(0,5,allTracks05, temp_Align_x, temp_Align_y, temp_Align_theta);
  int iret10 = fitParam(1,0,allTracks10, temp_Align_x, temp_Align_y, temp_Align_theta);
  int iret11 = fitParam(1,1,allTracks11, temp_Align_x, temp_Align_y, temp_Align_theta);
  int iret12 = fitParam(1,2,allTracks12, temp_Align_x, temp_Align_y, temp_Align_theta);
  int iret13 = fitParam(1,3,allTracks13, temp_Align_x, temp_Align_y, temp_Align_theta);
  int iret14 = fitParam(1,4,allTracks14, temp_Align_x, temp_Align_y, temp_Align_theta);
  int iret15 = fitParam(1,5,allTracks15, temp_Align_x, temp_Align_y, temp_Align_theta);

  if (iret00 != 0 || iret01 != 0 || iret02 != 0 || iret03 != 0 || iret04 != 0 || iret05 != 0 || iret10 != 0 || iret11 != 0 || iret12 != 0 || iret13 != 0 || iret14 != 0 || iret15 != 0) {
    std::cout << "T1RecHitProducer:  MINIMIZATION FAILED! " << std::endl;

  }

  else{

// salvo i risultati dell'allineamento

// forcing not to align the following CSC
/*
    temp_Align_x[0][3][0]=0; temp_Align_y[0][3][0]=0;
    temp_Align_x[0][4][0]=0; temp_Align_y[0][4][0]=0;
    temp_Align_x[0][3][1]=0; temp_Align_y[0][3][1]=0;
    temp_Align_x[0][1][1]=0; temp_Align_y[0][1][1]=0;
    temp_Align_x[0][4][2]=0; temp_Align_y[0][4][2]=0;
    temp_Align_x[0][3][4]=0; temp_Align_y[0][3][4]=0;
*/


    FILE *file = fopen("fit_results.txt","a");
    FILE *fileB = fopen("NewT1Align.txt","w");
    for(unsigned int i=0; i<2; i++)
      for(unsigned int j=0; j<5; j++)
	for(unsigned int k=0; k<6; k++){
	  fprintf(file,"%hd %hd %hd", i,j,k);
	  fprintf(fileB,"%hd %hd %hd", i,j,k);

	  fprintf(file," %g %g %g\n",temp_Align_x[i][j][k],temp_Align_y[i][j][k], temp_Align_theta[i][j][k]);
	  fprintf(fileB," %g %g %g\n",temp_Align_x[i][j][k]+Align_x[i][j][k],temp_Align_y[i][j][k]+Align_y[i][j][k],temp_Align_theta[i][j][k]+Align_theta[i][j][k]);

	}

    fclose(fileB);
    fclose(file);

  }











  delete allTracks00;
  delete allTracks01;
  delete allTracks02;
  delete allTracks03;
  delete allTracks04;
  delete allTracks05;

  delete allTracks10;
  delete allTracks11;
  delete allTracks12;
  delete allTracks13;
  delete allTracks14;
  delete allTracks15;

// theFile = new TFile("tracksFile.root","RECREATE");



  hDeltaX_CSC0_Pl0 -> Write();
  hDeltaY_CSC0_Pl0 -> Write();

  hDeltaX_CSC0_Pl1 -> Write();
  hDeltaY_CSC0_Pl1 -> Write();

  hDeltaX_CSC0_Pl2 -> Write();
  hDeltaY_CSC0_Pl2 -> Write();

  hDeltaX_CSC0_Pl3 -> Write();
  hDeltaY_CSC0_Pl3 -> Write();

  hDeltaX_CSC0_Pl4 -> Write();
  hDeltaY_CSC0_Pl4 -> Write();





  hDeltaX_CSC1_Pl0 -> Write();
  hDeltaY_CSC1_Pl0 -> Write();

  hDeltaX_CSC1_Pl1 -> Write();
  hDeltaY_CSC1_Pl1 -> Write();

  hDeltaX_CSC1_Pl2 -> Write();
  hDeltaY_CSC1_Pl2 -> Write();

  hDeltaX_CSC1_Pl3 -> Write();
  hDeltaY_CSC1_Pl3 -> Write();

  hDeltaX_CSC1_Pl4 -> Write();
  hDeltaY_CSC1_Pl4 -> Write();





  hDeltaX_CSC2_Pl0 -> Write();
  hDeltaY_CSC2_Pl0 -> Write();

  hDeltaX_CSC2_Pl1 -> Write();
  hDeltaY_CSC2_Pl1 -> Write();

  hDeltaX_CSC2_Pl2 -> Write();
  hDeltaY_CSC2_Pl2 -> Write();

  hDeltaX_CSC2_Pl3 -> Write();
  hDeltaY_CSC2_Pl3 -> Write();

  hDeltaX_CSC2_Pl4 -> Write();
  hDeltaY_CSC2_Pl4 -> Write();





  hDeltaX_CSC3_Pl0 -> Write();
  hDeltaY_CSC3_Pl0 -> Write();

  hDeltaX_CSC3_Pl1 -> Write();
  hDeltaY_CSC3_Pl1 -> Write();

  hDeltaX_CSC3_Pl2 -> Write();
  hDeltaY_CSC3_Pl2 -> Write();

  hDeltaX_CSC3_Pl3 -> Write();
  hDeltaY_CSC3_Pl3 -> Write();

  hDeltaX_CSC3_Pl4 -> Write();
  hDeltaY_CSC3_Pl4 -> Write();





  hDeltaX_CSC4_Pl0 -> Write();
  hDeltaY_CSC4_Pl0 -> Write();

  hDeltaX_CSC4_Pl1 -> Write();
  hDeltaY_CSC4_Pl1 -> Write();

  hDeltaX_CSC4_Pl2 -> Write();
  hDeltaY_CSC4_Pl2 -> Write();

  hDeltaX_CSC4_Pl3 -> Write();
  hDeltaY_CSC4_Pl3 -> Write();

  hDeltaX_CSC4_Pl4 -> Write();
  hDeltaY_CSC4_Pl4 -> Write();





  hDeltaX_CSC5_Pl0 -> Write();
  hDeltaY_CSC5_Pl0 -> Write();

  hDeltaX_CSC5_Pl1 -> Write();
  hDeltaY_CSC5_Pl1 -> Write();

  hDeltaX_CSC5_Pl2 -> Write();
  hDeltaY_CSC5_Pl2 -> Write();

  hDeltaX_CSC5_Pl3 -> Write();
  hDeltaY_CSC5_Pl3 -> Write();

  hDeltaX_CSC5_Pl4 -> Write();
  hDeltaY_CSC5_Pl4 -> Write();


// minus

  hDeltaMENOX_CSC0_Pl0 -> Write();
  hDeltaMENOY_CSC0_Pl0 -> Write();

  hDeltaMENOX_CSC0_Pl1 -> Write();
  hDeltaMENOY_CSC0_Pl1 -> Write();

  hDeltaMENOX_CSC0_Pl2 -> Write();
  hDeltaMENOY_CSC0_Pl2 -> Write();

  hDeltaMENOX_CSC0_Pl3 -> Write();
  hDeltaMENOY_CSC0_Pl3 -> Write();

  hDeltaMENOX_CSC0_Pl4 -> Write();
  hDeltaMENOY_CSC0_Pl4 -> Write();





  hDeltaMENOX_CSC1_Pl0 -> Write();
  hDeltaMENOY_CSC1_Pl0 -> Write();

  hDeltaMENOX_CSC1_Pl1 -> Write();
  hDeltaMENOY_CSC1_Pl1 -> Write();

  hDeltaMENOX_CSC1_Pl2 -> Write();
  hDeltaMENOY_CSC1_Pl2 -> Write();

  hDeltaMENOX_CSC1_Pl3 -> Write();
  hDeltaMENOY_CSC1_Pl3 -> Write();

  hDeltaMENOX_CSC1_Pl4 -> Write();
  hDeltaMENOY_CSC1_Pl4 -> Write();





  hDeltaMENOX_CSC2_Pl0 -> Write();
  hDeltaMENOY_CSC2_Pl0 -> Write();

  hDeltaMENOX_CSC2_Pl1 -> Write();
  hDeltaMENOY_CSC2_Pl1 -> Write();

  hDeltaMENOX_CSC2_Pl2 -> Write();
  hDeltaMENOY_CSC2_Pl2 -> Write();

  hDeltaMENOX_CSC2_Pl3 -> Write();
  hDeltaMENOY_CSC2_Pl3 -> Write();

  hDeltaMENOX_CSC2_Pl4 -> Write();
  hDeltaMENOY_CSC2_Pl4 -> Write();





  hDeltaMENOX_CSC3_Pl0 -> Write();
  hDeltaMENOY_CSC3_Pl0 -> Write();

  hDeltaMENOX_CSC3_Pl1 -> Write();
  hDeltaMENOY_CSC3_Pl1 -> Write();

  hDeltaMENOX_CSC3_Pl2 -> Write();
  hDeltaMENOY_CSC3_Pl2 -> Write();

  hDeltaMENOX_CSC3_Pl3 -> Write();
  hDeltaMENOY_CSC3_Pl3 -> Write();

  hDeltaMENOX_CSC3_Pl4 -> Write();
  hDeltaMENOY_CSC3_Pl4 -> Write();





  hDeltaMENOX_CSC4_Pl0 -> Write();
  hDeltaMENOY_CSC4_Pl0 -> Write();

  hDeltaMENOX_CSC4_Pl1 -> Write();
  hDeltaMENOY_CSC4_Pl1 -> Write();

  hDeltaMENOX_CSC4_Pl2 -> Write();
  hDeltaMENOY_CSC4_Pl2 -> Write();

  hDeltaMENOX_CSC4_Pl3 -> Write();
  hDeltaMENOY_CSC4_Pl3 -> Write();

  hDeltaMENOX_CSC4_Pl4 -> Write();
  hDeltaMENOY_CSC4_Pl4 -> Write();





  hDeltaMENOX_CSC5_Pl0 -> Write();
  hDeltaMENOY_CSC5_Pl0 -> Write();

  hDeltaMENOX_CSC5_Pl1 -> Write();
  hDeltaMENOY_CSC5_Pl1 -> Write();

  hDeltaMENOX_CSC5_Pl2 -> Write();
  hDeltaMENOY_CSC5_Pl2 -> Write();

  hDeltaMENOX_CSC5_Pl3 -> Write();
  hDeltaMENOY_CSC5_Pl3 -> Write();

  hDeltaMENOX_CSC5_Pl4 -> Write();
  hDeltaMENOY_CSC5_Pl4 -> Write();








  theFile->Close();
}


int T1InternalAlignment2::fitParam(int braccio, int sestante, T1T2TrackCollection *tracce, double allineo_x[2][5][6], double allineo_y[2][5][6], double allineo_theta[2][5][6]){

  enum {maxTracks = 100000};
  const unsigned int minTracks = 100;

  assert( (braccio == 0 || braccio == 1) && (sestante == 0 || sestante == 1 || sestante == 2 || sestante == 3 || sestante == 4 || sestante == 5 ));


  if(tracce->size()<minTracks){
    if(braccio == 0)
      std::cout << "WARNING: less than " << minTracks << " tracks in arm PLUS, sextant " << sestante << ". Alignment of this sextant aborted." << std::endl;
    if(braccio == 1)
      std::cout << "WARNING: less than " << minTracks << " tracks in arm MINUS, sextant "<<sestante << ". Alignment of this sextant aborted." << std::endl;

    for(unsigned int j=0; j<5; j++)
      {
	allineo_x[braccio][j][sestante]=0;
	allineo_y[braccio][j][sestante]=0;
	allineo_theta[braccio][j][sestante]=0;
      }
    return 0;
  }

  std::cout << "NOTE:  Tracks # limited to "<< maxTracks << std::endl;

  if(braccio == 0){
    std::cout << "Fitting sextant " << sestante << " in arm PLUS with " << tracce->size() << " tracks" <<std::endl;
  }else{
    std::cout << "Fitting sextant " << sestante << " in arm MINUS with " << tracce->size() << " tracks" << std::endl;
  }



  T1T2TrackCollection::const_iterator TC_it;





  double *eeX=NULL;
  double *eeY=NULL;
  double *x=NULL;
  double *y=NULL;
  double *z=NULL;
  double *x_prev=NULL;
  double *y_prev=NULL;

  x=new double[5*maxTracks];
  y=new double[5*maxTracks];
  z=new double[5*maxTracks];
  x_prev=new double[5*maxTracks];
  y_prev=new double[5*maxTracks];
  eeX=new double[5*maxTracks];
  eeY=new double[5*maxTracks];


  vector<TVector3> _Pl0;
  vector<TVector2> _Pl0_prev;
  vector<TVector3> _Pl1;
  vector<TVector2> _Pl1_prev;
  vector<TVector3> _Pl2;
  vector<TVector2> _Pl2_prev;
  vector<TVector3> _Pl3;
  vector<TVector2> _Pl3_prev;
  vector<TVector3> _Pl4;
  vector<TVector2> _Pl4_prev;
  vector<TVector3> _all;
  vector<TVector2> _all_prev;
  vector<TVector2> _all_error;



  int counter = 0;
  for(TC_it=tracce->begin(); TC_it!=tracce->end(); TC_it++){
    for(unsigned int iii = 0; iii<(*TC_it).GetHitEntries(); iii++){
      if(counter<(5*maxTracks)){
	x[counter]=(*TC_it).GetHitT1(iii).GlobalPosition().x();
	y[counter]=(*TC_it).GetHitT1(iii).GlobalPosition().y();
	z[counter]=(*TC_it).GetHitT1(iii).GlobalPosition().z();
	x_prev[counter]=(*TC_it).GetTx() * z[counter] + (*TC_it).X0();
	y_prev[counter]=(*TC_it).GetTy() * z[counter] + (*TC_it).Y0();
	  eeX[counter]=(*TC_it).GetTxSigma()*fabs(z[counter])+(*TC_it).X0Sigma();
	  eeY[counter]=(*TC_it).GetTySigma()*fabs(z[counter])+(*TC_it).Y0Sigma();
	if(_Verbosity>0)
	  cout << x[counter] << " " << x_prev[counter] << "   " << y[counter] << " " << y_prev[counter] << "   " << z[counter] <<endl;

	int sex = sextant(x[counter],y[counter],z[counter]);

	if(fabs(z[counter])<8000){

	  if(z[counter]>0)
	  switch(sex){
	  case 0:
	    hDeltaX_CSC0_Pl0 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC0_Pl0 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 1:
	    hDeltaX_CSC1_Pl0 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC1_Pl0 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 2:
	    hDeltaX_CSC2_Pl0 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC2_Pl0 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 3:
	    hDeltaX_CSC3_Pl0 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC3_Pl0 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 4:
	    hDeltaX_CSC4_Pl0 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC4_Pl0 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 5:
	    hDeltaX_CSC5_Pl0 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC5_Pl0 -> Fill(y[counter]-y_prev[counter]);;
	    //  cout << "0: " << Eta(x[counter],y[counter],z[counter]) << " " << Phi(x[counter],y[counter]) << endl;
	    break;

	  }
	  if(z[counter]<0)
	  switch(sex){
	  case 0:
	    hDeltaMENOX_CSC0_Pl0 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC0_Pl0 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 1:
	    hDeltaMENOX_CSC1_Pl0 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC1_Pl0 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 2:
	    hDeltaMENOX_CSC2_Pl0 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC2_Pl0 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 3:
	    hDeltaMENOX_CSC3_Pl0 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC3_Pl0 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 4:
	    hDeltaMENOX_CSC4_Pl0 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC4_Pl0 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 5:
	    hDeltaMENOX_CSC5_Pl0 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC5_Pl0 -> Fill(y[counter]-y_prev[counter]);;
	    //  cout << "0: " << Eta(x[counter],y[counter],z[counter]) << " " << Phi(x[counter],y[counter]) << endl;
	    break;

	  }






	  TVector3 V3(x[counter],y[counter],z[counter]);
	  TVector2 V2(x_prev[counter],y_prev[counter]);
	  TVector2 eV2(eeX[counter],eeY[counter]);


	  _Pl0.push_back(V3);
	  _Pl0_prev.push_back(V2);
	  _all.push_back(V3);
	  _all_prev.push_back(V2);
	  _all_error.push_back(eV2);

	}
	if(fabs(z[counter])<8400 && fabs(z[counter])>8000){
	  if(z[counter]>0)
	  switch(sex){
	  case 0:
	    hDeltaX_CSC0_Pl1 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC0_Pl1 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 1:
	    hDeltaX_CSC1_Pl1 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC1_Pl1 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 2:
	    hDeltaX_CSC2_Pl1 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC2_Pl1 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 3:
	    hDeltaX_CSC3_Pl1 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC3_Pl1 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 4:
	    hDeltaX_CSC4_Pl1 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC4_Pl1 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 5:
	    hDeltaX_CSC5_Pl1 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC5_Pl1 -> Fill(y[counter]-y_prev[counter]);;
	      //  cout <<  "1: " << Eta(x[counter],y[counter],z[counter]) << " " << Phi(x[counter],y[counter]) << endl;
	    break;

	  }
	  if(z[counter]<0)
	  switch(sex){
	  case 0:
	    hDeltaMENOX_CSC0_Pl1 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC0_Pl1 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 1:
	    hDeltaMENOX_CSC1_Pl1 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC1_Pl1 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 2:
	    hDeltaMENOX_CSC2_Pl1 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC2_Pl1 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 3:
	    hDeltaMENOX_CSC3_Pl1 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC3_Pl1 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 4:
	    hDeltaMENOX_CSC4_Pl1 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC4_Pl1 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 5:
	    hDeltaMENOX_CSC5_Pl1 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC5_Pl1 -> Fill(y[counter]-y_prev[counter]);;
	      //  cout <<  "1: " << Eta(x[counter],y[counter],z[counter]) << " " << Phi(x[counter],y[counter]) << endl;
	    break;

	  }



	  TVector3 V3(x[counter],y[counter],z[counter]);
	  TVector2 V2(x_prev[counter],y_prev[counter]);
	  TVector2 eV2(eeX[counter],eeY[counter]);



	  _all_error.push_back(eV2);

	  _Pl1.push_back(V3);
	  _Pl1_prev.push_back(V2);
	  _all.push_back(V3);
	  _all_prev.push_back(V2);
	  _all_error.push_back(eV2);

 	}
	if(fabs(z[counter])<9000 && fabs(z[counter])>8400){
	  if(z[counter]>0)
	  switch(sex){
	  case 0:
	    hDeltaX_CSC0_Pl2 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC0_Pl2 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 1:
	    hDeltaX_CSC1_Pl2 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC1_Pl2 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 2:
	    hDeltaX_CSC2_Pl2 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC2_Pl2 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 3:
	    hDeltaX_CSC3_Pl2 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC3_Pl2 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 4:
	    hDeltaX_CSC4_Pl2 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC4_Pl2 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 5:
	    hDeltaX_CSC5_Pl2 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC5_Pl2 -> Fill(y[counter]-y_prev[counter]);;
	      //  cout <<  "2: " << Eta(x[counter],y[counter],z[counter]) << " " << Phi(x[counter],y[counter]) << endl;
	    break;

	  }
	  if(z[counter]<0)
	  switch(sex){
	  case 0:
	    hDeltaMENOX_CSC0_Pl2 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC0_Pl2 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 1:
	    hDeltaMENOX_CSC1_Pl2 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC1_Pl2 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 2:
	    hDeltaMENOX_CSC2_Pl2 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC2_Pl2 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 3:
	    hDeltaMENOX_CSC3_Pl2 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC3_Pl2 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 4:
	    hDeltaMENOX_CSC4_Pl2 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC4_Pl2 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 5:
	    hDeltaMENOX_CSC5_Pl2 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC5_Pl2 -> Fill(y[counter]-y_prev[counter]);;
	      //  cout <<  "2: " << Eta(x[counter],y[counter],z[counter]) << " " << Phi(x[counter],y[counter]) << endl;
	    break;

	  }



	  TVector3 V3(x[counter],y[counter],z[counter]);
	  TVector2 V2(x_prev[counter],y_prev[counter]);
	  TVector2 eV2(eeX[counter],eeY[counter]);


	  _Pl2.push_back(V3);
	  _Pl2_prev.push_back(V2);
	  _all.push_back(V3);
	  _all_prev.push_back(V2);
	  _all_error.push_back(eV2);


 	}
	if(fabs(z[counter])<9500 && fabs(z[counter])>9000){
	  if(z[counter]>0)
	  switch(sex){
	  case 0:
	    hDeltaX_CSC0_Pl3 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC0_Pl3 -> Fill(y[counter]-y_prev[counter]);;
/*
	    if( (x[counter]-x_prev[counter]) > 1.5 ){
	      cout << "BAD:" <<endl;
	      cout << x[counter] << " " << x_prev[counter] << endl;
	      cout << y[counter] << " " << y_prev[counter] << endl;
	      cout << Eta(x[counter],y[counter],z[counter]) << " " << Phi(x[counter],y[counter]) << endl;
	      cout << (*TC_it).GetHitEntries() << " " << (*TC_it).ChiSquared() << endl;

	    }
	    if( fabs(x[counter]-x_prev[counter]) < 1. ){
	      cout << "GOOD:" <<endl;
	      cout << x[counter] << " " << x_prev[counter] << endl;
	      cout << y[counter] << " " << y_prev[counter] << endl;
	      cout << Eta(x[counter],y[counter],z[counter]) << " " << Phi(x[counter],y[counter]) << endl;
	      cout << (*TC_it).GetHitEntries() << " " << (*TC_it).ChiSquared() << endl;

	    }
*/

	    break;
	  case 1:
	    hDeltaX_CSC1_Pl3 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC1_Pl3 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 2:
	    hDeltaX_CSC2_Pl3 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC2_Pl3 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 3:
	    hDeltaX_CSC3_Pl3 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC3_Pl3 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 4:
	    hDeltaX_CSC4_Pl3 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC4_Pl3 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 5:
	    hDeltaX_CSC5_Pl3 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC5_Pl3 -> Fill(y[counter]-y_prev[counter]);;
	      //  cout <<  "3: " << Eta(x[counter],y[counter],z[counter]) << " " << Phi(x[counter],y[counter]) << endl;
	    break;

	  }
	  if(z[counter]<0)
	  switch(sex){
	  case 0:
	    hDeltaMENOX_CSC0_Pl3 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC0_Pl3 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 1:
	    hDeltaMENOX_CSC1_Pl3 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC1_Pl3 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 2:
	    hDeltaMENOX_CSC2_Pl3 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC2_Pl3 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 3:
	    hDeltaMENOX_CSC3_Pl3 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC3_Pl3 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 4:
	    hDeltaMENOX_CSC4_Pl3 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC4_Pl3 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 5:
	    hDeltaMENOX_CSC5_Pl3 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC5_Pl3 -> Fill(y[counter]-y_prev[counter]);;
	      //  cout <<  "3: " << Eta(x[counter],y[counter],z[counter]) << " " << Phi(x[counter],y[counter]) << endl;
	    break;

	  }



	  TVector3 V3(x[counter],y[counter],z[counter]);
	  TVector2 V2(x_prev[counter],y_prev[counter]);
	  TVector2 eV2(eeX[counter],eeY[counter]);


	  _Pl3.push_back(V3);
	  _Pl3_prev.push_back(V2);
	  _all.push_back(V3);
	  _all_prev.push_back(V2);
	  _all_error.push_back(eV2);


 	}
	if(fabs(z[counter])>9500){
	  if(z[counter]>0)
	  switch(sex){
	  case 0:
	    hDeltaX_CSC0_Pl4 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC0_Pl4 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 1:
	    hDeltaX_CSC1_Pl4 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC1_Pl4 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 2:
	    hDeltaX_CSC2_Pl4 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC2_Pl4 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 3:
	    hDeltaX_CSC3_Pl4 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC3_Pl4 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 4:
	    hDeltaX_CSC4_Pl4 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC4_Pl4 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 5:
	    hDeltaX_CSC5_Pl4 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaY_CSC5_Pl4 -> Fill(y[counter]-y_prev[counter]);;
	      //  cout << "4: " <<  Eta(x[counter],y[counter],z[counter]) << " " << Phi(x[counter],y[counter]) << endl;
	    break;

	  }
	  if(z[counter]<0)
	  switch(sex){
	  case 0:
	    hDeltaMENOX_CSC0_Pl4 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC0_Pl4 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 1:
	    hDeltaMENOX_CSC1_Pl4 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC1_Pl4 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 2:
	    hDeltaMENOX_CSC2_Pl4 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC2_Pl4 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 3:
	    hDeltaMENOX_CSC3_Pl4 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC3_Pl4 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 4:
	    hDeltaMENOX_CSC4_Pl4 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC4_Pl4 -> Fill(y[counter]-y_prev[counter]);;
	    break;
	  case 5:
	    hDeltaMENOX_CSC5_Pl4 -> Fill(x[counter]-x_prev[counter]);;
	    hDeltaMENOY_CSC5_Pl4 -> Fill(y[counter]-y_prev[counter]);;
	      //  cout << "4: " <<  Eta(x[counter],y[counter],z[counter]) << " " << Phi(x[counter],y[counter]) << endl;
	    break;

	  }



	  TVector3 V3(x[counter],y[counter],z[counter]);
	  TVector2 V2(x_prev[counter],y_prev[counter]);
	  TVector2 eV2(eeX[counter],eeY[counter]);

	  _Pl4.push_back(V3);
	  _Pl4_prev.push_back(V2);
	  _all.push_back(V3);
	  _all_prev.push_back(V2);
	  _all_error.push_back(eV2);
 	}

      }
      counter++;
    }
  }

  delete [] x;
  delete [] y;
  delete [] z;
  delete [] x_prev;
  delete [] y_prev;


  TFitterMinuit * minuit = new TFitterMinuit();



  MyFCN22 fcn2(_all,_all_prev,_all_error);



  minuit->SetMinuitFCN(&fcn2);
  double startX = 0;
  double startY = 0;
// SetParameter: name, param value, param error, min, max
// if not limited (vhigh <= vlow)
  minuit->SetParameter(0,"x0",startX,0.01,-5,5);
  minuit->SetParameter(1,"y0",startY,0.01,-5,5);
  minuit->SetParameter(2,"x1",startX,0.01,-5,5);
  minuit->SetParameter(3,"y1",startY,0.01,-5,5);
  minuit->SetParameter(4,"x2",startX,0.01,-5,5);
  minuit->SetParameter(5,"y2",startY,0.01,-5,5);
  minuit->SetParameter(6,"x3",startX,0.01,-5,5);
  minuit->SetParameter(7,"y3",startY,0.01,-5,5);
  minuit->SetParameter(8,"x4",startX,0.01,-5,5);
  minuit->SetParameter(9,"y4",startY,0.01,-5,5);
  minuit->SetParameter(10,"ang0",0,0.001,-2*_3DEG,2*_3DEG);
  minuit->SetParameter(11,"ang1",0,0.001,-2*_3DEG,2*_3DEG);
  minuit->SetParameter(12,"ang2",0,0.001,-2*_3DEG,2*_3DEG);
  minuit->SetParameter(13,"ang3",0,0.001,-2*_3DEG,2*_3DEG);
  minuit->SetParameter(14,"ang4",0,0.001,-2*_3DEG,2*_3DEG);

//		  minuit->SetParameter(2,"theta",0,0.1,0,0);
  minuit->SetPrintLevel(5);
// create Minimizer (default is Migrad)
  minuit->CreateMinimizer();
  int iret = minuit->Minimize();



  for(unsigned int j=0; j<5; j++)
    {
      unsigned int aj = j*2;
      if(iret == 0){
	allineo_x[braccio][j][sestante]=minuit->GetParameter(aj);
	allineo_y[braccio][j][sestante]=minuit->GetParameter(aj+1);
	allineo_theta[braccio][j][sestante]=minuit->GetParameter(10+j);
      }
      else{
	allineo_x[braccio][j][sestante]=0;
	allineo_y[braccio][j][sestante]=0;
	allineo_theta[braccio][j][sestante]=0;
      }
    }

  double chi2, edm, errdef;
  int nvpar, nparx;
  if(iret == 0){
    minuit->GetStats(chi2,edm,errdef,nvpar,nparx);
    std::cout << "                         Chi2 Fit = " << chi2 <<std::endl;
  }


  return iret;


}



float T1InternalAlignment2::Eta(float x,float y,float z){
  float xyt;
  float c=0;
  float eta2;
  xyt = sqrt(x*x + y*y);
  //theta
  if(z>0) c = atan(xyt/z);
  if(z<0) c = atan(xyt/z)+3.14159;
  if(z==0) {c = 3.14159;}
  //pseudorapidity
  eta2 = -log(tan(c/2.));
  return eta2;
}
float T1InternalAlignment2::Phi(float x,float y){
  float c=0;
  if(x>0 && y>0) c = atan(y/x);
  if(x<0) c = atan(y/x)+3.14159;
  if(x>0 && y<0) c = atan(y/x)+6.28318;
  return c;
}

int T1InternalAlignment2::sextant(float x, float y, float z){
    int sestante = -1;
    float fi = Phi(x,y);
//  cout << " FI " << fi << endl;
    if(z>0){
      if(plane(z)==4){
	if( fi < _30DEG || fi > (2*PI - _30DEG) ) sestante = 5;
	if(fi > (_60DEG -_30DEG) && fi < (_60DEG +_30DEG) ) sestante = 0;
	if(fi > (_120DEG -_30DEG) && fi < (_120DEG +_30DEG) ) sestante = 1;
	if(fi > (PI -_30DEG) && fi < (PI +_30DEG) ) sestante = 2;
	if(fi > (_240DEG -_30DEG) && fi < (_240DEG +_30DEG) ) sestante = 3;
	if(fi > (_300DEG -_30DEG) && fi < (_300DEG +_30DEG) ) sestante = 4;
      }
      if(plane(z)==3){
	if( fi < (_30DEG-_3DEG) || fi > (2*PI - _30DEG-_3DEG) ) sestante = 5;
	if(fi > (_60DEG -_30DEG-_3DEG) && fi < (_60DEG +_30DEG-_3DEG) ) sestante = 0;
	if(fi > (_120DEG -_30DEG-_3DEG) && fi < (_120DEG +_30DEG-_3DEG) ) sestante = 1;
	if(fi > (PI -_30DEG-_3DEG) && fi < (PI +_30DEG-_3DEG) ) sestante = 2;
	if(fi > (_240DEG -_30DEG-_3DEG) && fi < (_240DEG +_30DEG-_3DEG) ) sestante = 3;
	if(fi > (_300DEG -_30DEG-_3DEG) && fi < (_300DEG +_30DEG-_3DEG) ) sestante = 4;
      }
      if(plane(z)==2){
	if( fi < (_30DEG-_3DEG-_3DEG) || fi > (2*PI - _30DEG-_3DEG-_3DEG) ) sestante = 5;
	if(fi > (_60DEG -_30DEG-_3DEG-_3DEG) && fi < (_60DEG +_30DEG-_3DEG-_3DEG) ) sestante = 0;
	if(fi > (_120DEG -_30DEG-_3DEG-_3DEG) && fi < (_120DEG +_30DEG-_3DEG-_3DEG) ) sestante = 1;
	if(fi > (PI -_30DEG-_3DEG-_3DEG) && fi < (PI +_30DEG-_3DEG-_3DEG) ) sestante = 2;
	if(fi > (_240DEG -_30DEG-_3DEG-_3DEG) && fi < (_240DEG +_30DEG-_3DEG-_3DEG) ) sestante = 3;
	if(fi > (_300DEG -_30DEG-_3DEG-_3DEG) && fi < (_300DEG +_30DEG-_3DEG-_3DEG) ) sestante = 4;
      }
      if(plane(z)==1){
	if( fi < (_30DEG+_3DEG) || fi > (2*PI - _30DEG+_3DEG) ) sestante = 5;
	if(fi > (_60DEG -_30DEG+_3DEG) && fi < (_60DEG +_30DEG+_3DEG) ) sestante = 0;
	if(fi > (_120DEG -_30DEG+_3DEG) && fi < (_120DEG +_30DEG+_3DEG) ) sestante = 1;
	if(fi > (PI -_30DEG+_3DEG) && fi < (PI +_30DEG+_3DEG) ) sestante = 2;
	if(fi > (_240DEG -_30DEG+_3DEG) && fi < (_240DEG +_30DEG+_3DEG) ) sestante = 3;
	if(fi > (_300DEG -_30DEG+_3DEG) && fi < (_300DEG +_30DEG+_3DEG) ) sestante = 4;
      }
      if(plane(z)==0){
	if( fi < (_30DEG+_3DEG+_3DEG) || fi > (2*PI - _30DEG+_3DEG+_3DEG) ) sestante = 5;
	if(fi > (_60DEG -_30DEG+_3DEG+_3DEG) && fi < (_60DEG +_30DEG+_3DEG+_3DEG) ) sestante = 0;
	if(fi > (_120DEG -_30DEG+_3DEG+_3DEG) && fi < (_120DEG +_30DEG+_3DEG+_3DEG) ) sestante = 1;
	if(fi > (PI -_30DEG+_3DEG+_3DEG) && fi < (PI +_30DEG+_3DEG+_3DEG) ) sestante = 2;
	if(fi > (_240DEG -_30DEG+_3DEG+_3DEG) && fi < (_240DEG +_30DEG+_3DEG+_3DEG) ) sestante = 3;
	if(fi > (_300DEG -_30DEG+_3DEG+_3DEG) && fi < (_300DEG +_30DEG+_3DEG+_3DEG) ) sestante = 4;
      }

    }


    if(z<0){
      if(plane(z)==4){
	if( fi < _30DEG || fi > (2*PI - _30DEG) ) sestante = 2;
	if(fi > (_60DEG -_30DEG) && fi < (_60DEG +_30DEG) ) sestante = 1;
	if(fi > (_120DEG -_30DEG) && fi < (_120DEG +_30DEG) ) sestante = 0;
	if(fi > (PI -_30DEG) && fi < (PI +_30DEG) ) sestante = 5;
	if(fi > (_240DEG -_30DEG) && fi < (_240DEG +_30DEG) ) sestante = 4;
	if(fi > (_300DEG -_30DEG) && fi < (_300DEG +_30DEG) ) sestante = 3;
      }
      if(plane(z)==3){
	if( fi < (_30DEG+_3DEG) || fi > (2*PI - _30DEG+_3DEG) ) sestante = 2;
	if(fi > (_60DEG -_30DEG+_3DEG) && fi < (_60DEG +_30DEG+_3DEG) ) sestante = 1;
	if(fi > (_120DEG -_30DEG+_3DEG) && fi < (_120DEG +_30DEG+_3DEG) ) sestante = 0;
	if(fi > (PI -_30DEG+_3DEG) && fi < (PI +_30DEG+_3DEG) ) sestante = 5;
	if(fi > (_240DEG -_30DEG+_3DEG) && fi < (_240DEG +_30DEG+_3DEG) ) sestante = 4;
	if(fi > (_300DEG -_30DEG+_3DEG) && fi < (_300DEG +_30DEG+_3DEG) ) sestante = 3;
      }
      if(plane(z)==2){
	if( fi < (_30DEG+_3DEG+_3DEG) || fi > (2*PI - _30DEG+_3DEG+_3DEG) ) sestante = 2;
	if(fi > (_60DEG -_30DEG+_3DEG+_3DEG) && fi < (_60DEG +_30DEG+_3DEG+_3DEG) ) sestante = 1;
	if(fi > (_120DEG -_30DEG+_3DEG+_3DEG) && fi < (_120DEG +_30DEG+_3DEG+_3DEG) ) sestante = 0;
	if(fi > (PI -_30DEG+_3DEG+_3DEG) && fi < (PI +_30DEG+_3DEG+_3DEG) ) sestante = 5;
	if(fi > (_240DEG -_30DEG+_3DEG+_3DEG) && fi < (_240DEG +_30DEG+_3DEG+_3DEG) ) sestante = 4;
	if(fi > (_300DEG -_30DEG+_3DEG+_3DEG) && fi < (_300DEG +_30DEG+_3DEG+_3DEG) ) sestante = 3;
      }
      if(plane(z)==1){
	if( fi < (_30DEG-_3DEG) || fi > (2*PI - _30DEG-_3DEG) ) sestante = 2;
	if(fi > (_60DEG -_30DEG-_3DEG) && fi < (_60DEG +_30DEG-_3DEG) ) sestante = 1;
	if(fi > (_120DEG -_30DEG-_3DEG) && fi < (_120DEG +_30DEG-_3DEG) ) sestante = 0;
	if(fi > (PI -_30DEG-_3DEG) && fi < (PI +_30DEG-_3DEG) ) sestante = 5;
	if(fi > (_240DEG -_30DEG-_3DEG) && fi < (_240DEG +_30DEG-_3DEG) ) sestante = 4;
	if(fi > (_300DEG -_30DEG-_3DEG) && fi < (_300DEG +_30DEG-_3DEG) ) sestante = 3;
      }
      if(plane(z)==0){
	if( fi < (_30DEG-_3DEG-_3DEG) || fi > (2*PI - _30DEG-_3DEG-_3DEG) ) sestante = 2;
	if(fi > (_60DEG -_30DEG-_3DEG-_3DEG) && fi < (_60DEG +_30DEG-_3DEG-_3DEG) ) sestante = 1;
	if(fi > (_120DEG -_30DEG-_3DEG-_3DEG) && fi < (_120DEG +_30DEG-_3DEG-_3DEG) ) sestante = 0;
	if(fi > (PI -_30DEG-_3DEG-_3DEG) && fi < (PI +_30DEG-_3DEG-_3DEG) ) sestante = 5;
	if(fi > (_240DEG -_30DEG-_3DEG-_3DEG) && fi < (_240DEG +_30DEG-_3DEG-_3DEG) ) sestante = 4;
	if(fi > (_300DEG -_30DEG-_3DEG-_3DEG) && fi < (_300DEG +_30DEG-_3DEG-_3DEG) ) sestante = 3;
      }

    }


    return sestante;
  }
  int T1InternalAlignment2::plane(float z) {

    int piano = -1;
    if(fabs(z)<8000){
      piano = 0;

    }
    if(fabs(z)<8400 && fabs(z)>8000){
      piano = 1;
    }
    if(fabs(z)<9000 && fabs(z)>8400){
      piano = 2;
    }
    if(fabs(z)<9500 && fabs(z)>9000){
      piano = 3;
    }
    if(fabs(z)>9500){
      piano = 4;
    }
    return piano;

  }

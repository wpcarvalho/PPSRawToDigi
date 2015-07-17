/*
  Created by Fabrizio Ferro - INFN Genova for TOTEM
*/
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2D.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2DCollection.h"
#include "RecoTotemT1T2/T1RoadProducer/interface/T1RoadAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"


#include <vector>
// for pair class
#include <utility>

#include "TFile.h"
#include "TH1.h"

#include <CLHEP/Vector/LorentzVector.h>

#ifndef PI
#define PI 3.141592653589793
#endif 

using std::vector;
using std::pair;

//#define _DEBUG_

T1RoadAnalyzer::T1RoadAnalyzer(const edm::ParameterSet& iConfig)
{
  t1RecHit2DCollectionMyRecoCollLabel = iConfig.getParameter<edm::InputTag>("T1RecHit2DCollectionMyRecoCollLabel");
  t1RoadCollectionLabel = iConfig.getParameter<edm::InputTag>("T1RoadCollectionLabel");
  simVertexContainerLabel = iConfig.getParameter<edm::InputTag>("SimVertexContainerLabel");
  simTrackContainerLabel = iConfig.getParameter<edm::InputTag>("SimTrackContainerLabel");
  std::string thegeometryfile = "Geometry/TotemGeometry/data/T1_data_geometry.dat";
  layer = new T1Geometry(thegeometryfile);
}


T1RoadAnalyzer::~T1RoadAnalyzer()
{
 
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
T1RoadAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;




  std::cout <<"Evento "<<iEvent.id().event() <<std::endl;


  edm::Handle<T1RecHit2DCollection>myRecoColl;
  iEvent.getByLabel(t1RecHit2DCollectionMyRecoCollLabel, myRecoColl);
  vector<pair<T1RecHitGlobal,int> > * global_copy_of_myRecoColl = new vector<pair<T1RecHitGlobal,int> >;

  T1RecHit2DCollection::const_iterator T1Reco_it;
  vector<pair<T1RecHitGlobal,int> >::iterator T1RecoGlobal_it;
  vector<pair<T1RecHitGlobal,int> >::iterator T1RecoGlobal_it2;
  vector<T1RecHitGlobal>::iterator T1RecoGlobal_it3;
  vector<pair<T1RecHitGlobal,int> >::iterator T1RecoGlobal_it4;

  
  T1ChamberSpecs parametri;
 

  int numHits=0;
  nHits =0;

  for(T1Reco_it=myRecoColl->begin();T1Reco_it!=myRecoColl->end();T1Reco_it++){
    ++numHits; 
    //store rec hits in temporary collection
#ifdef _DEBUG_
    std::cout << "                                                             " << numHits<<std::endl;
#endif


    float xxxR = layer->xFromLocal2BeamSystem((*T1Reco_it).t1DetId().rawId(),(*T1Reco_it).localPosition().x() );
    float yyyR = layer->yFromLocal2BeamSystem((*T1Reco_it).t1DetId().rawId(),(*T1Reco_it).localPosition().y() );

    float xxxxR=0;
    float yyyyR=0;
    layer->RotationLocal2Global((*T1Reco_it).t1DetId().rawId(),xxxR,yyyR,xxxxR,yyyyR);

    float errorZ = 10.; // 10mm = 1cm;

    GlobalError theError ( (*T1Reco_it).localPositionError().xx(),0.,(*T1Reco_it).localPositionError().yy(),0.,0.,errorZ);
    GlobalPoint thePos( xxxxR, yyyyR, layer->Zeta( (*T1Reco_it).t1DetId().rawId() ) );
	  
    //	  T1RecHitGlobal *myGlobalHit = new T1RecHitGlobal(thePos, theError);
    T1RecHitGlobal myGlobalHit(thePos, theError);
    //	  pair<T1RecHitGlobal, int> *myGlobalHitPair=new pair<T1RecHitGlobal, int>(myGlobalHit,0);
    //	  global_copy_of_myRecoColl->push_back(*myGlobalHitPair);
   
    global_copy_of_myRecoColl->push_back(pair<T1RecHitGlobal, int>(myGlobalHit,0) );

  }


  for(T1RecoGlobal_it=global_copy_of_myRecoColl->begin();T1RecoGlobal_it!=global_copy_of_myRecoColl->end();T1RecoGlobal_it++){
#ifdef _DEBUG_
    std::cout << (*T1RecoGlobal_it).first <<std::endl;
#endif

    nHits++;
    rx[nHits]=(*T1RecoGlobal_it).first.GlobalPosition().x();
    ry[nHits]=(*T1RecoGlobal_it).first.GlobalPosition().y();
    rz[nHits]=(*T1RecoGlobal_it).first.GlobalPosition().z();
    ee[nHits]=Eta(rx[nHits],ry[nHits],rz[nHits]);
    ff[nHits]=Phi(rx[nHits],ry[nHits]);
    evto[nHits]=iEvent.id().event();
  }



  int RoadNum=0;
  int PriRoadNum=0;
  int PriTrackNum=0;

  edm::Handle<T1RoadCollection> roadCollection;
  iEvent.getByLabel(t1RoadCollectionLabel, roadCollection);
  //  iEvent.getByLabel("t1roads",roadCollection);

  std::cout << " Taglia della road collection: " << roadCollection->size() <<std::endl;


  vector< pair<double,double> > PrimarySimTracks;
  vector< pair<double,double> >::iterator PrimarySimTracks_it;


  // get G4 Vertexes
  std::vector<SimVertex> theSimVertexes;

  Handle<SimVertexContainer> SimVtx;
  iEvent.getByLabel(simVertexContainerLabel, SimVtx);
  theSimVertexes.insert(theSimVertexes.end(),SimVtx->begin(),SimVtx->end());


  int hjk=0;
  for(std::vector<SimVertex>::iterator isimvtx = theSimVertexes.begin(); isimvtx != theSimVertexes.end(); ++isimvtx){
#ifdef _DEBUG_
    cout << " Z Vertice["<<hjk<<"] = "<<isimvtx->position().z() << "  Parent indec = "<<isimvtx->parentIndex()<<endl;
#endif
    hjk++;
  }


  //get G4 tracks

  edm::Handle<edm::SimTrackContainer> G4TrkContainer;
  iEvent.getByLabel(simTrackContainerLabel, G4TrkContainer);
  if (!G4TrkContainer.isValid()) {
    edm::LogError("TrackerHitAnalyzer::analyze")
      << "Unable to find SimTrack in event!";
    return;
  }
 
  edm::SimTrackContainer::const_iterator itTrk;
  for (itTrk = G4TrkContainer->begin(); itTrk != G4TrkContainer->end(); 
       ++itTrk) {

    //    cout << "itTrk = "<< itTrk << endl;
    double eta =0, phi =0, p =0;
    //    const HepLorentzVector& G4Trk = itTrk->momentum();
    const CLHEP::HepLorentzVector G4Trk(itTrk->momentum().x(), itTrk->momentum().y(),itTrk->momentum().z(),itTrk->momentum().e() );
    p =sqrt(G4Trk[0]*G4Trk[0]+G4Trk[1]*G4Trk[1]+G4Trk[2]*G4Trk[2]);
    if ( p == 0) 
      edm::LogError("TrackerHitAnalyzer::analyze") 
	<< "TrackerTest::INFO: Primary has p = 0 ";
    else {
      double costheta  = G4Trk[2]/p;
      double theta = acos(TMath::Min(TMath::Max(costheta, -1.),1.));
      eta = -log(tan(theta/2));          
      if ( G4Trk[0] != 0 || G4Trk[1] != 0) phi = atan2(G4Trk[1],G4Trk[0]);
	
      if(phi<0)phi = 2*PI + phi;


      float carica =  itTrk->charge();

      if(fabs(eta)>3.1 && fabs(eta)<4.7 && carica != 0){
#ifdef _DEBUG_
      int IndiceVertice = itTrk->vertIndex();
	std::cout << "TRACCIA G4: Eta = "<< eta<< "  Phi = " << phi << " Vert num. = "<<IndiceVertice <<  "   V = (" << theSimVertexes[IndiceVertice].position().x()<<","<<  theSimVertexes[IndiceVertice].position().y()<<","<<  theSimVertexes[IndiceVertice].position().z()<<")  PDG = " <<itTrk->type() << "  Charge = "<< itTrk->charge() << "   P = " << p
		  << std::endl;
#endif
	pair<double,double> *myprimarytrack = new pair<double,double>(eta,phi);
	PrimarySimTracks.push_back(*myprimarytrack);
	PriTrackNum++;
      }
    }
  }


  /////////////////////////////////////////////////////////////



  T1RoadCollection::const_iterator RC_it;

  for(RC_it=roadCollection->begin(); RC_it!=roadCollection->end(); RC_it++){

#ifdef _DEBUG_
    std::cout << "Road Size: " << (*RC_it).size() << std::endl;
#endif

    hRoadSize->Fill( (*RC_it).size() );

    bool FlagPriRoadNum=false;

    for(PrimarySimTracks_it=PrimarySimTracks.begin(); PrimarySimTracks_it!=PrimarySimTracks.end(); PrimarySimTracks_it++){
      float de = fabs( (*PrimarySimTracks_it).first -    Eta((*RC_it)[0].GlobalPosition().x(),(*RC_it)[0].GlobalPosition().y(),(*RC_it)[0].GlobalPosition().z()) );
      float df = fabs( (*PrimarySimTracks_it).second - Phi((*RC_it)[0].GlobalPosition().x(),(*RC_it)[0].GlobalPosition().y()) );
      if(df>3.14159)df=2*3.14159-df;

#ifdef _DEBUG_
      
      std::cout << " Deta e Dfi " << de << " " << df << std::endl; 
#endif

      if(de < 0.2 && df < 0.2)
	FlagPriRoadNum=true;
    }
      
    RoadNum++;
    if(FlagPriRoadNum)PriRoadNum++;

    float de =0;
    float df =0;
    float de_temp =0;
    float df_temp =0;

    int piano[5]={0,0,0,0,0};
    int somma =0;
      

    for(unsigned int iy = 0; iy < (*RC_it).size(); iy++){
      de_temp = fabs( Eta((*RC_it)[0].GlobalPosition().x(),(*RC_it)[0].GlobalPosition().y(),(*RC_it)[0].GlobalPosition().z())-Eta((*RC_it)[iy].GlobalPosition().x(),(*RC_it)[iy].GlobalPosition().y(),(*RC_it)[iy].GlobalPosition().z()) );
      df = fabs( Phi((*RC_it)[0].GlobalPosition().x(),(*RC_it)[0].GlobalPosition().y())-Phi((*RC_it)[iy].GlobalPosition().x(),(*RC_it)[iy].GlobalPosition().y()) );
      if(df>3.14159)df=2*3.14159-df;
      
      if(de_temp>de)de = de_temp;
      if(df_temp>df)df = df_temp;
      
      double za = fabs( (*RC_it)[iy].GlobalPosition().z() );

      if( za < 7600 && za >7400)piano[0]=1;
      if( za < 8300 && za >8000)piano[1]=1;
      if( za < 9000 && za >8500)piano[2]=1;
      if( za < 9600 && za >9000)piano[3]=1;
      if( za < 10300 && za >9600)piano[4]=1;

      

#ifdef _DEBUG_
      std::cout << (*RC_it)[iy] <<std::endl;
#endif
    }
    //maximum de and df inside a track
    for(int iop=0; iop<5;iop++)somma= somma + piano[iop];
    hRoadPlaneSize->Fill(somma);
    hDeltaEta->Fill(de);
    hDeltaPhi->Fill(df);

    if(FlagPriRoadNum){
      hDeltaEtaPri->Fill(de);
      hDeltaPhiPri->Fill(df);
    }
      

  }
    
  

  hRoadNumber->Fill(RoadNum);
  hPrimaryRoadNumber->Fill(PriRoadNum);
  hPrimaryTracks->Fill(PriTrackNum);
  if(PriTrackNum>0)
    hPrimaryRoadTrackRatio->Fill((float)PriRoadNum/(float)PriTrackNum);
  if(PriTrackNum>0)
    hRoadTrackRatio->Fill((float)RoadNum/(float)PriTrackNum);
  if(RoadNum>0)
    hPrimaryRoadRoadRatio->Fill((float)PriRoadNum/(float)RoadNum);
  hLostPrimaryTracks->Fill(PriTrackNum-PriRoadNum);
  if((PriTrackNum-PriRoadNum)>=0 && PriTrackNum>0)
    hLostPrimaryTrackRatio->Fill((float)(PriTrackNum-PriRoadNum)/(float)PriTrackNum);
  if((PriTrackNum-PriRoadNum)<0 && PriTrackNum>0)
    hLostPrimaryTrackRatio->Fill(1);

  std::cout <<"Evt "<<iEvent.id().event() <<"  Primary tracks in range = " << PriTrackNum << " ##  Primary Roads = " << PriRoadNum << " ## Roads = "<< RoadNum <<std::endl;  
  std::cout << std::endl;



  tree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
T1RoadAnalyzer::beginJob()
{

  theFile = new TFile("roadsFile.root","RECREATE");
  hDeltaEta = new TH1D("deltaeta","deltaeta",1000,0,1);
  hDeltaPhi = new TH1D("deltaphi","deltaphi",1000,0,1);
  hDeltaEtaPri = new TH1D("deltaetaPri","deltaetaPri",1000,0,1);
  hDeltaPhiPri = new TH1D("deltaphiPri","deltaphiPri",1000,0,1);
  hRoadNumber = new TH1D("RoadNum","RoadNum",100,-0.5,99.5);
  hRoadSize = new TH1D("RoadSize","RoadSize",40,-0.5,39.5);
  hRoadPlaneSize = new TH1D("RoadPlaneSize","RoadPlaneSize",10,-0.5,9.5);

  hPrimaryRoadNumber = new TH1D("PrimaryRoadNum","PrimaryRoadNum",100,-0.5,99.5);
 
  hPrimaryTracks = new TH1D("PrimaryTracks","PrimaryTracks",100,-0.5,99.5);
  
  hPrimaryRoadTrackRatio = new TH1D("PrimaryRoadTrackRatio","PrimaryRoadTrackRatio",11,-0.05,1.05);
  hRoadTrackRatio = new TH1D("RoadTrackRatio","RoadTrackRatio",11,-0.05,1.05);
  hPrimaryRoadRoadRatio = new TH1D("PrimaryRoadRoadRatio","PrimaryRoadRoadRatio",11,-0.05,1.05);
  hLostPrimaryTracks = new TH1D("LostPrimaryTracks","LostPrimaryTracks",100,-0.5,99.5);
  hLostPrimaryTrackRatio = new TH1D("LostPrimaryTrackRatio","LostPrimaryTrackRatio",11,-0.05,1.05);


  tree = new TTree("tree", "CSC hits");
  tree->Branch("nHits", &nHits, "nHits/i");
  tree->Branch("rx", rx, "rx[nHits]/D");
  tree->Branch("ry", ry, "ry[nHits]/D");
  tree->Branch("rz", rz, "rz[nHits]/D");
  tree->Branch("ee", ee, "ee[nHits]/D");
  tree->Branch("ff", ff, "ff[nHits]/D");
  tree->Branch("evto", evto, "evto[nHits]/I");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
T1RoadAnalyzer::endJob() {
  theFile->Write();
  theFile->Close();
}

float T1RoadAnalyzer::Eta(float x,float y,float z){
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
float T1RoadAnalyzer::Phi(float x,float y){
  float c=0;
  if(x>0 && y>0) c = atan(y/x);
  if(x<0) c = atan(y/x)+3.14159;
  if(x>0 && y<0) c = atan(y/x)+6.28318;
  return c;
}

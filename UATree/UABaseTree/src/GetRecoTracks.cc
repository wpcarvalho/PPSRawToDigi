// Genaral Tracks and Vertex
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/src/Vertex.cc"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// UABaseTree Analysis class decleration
#include "UATree/UABaseTree/interface/UABaseTree.h"

bool TrackDebug = false;

void UABaseTree::GetRecoTracks(const edm::Event& iEvent, const string TrackCollName , vector<MyTracks>& TrackVector )
{


  TrackVector.clear(); 

  MyTracks mytrack;

  Handle<TrackCollection> tracks;
  iEvent.getByLabel(TrackCollName,tracks);

  for(TrackCollection::const_iterator tr = tracks->begin(); tr!=tracks->end(); ++tr){

    this->FillTrack(*tr , mytrack);
    TrackVector.push_back(mytrack);
    
    if(TrackDebug) mytrack.Print();
    
  }

}


void UABaseTree::GetAllTracks( const edm::Event& iEvent ){
  for(vector<InputTag>::iterator icoll = tracks_.begin() ; icoll!= tracks_.end() ; ++icoll)
    this->GetRecoTracks(iEvent , icoll->label() , allTracks[icoll->label()] );
}



void UABaseTree::FillTrack(const Track& intrack , MyTracks& outtrack , Double_t mass){

  outtrack.Reset();

  outtrack.charge = intrack.charge();
  outtrack.SetPtEtaPhiM(intrack.pt(),intrack.eta(),intrack.phi(), mass );  

  outtrack.nhit	              =  intrack.found() + intrack.lost();
  outtrack.nValidPixelHits    =  intrack.hitPattern().numberOfValidPixelHits  ();
  outtrack.nValidStripHits    =  intrack.hitPattern().numberOfValidStripHits  ();
  outtrack.nValidMuonCSCHits  =  intrack.hitPattern().numberOfValidMuonCSCHits();
  outtrack.nValidMuonDTHits   =  intrack.hitPattern().numberOfValidMuonDTHits ();
  outtrack.nValidMuonRPCHits  =  intrack.hitPattern().numberOfValidMuonRPCHits();
  
  outtrack.chi2n   =  intrack.normalizedChi2();
  outtrack.dz	   =  intrack.dz();
  outtrack.d0	   =  intrack.d0();
  outtrack.edz	   =  intrack.dzError();
  outtrack.ed0	   =  intrack.d0Error();
  outtrack.ept	   =  intrack.ptError();

  outtrack.vx    =  intrack.referencePoint().x();
  outtrack.vy    =  intrack.referencePoint().y();
  outtrack.vz    =  intrack.referencePoint().z();

  outtrack.quality[0] = intrack.quality(TrackBase::qualityByName("loose"));
  outtrack.quality[1] = intrack.quality(TrackBase::qualityByName("tight"));
  outtrack.quality[2] = intrack.quality(TrackBase::qualityByName("highPurity"));
  outtrack.quality[3] = intrack.quality(TrackBase::qualityByName("confirmed"));
  outtrack.quality[4] = intrack.quality(TrackBase::qualityByName("goodIterative"));

  outtrack.trackAlgo = intrack.algo();

  // BeamSpot (id=0) and Vertex (id>0) Links
  outtrack.vtxid.clear();
  outtrack.vtxdxy.clear();
  outtrack.vtxdz.clear();

  for ( int i = 0 ; i != vtxid ; i++ ){
    outtrack.vtxid.push_back( i ); 
    outtrack.vtxdxy.push_back( intrack.dxy( vtxid_xyz[i] ) );
    outtrack.vtxdz.push_back(  intrack.dz( vtxid_xyz[i] )  );
  }
}

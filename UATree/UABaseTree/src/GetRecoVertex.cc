// Genaral Tracks and Vertex
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// UABaseTree Analysis class decleration
#include "UATree/UABaseTree/interface/UABaseTree.h"

bool RecoVtxDebug = false;

void UABaseTree::GetRecoVertex(const edm::Event& iEvent, const string VertexCollName , vector<MyVertex>& VertexVector )
{
   VertexVector.clear();
   MyVertex myvertex;

   Handle<reco::VertexCollection> vtxcoll ;
   iEvent.getByLabel(VertexCollName,vtxcoll);

   if (RecoVtxDebug) cout<<"name of the requested vertex collection: "<<VertexCollName<<endl;

   for(VertexCollection::const_iterator p=vtxcoll->begin(); p!= vtxcoll->end() ; ++p){

     myvertex.Reset();

     myvertex.id	= vtxid++;
     myvertex.x 	= p->x()  ;
     myvertex.y 	= p->y()  ;
     myvertex.z 	= p->z()  ;


     myvertex.ex	= p->xError()  ;
     myvertex.ey	= p->yError()  ;
     myvertex.ez	= p->zError()  ;
     
     myvertex.validity  = p->isValid() ;
     myvertex.fake	= p->isFake()  ;
     
     myvertex.chi2        = p->chi2();
     myvertex.ndof        = p->ndof();

     myvertex.ntracks	= p->tracksSize() ;
     
     //SumPt from TrackRef
     for(reco::Vertex::trackRef_iterator iTrack = p->tracks_begin(); iTrack!=p->tracks_end(); ++iTrack)
       myvertex.SumPtTracks += (*iTrack)->pt();
       
     
     VertexVector.push_back(myvertex);
     vtxid_xyz.push_back(p->position());
	   
     if (RecoVtxDebug) myvertex.Print();
   }   

}



void UABaseTree::GetAllVertices( const edm::Event& iEvent ){
  for(vector<InputTag>::iterator icoll = vertices_.begin() ; icoll!= vertices_.end() ; ++icoll)
    GetRecoVertex(iEvent , icoll->label() , allVertices[icoll->label()] );
}



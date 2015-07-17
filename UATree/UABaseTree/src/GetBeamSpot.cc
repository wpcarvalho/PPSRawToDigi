//-- Beam Spot
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "UATree/UABaseTree/interface/UABaseTree.h"

bool BeamSpotDebug = false;

void UABaseTree::GetBeamSpot(const edm::Event& iEvent , const string& inputTag , MyBeamSpot& beamSpot ) {

   Handle<reco::BeamSpot> BShandle;
   iEvent.getByLabel(inputTag,BShandle);

   const reco::BeamSpot* BS = BShandle.product();

   beamSpot.x           = BS->position().x();    
   beamSpot.y           = BS->position().y();    
   beamSpot.z           = BS->position().z();    

   beamSpot.ex          = BS->x0Error();
   beamSpot.ey          = BS->y0Error();
   beamSpot.ez          = BS->z0Error();

   beamSpot.sigmaZ      = BS->sigmaZ();
   beamSpot.dxdz        = BS->dxdz();
   beamSpot.dydz        = BS->dydz();

   beamSpot.esigmaZ     = BS->sigmaZ0Error();
   beamSpot.edxdz       = BS->dxdzError();
   beamSpot.edydz       = BS->dydzError();

   beamSpot.BeamWidthX  = BS->BeamWidthX() ;
   beamSpot.BeamWidthY  = BS->BeamWidthY() ;

   beamSpot.eBeamWidthX = BS->BeamWidthXError();
   beamSpot.eBeamWidthY = BS->BeamWidthYError();
   
   //Usefull for tracks to associate d0 & dz to beamspot/vertices
   vtxid++;
   vtxid_xyz.push_back(BS->position());
   

   if(BeamSpotDebug) beamSpot.Print();
}

void UABaseTree::GetAllBeamSpots(const edm::Event& iEvent){
  for(vector<InputTag>::iterator it = beamspots_.begin() ; it != beamspots_.end() ; ++it)
    this->GetBeamSpot(iEvent , it->label() , allBeamSpots[it->label()] );
}


/** 
 * Class T2Digitizer
 * 
 * Author: Erik Brücken / University of Helsinki
 *         Mirko Berretti / University of Siena & Pisa INFN
 * Email:  brucken@cc.helsinki.fi
 *         mirko.berretti@gmail.com  
 * Date:   2007-11-26
 */

#include "Utilities/Timing/interface/TimingReport.h" 
#include "SimTotem/T2Digitizer/interface/T2Digitizer.h"
#include "SimTotem/T2Digitizer/interface/T2DetectorHit.h"
#include "SimTotem/T2Digitizer/interface/T2StripHitSim.h"
#include "SimTotem/T2Digitizer/interface/T2PadHitSim.h"
#include "SimTotem/T2Digitizer/interface/T2StripElectronicsSim.h"
#include "SimTotem/T2Digitizer/interface/T2PadElectronicsSim.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include <iostream>
#include <vector>

T2Digitizer::T2Digitizer(const edm::ParameterSet & paraSet,const edm::EventSetup& iSetup,CLHEP::HepRandomEngine* rndEngine) 
  : theStripHitSim(new T2StripHitSim(paraSet.getParameter<edm::ParameterSet>("stripHit"))), 
    thePadHitSim(new T2PadHitSim(paraSet.getParameter<edm::ParameterSet>("padHit"))),
    theStripElectronicsSim(new T2StripElectronicsSim(paraSet.getParameter<edm::ParameterSet>("stripElec"),iSetup,rndEngine)),
    thePadElectronicsSim(new T2PadElectronicsSim(paraSet.getParameter<edm::ParameterSet>("padElec"),iSetup))   
{

  // Simulate_Misalignment= thePSimHitMisaligner->SimulationActivated();
  thePSimHitMisaligner=boost::shared_ptr<PSimHitMisaligner>(new PSimHitMisaligner(paraSet.getParameter<edm::ParameterSet>("Misalignment"),iSetup));
  theT2Geometry = new T2Geometry();
  createStripChargeMap();
  createPadChargeMap();

  TakeOnlyPrim=false;
  TakeOnlySec=false;
  TakeOnlyPrim=paraSet.getParameter<bool>("TakeOnlyPrim");
  TakeOnlySec=paraSet.getParameter<bool>("TakeOnlySec");
  
  
  AllowedDecayToCount.push_back(15); //10-12s 10-13s (Tau), will decay before T2.
  //AllowedDecayToCount.push_back(16); AllowedDecayToCount.push_back(17); AllowedDecayToCount.push_back(18); Tau and its neutrino exclude
  AllowedDecayToCount.push_back(23); //W 10-25s
  AllowedDecayToCount.push_back(24); //Z 
  AllowedDecayToCount.push_back(130); AllowedDecayToCount.push_back(310); AllowedDecayToCount.push_back(311); 
  //K0L, K0S,K0. 10-8 and 10-11 sec can be seen in T2. Have pi0 and charged decay mode: To review
   //K+ lifetime 10-8
  //AllowedDecayToCount.push_back(111); Pi0 Excluded   
  AllowedDecayToCount.push_back(211); //Pi+- lifetime is 10-8 but gamma~10 for E= 1GeV Can be seen in T2.
  AllowedDecayToCount.push_back(113); //pi+ + pi- decay
  //WARNING: to rewievw: pi+ + pi- and 2pi0 decay
  AllowedDecayToCount.push_back(213);  //rho+
  AllowedDecayToCount.push_back(215); //a2 meson ? 
  //Other quite rare I=1 mesons
  AllowedDecayToCount.push_back(117); AllowedDecayToCount.push_back(217); AllowedDecayToCount.push_back(119); AllowedDecayToCount.push_back(219); 
  //WARNING: Neutral decay mode is 72%. At the moment I drop it. Very short lifetime
  //AllowedDecayToCount.push_back(221); AllowedDecayToCount.push_back(331);//eta, eta prime   
  AllowedDecayToCount.push_back(223); AllowedDecayToCount.push_back(333); //w phi
  AllowedDecayToCount.push_back(225); AllowedDecayToCount.push_back(335); //f2 f2'
  AllowedDecayToCount.push_back(227); //w3 phi3
  AllowedDecayToCount.push_back(337); AllowedDecayToCount.push_back(229);  
  AllowedDecayToCount.push_back(411); AllowedDecayToCount.push_back(421); AllowedDecayToCount.push_back(413); 
  AllowedDecayToCount.push_back(423); AllowedDecayToCount.push_back(415); AllowedDecayToCount.push_back(425); AllowedDecayToCount.push_back(431); AllowedDecayToCount.push_back(433); 
  AllowedDecayToCount.push_back(435); AllowedDecayToCount.push_back(511); AllowedDecayToCount.push_back(521); AllowedDecayToCount.push_back(513); AllowedDecayToCount.push_back(523); 
  AllowedDecayToCount.push_back(515); AllowedDecayToCount.push_back(525); AllowedDecayToCount.push_back(531); AllowedDecayToCount.push_back(533); AllowedDecayToCount.push_back(535); 
  AllowedDecayToCount.push_back(541); AllowedDecayToCount.push_back(543); AllowedDecayToCount.push_back(545); AllowedDecayToCount.push_back(441); AllowedDecayToCount.push_back(443); 
  AllowedDecayToCount.push_back(445); AllowedDecayToCount.push_back(551); AllowedDecayToCount.push_back(553); AllowedDecayToCount.push_back(555); AllowedDecayToCount.push_back(557); 
  
  //AllowedDecayToCount.push_back(2212); AllowedDecayToCount.push_back(2112); //p and n excluded

  AllowedDecayToCount.push_back(2224); AllowedDecayToCount.push_back(2214); AllowedDecayToCount.push_back(2114); 
  AllowedDecayToCount.push_back(1114); AllowedDecayToCount.push_back(3122); AllowedDecayToCount.push_back(3222); AllowedDecayToCount.push_back(3212); 
  AllowedDecayToCount.push_back(3112); AllowedDecayToCount.push_back(3224); AllowedDecayToCount.push_back(3214); AllowedDecayToCount.push_back(3114); AllowedDecayToCount.push_back(3322); 
  AllowedDecayToCount.push_back(3312); AllowedDecayToCount.push_back(3324); AllowedDecayToCount.push_back(3314); AllowedDecayToCount.push_back(3334); AllowedDecayToCount.push_back(4122); 
  AllowedDecayToCount.push_back(4222); AllowedDecayToCount.push_back(4212); AllowedDecayToCount.push_back(4112); AllowedDecayToCount.push_back(4224); AllowedDecayToCount.push_back(4214); 
  AllowedDecayToCount.push_back(4114); AllowedDecayToCount.push_back(4232); AllowedDecayToCount.push_back(4132); AllowedDecayToCount.push_back(4322); AllowedDecayToCount.push_back(4312); 
  AllowedDecayToCount.push_back(4324); AllowedDecayToCount.push_back(4314); AllowedDecayToCount.push_back(4332); AllowedDecayToCount.push_back(4334); AllowedDecayToCount.push_back(4412); 
  AllowedDecayToCount.push_back(4422); AllowedDecayToCount.push_back(4414); AllowedDecayToCount.push_back(4424); AllowedDecayToCount.push_back(4432); AllowedDecayToCount.push_back(4434); 
  AllowedDecayToCount.push_back(4444); AllowedDecayToCount.push_back(5122); AllowedDecayToCount.push_back(5112); AllowedDecayToCount.push_back(5212); AllowedDecayToCount.push_back(5222); 
  AllowedDecayToCount.push_back(5114); AllowedDecayToCount.push_back(5214); AllowedDecayToCount.push_back(5224); AllowedDecayToCount.push_back(5132); AllowedDecayToCount.push_back(5232); 
  AllowedDecayToCount.push_back(5312); AllowedDecayToCount.push_back(5322); AllowedDecayToCount.push_back(5314); AllowedDecayToCount.push_back(5324); AllowedDecayToCount.push_back(5332); 
  AllowedDecayToCount.push_back(5334); AllowedDecayToCount.push_back(5142); AllowedDecayToCount.push_back(5242); AllowedDecayToCount.push_back(5412); AllowedDecayToCount.push_back(5422); 
  AllowedDecayToCount.push_back(5414); AllowedDecayToCount.push_back(5424); AllowedDecayToCount.push_back(5342); AllowedDecayToCount.push_back(5432); AllowedDecayToCount.push_back(5434); 
  AllowedDecayToCount.push_back(5442); AllowedDecayToCount.push_back(5444); AllowedDecayToCount.push_back(5512); AllowedDecayToCount.push_back(5522); AllowedDecayToCount.push_back(5514); 
  AllowedDecayToCount.push_back(5524); AllowedDecayToCount.push_back(5532); AllowedDecayToCount.push_back(5534); AllowedDecayToCount.push_back(5542); AllowedDecayToCount.push_back(5544); 
  AllowedDecayToCount.push_back(5554);





} // T2Digitizer


T2Digitizer::~T2Digitizer() {

  delete theStripHitSim;
  delete thePadHitSim;
  delete theStripElectronicsSim;
  delete thePadElectronicsSim;

  delete theT2Geometry;

  deleteStripChargeMap();
  deletePadChargeMap();

} // ~T2Digitizer


void T2Digitizer::process(MixCollection<PSimHit> & simHits, 
			  T2StripDigiCollection & stripDigis, 
			  T2PadDigiCollection & padDigis,
			  const HepMC::GenEvent* hepmcevt,
			  std::auto_ptr<edm::SimVertexContainer> theSimVertices,
			  std::auto_ptr<edm::SimTrackContainer> theSimTracks) {
  
  // map hits to T2 plane
  std::map<int, edm::PSimHitContainer> hitMap;
 

  for(MixCollection<PSimHit>::MixItr hitItr = simHits.begin(); hitItr != simHits.end(); ++hitItr) {
   

    if(TakeOnlyPrim || TakeOnlySec){

      //      for (edm::PSimHitContainer::const_iterator hitItr = simHits.begin(); hitItr != simHits.end();  ++hitItr) 

      int Hit_TrkId=hitItr->trackId(); 

      bool primaryhit=false;
      //Check now if this hit is from primary or secondary GEANT tracks

      SimTrack PythiaIPChParticleTrack;
      
      if(getTrack(Hit_TrkId,theSimTracks,PythiaIPChParticleTrack) == false){
	primaryhit=false;
	//std::cout<<" Trk Not found"<<std::endl;
      }
      else{	
	int motherbarcode=99; 
	 primaryhit=PrimaryTrackCondition(PythiaIPChParticleTrack,theSimVertices,theSimTracks,hepmcevt,motherbarcode);
	 //std::cout<<" Trk  found. Primary hit:"<<primaryhit<<"  TakeOnlyPrim is"<<TakeOnlyPrim<<std::endl;
      }
      
      if(primaryhit && TakeOnlyPrim){
	hitMap[hitItr->detUnitId()].push_back(*hitItr);
	//std::cout<<" Prim hit inserted"<<std::endl;
      }

      if((!primaryhit) && TakeOnlySec){
	hitMap[hitItr->detUnitId()].push_back(*hitItr);
	//std::cout<<" Sec hit inserted"<<std::endl;
      }

    }else{
       hitMap[hitItr->detUnitId()].push_back(*hitItr);
    }

  }

  // reset the chargecontainer to zero
  T2Digitizer::resetStripChargeMap();
  T2Digitizer::resetPadChargeMap();

  // now we run for each plane the digi-simulation
  for(std::map<int, edm::PSimHitContainer>::const_iterator hitMapItr = hitMap.begin();
      hitMapItr != hitMap.end(); ++hitMapItr) {
  
    theT2Geometry->setPlane((uint32_t)hitMapItr->first);
    const edm::PSimHitContainer & planeSimHits = hitMapItr->second;

    //std::cout<<"T2Digitizer: found " << planeSimHits.size() <<" hits in plane"<<std::endl;

    std::vector<T2DetectorHit> newStripHits, newPadHits;
    LogDebug("T2Digitizer") << "T2Digitizer: found " << planeSimHits.size() <<" hits in plane";
    
    // container of the percentage of charge taken by the strips 
    std::vector<double> chargeDivider;

    // Let's process the stripHits
    newStripHits.swap(theStripHitSim->simulate(theT2Geometry, planeSimHits, chargeDivider,thePSimHitMisaligner));

    // Let's process the padHits; no real simulation yet, just closest Pad is turned on
    newPadHits.swap(thePadHitSim->simulate(theT2Geometry, planeSimHits, chargeDivider,thePSimHitMisaligner));

    // here we put the charge of the strip signal into a container that 
    // exists as well for empty frames
    size_t nStripis = newStripHits.size();
      
    for( size_t i = 0; i < nStripis; ++i) {
		int sR = newStripHits[i].getRow();
		int sC = newStripHits[i].getCol();
		int adc = (int)newStripHits[i].getCharge();
		stripChargeMap[hitMapItr->first][sR + sC*256] += adc;
    }

    // here we put the charge of the pad signal into a container that 
    // exists as well for empty frames
    size_t nPadis = newPadHits.size();
      
    for( size_t i = 0; i < nPadis; ++i) {
		int pR = newPadHits[i].getRow();
		int pC = newPadHits[i].getCol();
		int adc = (int)newPadHits[i].getCharge();
		padChargeMap[hitMapItr->first][pR + pC*24] += adc;
    }
  }
  
  // Electronics simulation for strips including noise and threshold of vfat-chip
  theStripElectronicsSim->simulate(stripChargeMap);
  theStripElectronicsSim->fillDigis(stripDigis, stripChargeMap);

  // Electronics simulation for pads including noise and threshold of vfat-chip
  thePadElectronicsSim->simulate(padChargeMap);
  thePadElectronicsSim->fillDigis(padDigis, padChargeMap);
  
  //Now The Dead channel with neighbouring channel on will Be activated
  //in order to have the cluster properly computed at the clustering level.
  //theStripElectronicsSim->StripVFats contain already a list of dead
  //thePadElectronicsSim->PadVFats contain already a list of dead
  //T2PadDigis in T2PadDigiCollection will be updated.
  //T2StripDigis in T2StripDigiCollection will be updated.
  
  StripVFatsdig= &(theStripElectronicsSim->StripVFats);
  
  PadVFatsdig= &(thePadElectronicsSim->PadVFats);
  
} // process

void T2Digitizer::createStripChargeMap() {
  for(int i=0; i<2; i++) {
    for(int j=0; j<2; j++) {
      for(int k=0; k<5; k++) {
		for(int l=0; l<2; l++) {
		  int* tmp = new int[2*256];
		  int detId = T2DetId::calculateRawId(i, j, k, l);
		  stripChargeMap[detId] = tmp;
		}
      }
    }
  }
} // createStripChargeMap


void T2Digitizer::resetStripChargeMap() {
  std::map<int, int*>::iterator i;
  for (i=stripChargeMap.begin(); i != stripChargeMap.end(); ++i) {
    for (int j=0; j<2*256; j++) {
      (*i).second[j] = 0;
    }  
  }  
} // resetStripChargeMap

void T2Digitizer::deleteStripChargeMap() {
  std::map<int, int*>::iterator i;
  for (i=stripChargeMap.begin(); i != stripChargeMap.end(); ++i) {
    delete[] (*i).second;
  }  
  stripChargeMap.clear();  
} // resetStripChargeMap


void T2Digitizer::printoutStuff() {

  std::map<int, int*>::iterator i;
  for (i=stripChargeMap.begin(); i != stripChargeMap.end(); ++i) {
    std::cout<<std::endl;
    std::cout << (*i).first << std::endl;
    for(int k=0;k<2; k++){
      for(int l=0;l<256; l++) {
    	  std::cout << " " << (*i).second[l+k*256];
      }
    }
  }
} // printoutStuff

void T2Digitizer::createPadChargeMap() {
  
  for(int i=0; i<2; i++) {
    for(int j=0; j<2; j++) {
      for(int k=0; k<5; k++) {
		for(int l=0; l<2; l++) {
		  int* tmp = new int[65*24];
		  int detId = T2DetId::calculateRawId(i, j, k, l);
		  padChargeMap[detId] = tmp;
		}
      }
    }
  }
} // createPadChargeMap

void T2Digitizer::resetPadChargeMap() {
  std::map<int, int*>::iterator it;
  for (it=padChargeMap.begin(); it != padChargeMap.end(); ++it) {
    for (int j=0; j<65*24; j++) {
      (*it).second[j] = 0;
    }  
  }  
} // resetPadChargeMap

void T2Digitizer::deletePadChargeMap() {
  std::map<int, int*>::iterator it;
  for (it=padChargeMap.begin(); it != padChargeMap.end(); ++it) {
    delete[] (*it).second;
  }  
  padChargeMap.clear();  
} // resetPadChargeMap

unsigned int T2Digitizer::RawtoSymb(uint32_t thedet)
{
  T2DetId converter;
  unsigned int pl=converter.plane(thedet);
  unsigned int pls=converter.planeSide(thedet);
  unsigned int ht=converter.halfTelescope(thedet);
  unsigned int arm=converter.arm(thedet);
  unsigned int symbolic=pl*2+pls+ht*10+20*arm;
  return symbolic;
}


void T2Digitizer::HitNoisePosition(double &x, double &y,double RNoiseSlope)
{
  double StripMinR = 42.46;
  double StripMaxR =StripMinR +256*0.4; 
  bool found=false;
  //I Want large-R hits (dente di sega) and quite close to x initial 

  while(found==false){
     double ei=CLHEP::RandFlat::shoot();
     if(((sqrt(x*x+y*y))/StripMaxR)>ei){
    	 found=true;
     }else{
		 double ux=x*(-1.0);
		 double uy=y*(-1.0);
		 while((ux*x<0)||(uy*y<0)||((ux*ux+uy*uy)>StripMaxR*StripMaxR)||((ux*ux+uy*uy)<StripMinR*StripMinR)){
			 ux=CLHEP::RandGaussQ::shoot(x,2*StripMinR);
			 uy=CLHEP::RandGaussQ::shoot(y,2*StripMinR);
		 }
		 x=ux;
		 y=uy;
     }
  }
}





bool T2Digitizer::PrimaryTrackCondition(SimTrack atrk,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, const HepMC::GenEvent* evt,int &barcodeMother){


  bool primarycondition=false;
  
  int trackId=atrk.trackId();
  //SimTrack is stored and the generating vertex is at IP.
  
  bool isTrackPrimaryIP=false;
  if(fabs(atrk.charge())>0){
    isTrackPrimaryIP = isTrackPrimary(trackId, theSimTracks, theSimVertices);
  }
    
    //SimTrack PythiaIPChParticleTrack=(*ittrack);      
  int thePidGenerator=0;
  int thePidGeneratorOfThisTrk=0;
  
  thePidGenerator=GetTrkOlderMotherPid(atrk,theSimVertices,theSimTracks,evt,barcodeMother,thePidGeneratorOfThisTrk);

  if(thePidGeneratorOfThisTrk==0) //It means that Geant trk was not associated to a pythia particle
    thePidGeneratorOfThisTrk=atrk.type();

  //bool isTrackPrimaryIP = isTrackPrimary(trackId, theSimTracks, theSimVertices);     	        

  //if(verbosity>0)
  bool alloweddecayparticle=false;
  if(std::find(AllowedDecayToCount.begin(),AllowedDecayToCount.end(),fabs(thePidGenerator))!=AllowedDecayToCount.end())
    alloweddecayparticle=true;

  //I'm Excluding 211 also because I assume that a nonprimary older mother 211 have done Nuclear interaction
  if(fabs(atrk.charge())>0)
    {
      //&&(fabs(thePidGeneratorOfThisTrk)!=11)&&(fabs(thePidGeneratorOfThisTrk)!=22)
      if((isTrackPrimaryIP)||((isTrackPrimaryIP==false)&&(fabs(thePidGenerator)!=211)&&(thePidGeneratorOfThisTrk!=0)&&(alloweddecayparticle==true)))
	primarycondition=true;
    }

  //   std::cout<<"GeantTrkID:"<<trackId<<" PID_Mother:"<<thePidGenerator<<" Track PID: "<<thePidGeneratorOfThisTrk<<" Direct Primary_Condition: "<<isTrackPrimaryIP<<" Primary_Condition: "<<primarycondition<<std::endl;
  //Just a patch.
  primarycondition=false;
  if(isTrackPrimaryIP)
    primarycondition=true;
  return primarycondition;

}



bool T2Digitizer::getTrack(const unsigned int trackId,  const std::auto_ptr<edm::SimTrackContainer>& theSimTracks, SimTrack &track)
{
  if(trackId < 1){
    return false;
  }

  //Guesss
  if (trackId <= theSimTracks.get()->size()){
    track = (*theSimTracks.get())[trackId-1];
    if (track.trackId() == trackId){
      return true;
    }
  }

  //Iterate all
  for(edm::SimTrackContainer::const_iterator ittrack = theSimTracks->begin();ittrack != theSimTracks->end(); ++ittrack){     
    if (ittrack->trackId() == trackId){
      track  = (*ittrack);
      return true;
    }
  }
  
  return false;  
}

/*
  Returns vertex corresponding to vertexId 
  Returns true if found, false otherwise.
*/

bool
T2Digitizer::getVertex(const  int vertexId,  const std::auto_ptr<edm::SimVertexContainer>& theSimVertices, SimVertex &vertex)
{
  if (vertexId < (int) theSimVertices.get()->size() && vertexId >=0){
    vertex = (*theSimVertices.get())[vertexId];
    return true;
  }
  return false;  
}


bool T2Digitizer::isTrackPrimary(const int trackId,
					 const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,
					 const std::auto_ptr<edm::SimVertexContainer>& theSimVertices)
{
  SimTrack track;
  if (getTrack(trackId,theSimTracks,track) == true){
    SimVertex vertex;    
    if (getVertex(track.vertIndex(), theSimVertices,vertex) == true){
      if (vertex.parentIndex() < 0){
	return true;
      }
    }
  }
  return false;
}






int T2Digitizer::GetTrkOlderMotherPid(SimTrack aTrack,const std::auto_ptr<edm::SimVertexContainer>& theSimVertices,
						const std::auto_ptr<edm::SimTrackContainer>& theSimTracks,const HepMC::GenEvent* evt, int &motherbarcode,int &thePidGeneratorOfThisTrk)
{
  
  /////////////////////////////////
  int OlderMotherParticle_Id=0; //This will be returned.
  /////////////////////////////////

  unsigned int count=0;
  SimTrack DirectMotherTrk=aTrack;
  int DirectMotherTrk_Id=aTrack.trackId();
  
 
  for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
    if((*p)->barcode()==aTrack.genpartIndex())
      {
	// thePidGenerator=(*p)->pdg_id();
	((*p)->barcode());
	motherbarcode=((*p)->barcode());
	thePidGeneratorOfThisTrk=(*p)->pdg_id();	
	OlderMotherParticle_Id=(*p)->pdg_id();		
      }
  }
 //--  std::cout<<"GetTrkOlderMotherPid begin. Trk Id: "<<DirectMotherTrk_Id<<" Pythia pdg:"<<OlderMotherParticle_Id<<" Trk_PID:"<<(DirectMotherTrk.type())<<std::endl;

 //std::cout<<"Try to find mother for Particle_ID: "<<(aTrack.type())<<" with Trk_ID:"<<DirectMotherTrk_Id<<std::endl;

 while(DirectMotherTrk_Id!=(-1))//loop until you dont' find the stable mother or the source
   {
     count++;     

     SimVertex simuvertex;
 
     //getTrack: Returns track corresponding to DirectMotherTrk_Id. Returns true if found, false otherwise.
  
     if(getTrack(DirectMotherTrk_Id,theSimTracks,DirectMotherTrk) == true)
       {
	  //-- std::cout<<"Retrieved geant Trk with Particle_PID: "<<(DirectMotherTrk.type())<<" and Trk_ID:"<<DirectMotherTrk_Id<<std::endl;
	 OlderMotherParticle_Id=DirectMotherTrk.type();
	 //std::cout<<"Mother Trk id= "<<DirectMotherTrk_Id<<std::endl;
	 
	 //This look only if parentVtx is primary (Index<0) 
	 if(isTrackPrimary(DirectMotherTrk_Id, theSimTracks, theSimVertices)==false)
	   {	    
	     //-- std::cout<<"Not a primary track. Try to find originating vertex..."<<std::endl;
	     //Take the vtx generating the trk.
	     if(getVertex(DirectMotherTrk.vertIndex(), theSimVertices,simuvertex) == true) 
	       {	    
		 
		//--  std::cout<<"Originating vertex found. Try to found the ancestors of this track:"<<std::endl;
		 // Note that DirectMotherTrk.genpartIndex() is the barcode, 
		 // not the index within GenParticleCollection, so I have to search the particle
		  

		 //THIS IS NEVER WORKING ACTUALLY
		  for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
		    if((*p)->barcode()==DirectMotherTrk.genpartIndex())
		      {
			(*p)->pdg_id();
		   
			//-- std::cout<<"Pythia Ancestors Particle-IDs of the track: ";
			HepMC::GenVertex::particle_iterator ancestor;
			ancestor = (*p)->production_vertex()->particles_begin(HepMC::parents);
			for(;ancestor != (*p)->production_vertex()->particles_end(HepMC::parents); ++ancestor) {
			 
			  //std::cout<<"  "<<((*ancestor)->pdg_id());
			  //motherbarcode to return is saved.
			  motherbarcode=(*ancestor)->barcode();
			  
			  if((std::find(knowPart.begin(),knowPart.end(),fabs((*ancestor)->pdg_id()))!=knowPart.end()))
			    {
			      ((*ancestor)->barcode());
			      OlderMotherParticle_Id=(*ancestor)->pdg_id();
			      //--std::cout<<"Ancestor is a Known particle: Barcode:"<<OlderMotherId<<" PID:"<<OlderMotherParticle_Id<<std::endl;
			    }

			}
			//-- std::cout<<"OlderMotherParticle_Id is now: "<<	OlderMotherParticle_Id<<std::endl;		   		   		   
		      }
		    //std::cout<<"Begin-2"<<motherbarcode<<std::endl;
		  }
		 
		  //Find the mother track originating this vtx 
		  //This instruction allows to go back in the three to found the origin
		   //THIS IS THE ONLY IMPORTANT INSTRUCTION
		  DirectMotherTrk_Id=simuvertex.parentIndex(); 		  
		  //--std::cout<<"Mother Trk Id originating this vtx:"<<DirectMotherTrk_Id<<std::endl;
		
	       }
	    //--  else
	        //-- std::cout<<"Originating vertex NOT found:"<<std::endl;
	   }
	 else //Track is primary
	   { 
	    //--  std::cout<<"isTrackPrimary==true"<<std::endl;
	      //In this case DirectMotherTrk_Id should be <0 .
	     // std::cout<<"isTrackPrimary==true so Mother Trk is < 0, I want to store this info also. Track printing:"<<std::endl;
	     // Note that DirectMotherTrk.genpartIndex() is the barcode, not the index within GenParticleCollection, so I have to search the particle
	     for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
	       if((*p)->barcode()==DirectMotherTrk.genpartIndex())
		 {
		   (*p)->pdg_id();
		   
		   //   std::cout<<"Pythia Ancestors Particle-ID : ";
		   HepMC::GenVertex::particle_iterator ancestor;
		   ancestor = (*p)->production_vertex()->particles_begin(HepMC::parents);
		   for(;ancestor != (*p)->production_vertex()->particles_end(HepMC::parents); ++ancestor) {
		     
		     motherbarcode=(*ancestor)->barcode();
		     if((std::find(knowPart.begin(),knowPart.end(),fabs((*ancestor)->pdg_id()))!=knowPart.end()))
		       {
			((*ancestor)->barcode());
			 OlderMotherParticle_Id=(*ancestor)->pdg_id();
			 motherbarcode=(*ancestor)->barcode();
		       }
		     //otherwise keep what you already have, since pythia is not giving any true particle
		   }
		   		   		   		   
		 }
	       //std::cout<<"Begin-2"<<motherbarcode<<std::endl;
	     }
	     
	     //std::cout<<"Exit condition active since Mother Trk was primary"<<std::endl;
	     DirectMotherTrk_Id=-1;//Exit condition active
	   }
       }
     else
      DirectMotherTrk_Id=-1;//Exit condition active  but you have not recovered the trk  
     
   }//while end

 if(count==0) //The input trkId was already primary
   {
      for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
	if((*p)->barcode()==aTrack.genpartIndex())
	  {
	    (*p)->pdg_id();
	  ((*p)->barcode());
	    motherbarcode=((*p)->barcode());
	    OlderMotherParticle_Id=(*p)->pdg_id();		
	  }
      }

   }

 // std::cout<<"GetTrkOlderMotherPid Return Trk_Id:"<<OlderMotherId<<std::endl;
 return OlderMotherParticle_Id;

}

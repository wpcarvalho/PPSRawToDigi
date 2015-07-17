#ifndef HelpfulWatchers_SimTracer_h
#define HelpfulWatchers_SimTracer_h
// -*- C++ -*-
//
// Package:     HelpfulWatchers
// Class  :     SimTracer
// 
/**\class SimTracer SimTracer.h SimG4Core/HelpfulWatchers/interface/SimTracer.h

 Description: Prints a message for each Oscar signal

 Usage:
    <usage>

*/
//
// Original Author:  
//         Created:  Tue Nov 22 16:41:33 EST 2005
//

// system include files
#include <iostream>

// user include files
#include "SimG4Core/Notification/interface/BeginOfTrack.h"
#include "SimG4Core/Notification/interface/Observer.h"
#include "SimG4Core/Watcher/interface/SimWatcher.h"
#include "SimG4Core/Notification/interface/TrackInformation.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4String.hh"
#include "G4VProcess.hh"

// forward declarations
class DDDWorld;
class BeginOfJob;
class BeginOfRun;
class BeginOfEvent;
class BeginOfTrack;
class G4Step;

class EndOfRun;
class EndOfEvent;
class EndOfTrack;

#define OBSERVES(type) public Observer<const type*>
#define UPDATE(type) void update(const type*) { std::cout <<"++ signal " #type<<std::endl; }
class SimTracer : public SimWatcher, 
OBSERVES(DDDWorld),
OBSERVES(BeginOfJob),
OBSERVES(BeginOfRun),
OBSERVES(BeginOfEvent),
OBSERVES(BeginOfTrack),
OBSERVES(G4Step),
OBSERVES(EndOfRun),
OBSERVES(EndOfEvent),
OBSERVES(EndOfTrack)
{

   public:
   SimTracer(const edm::ParameterSet& pSet) :
   m_verbose(pSet.getUntrackedParameter<bool>("verbose",true)) {
   }
     //virtual ~SimTracer();

      // ---------- const member functions ---------------------

      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------
UPDATE(DDDWorld)
UPDATE(BeginOfJob)
UPDATE(BeginOfRun)
UPDATE(BeginOfEvent)
    
   void update(const G4Step* iStep) {
   std::cout <<"++ signal G4Step " ;
   if(m_verbose) {
      const G4StepPoint* pre = iStep->GetPreStepPoint();
      const G4StepPoint* post = iStep->GetPostStepPoint();
      const G4ThreeVector pre_pos = pre->GetPosition();
      const G4ThreeVector pos = post->GetPosition();
      const G4ThreeVector mom_pre = pre->GetMomentum();
      const G4ThreeVector mom_post = post->GetMomentum();
      const G4Material* mat_pre = pre->GetMaterial();
      const G4Material* mat_post = post->GetMaterial();

      std::cout.precision(10);
      if(pre->GetPhysicalVolume())
      {
        std::cout << pre->GetPhysicalVolume()->GetName();
      }
      if(mat_pre)
      {
        std::cout<<" Material="<<mat_pre->GetName()<<std::endl;
      }
      std::cout << "pre_pos=( "<<pre_pos.x()<<","<<pre_pos.y()<<","<<pre_pos.z()<<") ";
      std::cout << "pre_mom=( "<<mom_pre.x()<<","<<mom_pre.y()<<","<<mom_pre.z()<<") "<<std::endl;
      if(post->GetPhysicalVolume())
      {
        std::cout << post->GetPhysicalVolume()->GetName();
      }
      if(mat_post)
      {
        std::cout<<" Material="<<mat_post->GetName()<<std::endl;
      }
      std::cout << "post_pos=( "<<pos.x()<<","<<pos.y()<<","<<pos.z()<<") ";
      std::cout << "post_mom=( "<<mom_post.x()<<","<<mom_post.y()<<","<<mom_post.z()<<")"<<std::endl;
    }
    std::cout <<std::endl;
  }
  //UPDATE(BeginOfTrack)

  void update(const BeginOfTrack* beg_of_track) {
    std::cout <<"++ signal BeginOfTrack " ;
    if(m_verbose) {
      G4ParticleDefinition* pdef = (*beg_of_track)()->GetDefinition();
      int tr_id = (*beg_of_track)()->GetTrackID();
      int par_id = (*beg_of_track)()->GetParentID();
      const G4String part_name = pdef->GetParticleName();
      G4String proc_name="";
      G4ProcessType proc_type = fNotDefined;
      if((*beg_of_track)()->GetCreatorProcess())
      {
        proc_name=(*beg_of_track)()->GetCreatorProcess()->GetProcessName();
        proc_type = (*beg_of_track)()->GetCreatorProcess()->GetProcessType();
      }
      TrackInformation* info
      = dynamic_cast<TrackInformation*>( (*beg_of_track)()->GetUserInformation() );

      std::cout<<" track id:"<<tr_id<<" parent id:"<<par_id<<" particle:"<<part_name;
      std::cout<<" tot.energy[GeV]:"<<(*beg_of_track)()->GetTotalEnergy()/GeV<<std::endl;
      std::cout<<" prim vertex:"<<(*beg_of_track)()->GetVertexPosition()<<
      " vert.mom.dir:"<<(*beg_of_track)()->GetVertexMomentumDirection()<<
      " creat.p.type:" <<proc_type <<
      " creat.p.name:" << proc_name <<std::endl;
      if(info && info->isPrimary())
      {
        std::cout<<"Primary particle"<<std::endl;
      }
   }
   std::cout <<std::endl;
}
//UPDATE(G4Step)
UPDATE(EndOfRun)
UPDATE(EndOfEvent)
UPDATE(EndOfTrack)

   private:
     //SimTracer(const SimTracer&); // stop default

     //const SimTracer& operator=(const SimTracer&); // stop default

     // ---------- member data --------------------------------
     bool m_verbose;
};


#endif

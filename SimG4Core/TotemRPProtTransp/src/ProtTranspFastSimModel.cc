/****************************************************************************
*
* This is a part of the TOTEM offline software.
* Authors:
*    Hubert Niewiadomski
*    Jan Kašpar (jan.kaspar@gmail.com)
*
* $$RCSfile: ProtTranspFastSimModel.cc,v $: $
* $Revision: 1.5.2.3 $
* $Date: 2009/11/16 16:54:43 $
*
****************************************************************************/

//#define G4V7
#define DEBUG 0

#include "SimG4Core/TotemRPProtTransp/interface/ProtTranspFastSimModel.h"

#include <iostream>
#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleMomentum.hh"
#include "G4Proton.hh"
#include "G4Track.hh"

using namespace std;

const double ProtTranspFastSimModel::beampipe_aperture_radius = 0.040;  //in m

ProtTranspFastSimModel::ProtTranspFastSimModel(G4String ModelVolumeName, 
    G4Envelope* envelope, const LHCOpticsApproximator &approx, double zin, double zout,
    bool verbosity)
     : G4VFastSimulationModel(ModelVolumeName, envelope), modelName_(ModelVolumeName),
     approximator_(approx), zin_(zin), zout_(zout),
    verbosity_(verbosity)
{ 
  if(verbosity_)
    std::cout<<"ProtTranspFastSimModel created for volume "<<ModelVolumeName<<std::endl;
    
  if(ModelVolumeName == "Beam_IP_150_R")
    transport_region_ = BEAM_IP_150_R;
  else if (ModelVolumeName == "Beam_IP_150_L")
    transport_region_ = BEAM_IP_150_L;
  else
    transport_region_ = NOT_DEF;
}


void ProtTranspFastSimModel::DoIt(const G4FastTrack& track, G4FastStep& step)
{
  if(verbosity_)
  {
    std::cout<<">> ProtTranspFastSimModel::DoIt"<<std::endl;
    cout << "\tstep info: " << endl;
    step.DumpInfo();
  }

  step.Initialize(track);
  //step.ForceSteppingHitInvocation();
  step.KillPrimaryTrack();
  
  G4ThreeVector in_loc_pos = track.GetPrimaryTrackLocalPosition();
  G4ThreeVector in_loc_mom = track.GetPrimaryTrackLocalMomentum();
  const G4AffineTransform* loc_glob_tran = track.GetInverseAffineTransformation();
  
  G4ThreeVector in_glob_pos = loc_glob_tran->TransformPoint(in_loc_pos);
  G4ThreeVector in_glob_mom = loc_glob_tran->TransformPoint(in_loc_mom);
  
  if(verbosity_)
  {
    std::cout<<"\tin_loc_pos:"<<in_loc_pos<<"\n\tin_loc_mom:"<<in_loc_mom<<
      "\n\tin_glob_pos:"<<in_glob_pos<<"\n\tin_glob_mom:"<<in_glob_mom<<std::endl;
  }

  if(IsBeamPipeProton(track))
  {
    if(verbosity_)
      std::cout<<"!!!! Is a BeamPipe Proton"<<std::endl;
      
    if(!Transport(track, step))
    {
      if(verbosity_)
        std::cout<<"!!!! Transport failed, aperture limits"<<std::endl;
      KillParticleAndSecondaries(track, step);
    }
  }
  else
  {
    if(verbosity_)
      std::cout<<"!!!! Is not a BeamPipe Proton, to be killed"<<std::endl;
    KillParticleAndSecondaries(track, step);
  }
//  double step_length = TMath::Abs(zout_-zin_)*m;
//  double step_time = step_length/3e8*second;
//  double step_proper_time = step_time/7777;
//  const G4Track *prim_track = track.GetPrimaryTrack();
//  
//  step.ProposePrimaryTrackPathLength(step_length + prim_track->GetTrackLength());
//  step.ProposePrimaryTrackFinalTime(step_time + prim_track->GetGlobalTime());
//  step.ProposePrimaryTrackFinalProperTime(step_proper_time + prim_track->GetProperTime());
//  step.DumpInfo();
}


bool ProtTranspFastSimModel::IsBeamPipeProton(const G4FastTrack& track)
{
  G4ThreeVector pos = track.GetPrimaryTrackLocalPosition();
  G4ThreeVector dir = track.GetPrimaryTrackLocalDirection();
  
  //if is a proton
  if(track.GetPrimaryTrack()->GetDefinition()->GetPDGEncoding()!=2212)
    return false;
  
  //if inside of the beampipe
  if(pos.x()/m*pos.x()/m + pos.y()/m*pos.y()/m > beampipe_aperture_radius*beampipe_aperture_radius)
    return false;
  
  //if the momentum is in the correct direction
  if(dir.z()<=0 && (transport_region_ == BEAM_IP_150_R))
    return false;
    
  if(dir.z()>=0 && (transport_region_ == BEAM_IP_150_L))
    return false;
  
  return true;
}


void ProtTranspFastSimModel::KillParticleAndSecondaries(const G4FastTrack& fastTrack, 
      G4FastStep& fastStep)
{
  fastStep.KillPrimaryTrack();
}


bool ProtTranspFastSimModel::Transport(const G4FastTrack& track, G4FastStep& step)
{
  if(verbosity_) {
    cout.precision(25);
    cout<<"===== BEGIN Transport " << GetName() << "=================="<<endl;
  }

  // the global and local coordinates in case of the beam parametrized area are the same
  G4ThreeVector in_pos = Unsmear_z_position(track);
  G4ThreeVector in_mom = track.GetPrimaryTrackLocalMomentum();
  //in_mom = track.GetPrimaryTrack()->GetDynamicParticle()->GetPrimaryParticle()->GetMomentum()

  if(verbosity_) {
    cout << "input" << endl;
    cout << "\tposition: " << in_pos << endl;
    cout << "\tmomentum: " << in_mom << endl;
  }

  double in_position[3];
  double in_momentum[3];
  double out_position[3];
  double out_momentum[3];
  
  in_position[0] = in_pos.x()/meter;
  in_position[1] = in_pos.y()/meter;
  in_position[2] = zin_;
  in_momentum[0] = in_mom.x()/GeV;
  in_momentum[1] = in_mom.y()/GeV;
  in_momentum[2] = in_mom.z()/GeV;
  
  if(verbosity_) {
    std::cout<<"before transport"<<std::endl;
    std::cout<<"\tposition: " << in_position[0] << ", " << in_position[1] << ", " << in_position[2] <<std::endl;
    std::cout<<"\tmomentum: " << in_momentum[0] << ", " << in_momentum[1] << ", " << in_momentum[2] <<std::endl;
    
    /* TODO, E might not be 7TeV !!!
    std::cout<<"\tin state: MADThx="<<in_momentum[0]/7000
    <<", MADThy="<<in_momentum[1]/7000
    <<", xi="<<(fabs(in_momentum[2])-7000)/7000.0<<std::endl;
    */
  }

  bool tracked = approximator_.Transport_m_GeV(in_position, in_momentum, 
      out_position, out_momentum, true, zout_-zin_);

  if(!tracked) {
    if(verbosity_) {
      cout << "* proton not tracked" << endl;
    cout<<"===== END Transport " << GetName() << "===================="<<endl;
    }
    return false;
  }

  if(verbosity_) {
    std::cout<<"after transport"<<std::endl;
    std::cout<<"\tposition: " << out_position[0] << ", " << out_position[1] << ", " << out_position[2] <<std::endl;
    std::cout<<"\tmomentum: " << out_momentum[0] << ", " << out_momentum[1] << ", " << out_momentum[2] <<std::endl;
  }
      
  if(out_position[0]*out_position[0]+out_position[1]*out_position[1]>
      beampipe_aperture_radius*beampipe_aperture_radius) {
    if(verbosity_) {
      cout << "* proton ouside beampipe" << endl;
      cout<<"===== END Transport " << GetName() << "===================="<<endl;
    }
    return false;
  }
  
  G4ThreeVector out_pos(out_position[0]*meter, out_position[1]*meter, out_position[2]*meter);
  G4ThreeVector out_mom(out_momentum[0]*GeV, out_momentum[1]*GeV, out_momentum[2]*GeV);

  if(verbosity_) {
    cout << "output" << endl;
    cout << "\tposition: " << out_pos << endl;
    cout << "\tmomentum: " << out_mom << endl;
  }


//  step.ProposePrimaryTrackFinalPosition(out_pos, false);
//  step.ProposePrimaryTrackFinalMomentumDirection(out_mom, false);
  
  //---------------------------
  // Secondary:
  //   Adds one "secondary":
  //
  //---------------------------
  // -- First, user has to say how many secondaries will be created:
  step.SetNumberOfSecondaryTracks(1);
 
  //------------------------
  // -- Build the secondary:
  //------------------------
  // -- direction:
  G4ParticleMomentum direction(out_mom.unit());
 
  // -- dynamics (Note that many constructors exists for G4DynamicParticle
  // -- see prototype/particle+matter/particles/management/include/G4DynamicParticle.hh)
  G4DynamicParticle dynamique(G4Proton::ProtonDefinition(),
                              direction,
                              track.GetPrimaryTrack()->GetKineticEnergy());
  // -- position:   
  //------------------------------------
  //-- Creation of the secondary Track:
  //------------------------------------
  /*G4Track *out_track = */step.CreateSecondaryTrack(dynamique, out_pos,
                       track.GetPrimaryTrack()->GetGlobalTime());
  
//  G4int parent_id = track.GetPrimaryTrack()->GetParentID();
//  out_track->SetParentID(parent_id);
  
  if(verbosity_) {
    std::cout<<"* Proton transported successfully"<<std::endl;
    cout<<"===== END Transport " << GetName() << "===================="<<endl;
  }
  
  return true;
}


G4ThreeVector ProtTranspFastSimModel::Unsmear_z_position(const G4FastTrack& track)
{
  //if(verbosity_) std::cout<<">> ProtTranspFastSimModel::Unsmear_z_position"<<std::endl;
  
  G4ThreeVector registered_pos;
  if(transport_region_ == BEAM_IP_150_R || transport_region_ == BEAM_IP_150_L)
  { 
    //if(verbosity_) std::cout<<"\tBEAM_IP_150_R || BEAM_IP_150_L"<<std::endl;
    registered_pos = track.GetPrimaryTrack()->GetVertexPosition();
  }
  else
  {
    //if(verbosity_) std::cout<<"\tNOT A BEAM_IP_150_R || BEAM_IP_150_L"<<std::endl;
    registered_pos = track.GetPrimaryTrackLocalPosition();
  }
  
  double dist = registered_pos.z() - zin_*m;
  //if(verbosity_) std::cout<<"\tz distance to normalize position: "<<dist<<std::endl;
  
  G4ThreeVector dir = track.GetPrimaryTrackLocalDirection();
  dir = dir/dir.z();
  G4ThreeVector pos_for_param = registered_pos - dir*dist;
  
  return pos_for_param;
}


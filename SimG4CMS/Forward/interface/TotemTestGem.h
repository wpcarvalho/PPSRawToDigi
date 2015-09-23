#ifndef Forward_TotemTestGem_h
#define Forward_TotemTestGem_h 1
// -*- C++ -*-
//
// Package:     Forward
// Class  :     TotemTestGem
//
/**\class TotemTestGem TotemTestGem.h SimG4CMS/Forward/interface/TotemTestGem.h
 
 Description: Manages Root file creation for Totem Tests
 
 Usage:
    Used in testing Totem simulation
 
*/
//
// Original Author: 
//         Created:  Tue May 16 10:14:34 CEST 2006
//
 
// system include files
#include <iostream>
#include <memory>
#include <vector>
#include <string>

// user include files
#include "SimG4Core/Notification/interface/Observer.h"
#include "SimG4Core/Notification/interface/BeginOfJob.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/Notification/interface/EndOfEvent.h"
#include "SimG4Core/Watcher/interface/SimWatcher.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "SimG4CMS/Forward/interface/TotemTestHistoManager.h"
#include "SimG4CMS/Forward/interface/TotemTestHistoClass.h"
#include "SimG4CMS/Forward/interface/TotemG4Hit.h"
#include "SimG4CMS/Forward/interface/TotemHisto.h"

class G4Step;

class TotemTestGem : public SimWatcher,
                     public Observer<const BeginOfJob *>,
                     public Observer<const BeginOfEvent *>,
                     public Observer<const EndOfEvent *>,
                     public Observer<const G4Step *> {

public:
  TotemTestGem(const edm::ParameterSet &p);
  virtual ~TotemTestGem();

  virtual void produce(edm::Event&, const edm::EventSetup&);

private:
  // observer classes
  void update(const BeginOfJob * job);
  void update(const BeginOfEvent * evt);
  void update(const EndOfEvent * evt);
  void update(const G4Step * step);

  void clear();
  void fillEvent(TotemTestHistoClass&);

private:
  edm::ParameterSet                       parameters;

  //Keep parameters to instantiate TotemTestHistoManager later
  std::string                             fileName;
  std::vector<std::string>                names;

  // Private Tuples
  std::auto_ptr<TotemTestHistoManager>    tuplesManager;
  TotemTestHistoClass*                    tuples;

  TotemHisto*                             histos;
  std::string                             fileNameOld;
  std::vector<TotemG4Hit*>                hits;
  int                                     evtnum;
};

#endif

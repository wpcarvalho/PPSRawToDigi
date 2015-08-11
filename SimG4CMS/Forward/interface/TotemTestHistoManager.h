#ifndef Forward_TotemTestHistoManager_h
#define Forward_TotemTestHistoManager_h 1
// -*- C++ -*-
//
// Package:     Forward
// Class  :     TotemTestHistoManager
//
/**\class TotemTestHistoManager TotemTestHistoManager.h SimG4CMS/Forward/interface/TotemTestHistoManager.h
 
 Description: Manages Root file creation for Totem Tests
 
 Usage:
    Used in testing Totem simulation
 
*/
//
// Original Author: 
//         Created:  Tue May 16 10:14:34 CEST 2006
// $Id: TotemTestHistoManager.h,v 1.1.1.1 2007/05/16 15:44:30 hniewiad Exp $
//
 
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"

#include <string>

class TotemTestHistoClass;

class TotemTestHistoManager {

public: 

  TotemTestHistoManager(const std::string &);
  virtual ~TotemTestHistoManager();

  void fillTree(TotemTestHistoClass *  histos);

private:
  edm::Service<TFileService> fs;
  TTree                     *tree;
  TotemTestHistoClass       *h;
  int                       kount;
};

#endif

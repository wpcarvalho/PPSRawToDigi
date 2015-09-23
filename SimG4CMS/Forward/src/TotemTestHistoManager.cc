// -*- C++ -*-
//
// Package:     Forward
// Class  :     TotemTestHistoManager
//
// Implementation:
//     <Notes on implementation>
//
// Original Author: 
//         Created:  Tue May 16 10:14:34 CEST 2006
// $Id: TotemTestHistoManager.cc,v 1.1.1.1.6.1 2009/07/13 14:36:52 jkaspar Exp $
//

#include <iostream>
#include <cmath>

#include "SimG4CMS/Forward/interface/TotemTestHistoManager.h"
#include "SimG4CMS/Forward/interface/TotemTestHistoClass.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/PluginManager/interface/PluginManager.h"

TotemTestHistoManager::TotemTestHistoManager(const std::string & file):
	tree(0), h(0), kount(0) {
	if (fs.isAvailable()) {
		h    = new TotemTestHistoClass();
		tree = fs->make<TTree>("ForwardSim", "ForwardSim");
		tree->SetAutoSave(10000);
		tree->Branch("TotemTestHisto", "TotemTestHistoClass", &h);
		edm::LogInfo("ForwardSim") << "TotemTestHistoManager:===>>>  Book the Tree";
	} else
		edm::LogInfo("ForwardSim") << "TotemTestHistoManager:===>>> No file provided";
}

TotemTestHistoManager::~TotemTestHistoManager() {
	edm::LogInfo("ForwardSim") << "============================================="
	                        << "========================================\n"
	                        << "=== TotemTestHistoManager: Start writing user "
	                        << "histograms after " << kount << " events ";
	if (h) delete h;
}

void TotemTestHistoManager::fillTree(TotemTestHistoClass *  histos) {
	kount++;
	LogDebug("ForwardSim") << "TotemTestHistoManager: tree pointer for " << kount
	                    << " = " << histos;
	if (tree) {
	    h = histos;
	    tree->Fill();
	}
}

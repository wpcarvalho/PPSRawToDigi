// -*- C++ -*-
//
// Package:    L1TriggerTotem
// Class:      RPTriggerAnalyzer
// 
/**\class RPTriggerAnalyzer RPTriggerAnalyzer.cc L1TriggerTotem/RPTriggerAnalyzer/src/RPTriggerAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jiri Prochazka
//         Created:  Mon Mar  1 13:40:46 CET 2010
// $Id$
//
//
#ifndef _L1TriggerTotemRPTriggerAnalyzer_H_
#define _L1TriggerTotemRPTriggerAnalyzer_H_

// system 
#include <memory>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

// ROOT 
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"

// CMSSW 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
//#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/DetSetVector.h"
//#include "DataFormats/Common/interface/Wrapper.h"
//#include "DataFormats/Common/interface/TriggerResults.h"

//#include "FWCore/Utilities/interface/InputTag.h"

// TOTEM 
//#include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDetTrigger.h"
#include "DataFormats/TotemL1Trigger/interface/RPCCBits.h"

#include "L1TriggerTotem/CoincidenceChip/interface/RPTriggerAnalyzerInfoCollector.h"
#include "L1TriggerTotem/CoincidenceChip/interface/PotCollection.h"



class RPTriggerAnalyzer : public edm::EDAnalyzer {
    public:
        explicit RPTriggerAnalyzer(const edm::ParameterSet&);
        ~RPTriggerAnalyzer();
        virtual std::string ClassName() const { return "RPTriggerAnalyzer"; }

    private:
        virtual void beginRun(edm::Run const&, edm::EventSetup const&);
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob();

        // void PrintText(double val, const TString& separator = " & ") const;

        // ----------member data ---------------------------
        unsigned int verbosity;


        
        RPTriggerAnalyzerInfoCollector fInfoCollector;
        PotCollection  fPotCollection;
        const edm::ParameterSet fConfig;
        edm::InputTag detTriggerLabel;
  //      boost::ptr_vector<TH1D> fHist1D;
  //      boost::ptr_vector<TH2D> fHist2D;
};




#endif

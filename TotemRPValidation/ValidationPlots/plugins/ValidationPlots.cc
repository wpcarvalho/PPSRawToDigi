/*
 * ValidationPlots.h
 *
 *  Created on: Oct 29, 2008
 *      Author: Leszek Grzanka
 */

#include <stdio.h>
#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPTrackCandidateCollection.h"
#include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDetTrigger.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "DataFormats/TotemRPDataTypes/interface/RPRecoHit.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecoElasticEvent.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TSystem.h"
#include "TRint.h"
#include "TStyle.h"

#include "TotemRPValidation/ValidationPlots/interface/ValidationPlots.h"

//----------------------------------------------------------------------------------------------------

ValidationPlots::ValidationPlots(const edm::ParameterSet& conf)
{
  using namespace std;
  verbosity = conf.getUntrackedParameter<unsigned int>("verbosity", 0);
  outputFile = conf.getParameter<std::string>("outputFile");
  elasticHists = new ElasticRecoValLibrary(conf);
  hitDistributions = new HitDistributionsLibrary(conf);
}

//----------------------------------------------------------------------------------------------------

ValidationPlots::~ValidationPlots()
{
  delete elasticHists;
  delete hitDistributions;
}

//----------------------------------------------------------------------------------------------------

void ValidationPlots::beginRun(edm::Run const&, edm::EventSetup const& es)
{
  TH1::AddDirectory(kFALSE);
  elasticHists->initialize(es);
}

//----------------------------------------------------------------------------------------------------

void ValidationPlots::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  elasticHists->analyze(event,eSetup);
  hitDistributions->analyze(event,eSetup);
}

//----------------------------------------------------------------------------------------------------

void ValidationPlots::endJob()
{
  elasticHists->finalize();
  elasticHists->ExportAllHistograms();
  hitDistributions->ExportAllHistograms();
  std::vector<TCanvas*>  elHists = elasticHists->getHistograms();
  std::vector<TCanvas*>  histGraphs = hitDistributions->getHistograms();


  TFile *f = TFile::Open(outputFile.c_str(), "recreate");
  if(!f || !f->IsWritable())
    {
      std::cout<<"Output file not opened correctly!!"<<std::endl;
    }
  if(verbosity)
    std::cout<<"Writting histograms..."<<std::endl;

  TH1::AddDirectory(kFALSE);
  gDirectory->cd("/");
  if(!gDirectory->cd("elastic"))
    {
      gDirectory->mkdir("elastic");
      gDirectory->cd("elastic");
    }

  for(unsigned int i=0; i<elHists.size(); ++i)
    {
      elHists[i]->Write("");
    }

  gDirectory->cd("/");
  if(!gDirectory->cd("hits"))
    {
      gDirectory->mkdir("hits");
      gDirectory->cd("hits");
    }

  for(unsigned int i=0; i<histGraphs.size(); ++i)
    {
      histGraphs[i]->Write("");
    }

  if(verbosity)
    std::cout<<"Writting histograms finnished."<<std::endl;
  f->Close();

}

DEFINE_FWK_MODULE(ValidationPlots);

#ifndef _TotemRPValidationHitDistributionsHitDistributionsLibrary_H_
#define _TotemRPValidationHitDistributionsHitDistributionsLibrary_H_

/*
 * HitDistributionsLibrary.h
 *
 *  Created on: Oct 29, 2008
 *      Author: Leszek Grzanka
 */

#include <stdio.h>
#include <iostream>
#include <map>

#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <TH2F.h>
#include <TH1.h>
#include <TStyle.h>
#include <TFile.h>
#include <TKey.h>
#include <TLine.h>
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TSystem.h"
#include "TRint.h"
#include "TStyle.h"
#include "TMath.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPTrackCandidateCollection.h"
#include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDetTrigger.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "DataFormats/TotemRPDataTypes/interface/RPRecoHit.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecoElasticEvent.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "HepMC/GenEvent.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "Geometry/TotemRecords/interface/RealGeometryRecord.h"
#include "Geometry/TotemRPGeometryBuilder/interface/DetGeomDesc.h"
#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"


namespace edm {
  class ParameterSet;
  class EventSetup;
  class Event;
}

class HitDistributionsLibrary
{
 public:
  HitDistributionsLibrary(const edm::ParameterSet&);
  ~HitDistributionsLibrary();

  void analyze(const edm::Event&, const edm::EventSetup&);
  void writeToFile();

  std::vector<TCanvas*> getHistograms();
  void ExportAllHistograms();

 private:
  unsigned int verbosity;
  unsigned int validTracksOnly;
  std::vector<TCanvas*> hitsGraphs;

  std::string outputFile;
  std::map< unsigned int, TGraph *> hitDists;
  std::set< unsigned int> unitIDs;
  std::map< unsigned int, CLHEP::Hep3Vector> rpTranslations;
  std::map< unsigned int, CLHEP::HepRotation> rpRotations;
  edm::InputTag rpFittedTrackCollectionLabel;

  void drawRPPlaneEnvelopes(unsigned int rpId);
  void setupCanvas();

};

#endif /* _TotemRPValidationHitDistributionsHitDistributionsLibrary_H_ */

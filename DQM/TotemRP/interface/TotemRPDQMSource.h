#ifndef TotemRPDQMSource_H
#define TotemRPDQMSource_H

// Framework
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// event
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

// DQM
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

// Electron
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

// PFMET
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

// Vertex utilities
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// CaloJets
#include "DataFormats/JetReco/interface/CaloJet.h"

// Conversions
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

// Trigger
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDigCluster.h"
#include "DataFormats/TotemRPDataTypes/interface/RPRecoHit.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecognizedPatternsCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrack.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPTrackCandidateCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPMulFittedTrackCollection.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "Geometry/TotemRPDetTopology/interface/RPTopology.h"


 
class TotemRPDQMSource: public DQMEDAnalyzer
{
	public:
		TotemRPDQMSource(const edm::ParameterSet& ps);
		virtual ~TotemRPDQMSource();
  
	protected:
		void dqmBeginRun(edm::Run const &, edm::EventSetup const &) override;
		void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
		void analyze(edm::Event const& e, edm::EventSetup const& eSetup);
		void beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& eSetup);
		void endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& eSetup);
		void endRun(edm::Run const& run, edm::EventSetup const& eSetup);

	private:
  	  edm::EDGetTokenT< edm::DetSetVector<RPStripDigi> > tokenStripDigi;
  	  edm::EDGetTokenT< edm::DetSetVector<RPDigCluster> > tokenDigiCluster;
  	  edm::EDGetTokenT< edm::DetSetVector<RPRecoHit> > tokenRecoHit;
  	  edm::EDGetTokenT< RPRecognizedPatternsCollection > tokenPatternColl;
  	  edm::EDGetTokenT< RPTrackCandidateCollection > tokenTrackCandColl;
  	  edm::EDGetTokenT< RPFittedTrackCollection > tokenTrackColl;
  	  edm::EDGetTokenT< RPMulFittedTrackCollection > tokenMultiTrackColl;

      /*
      bool buildCorrelationPlots;                           ///< decides wheather the correlation plots are created
      std::string correlationPlotsFilter;                   ///< decides which correlation plots are created
      unsigned int correlationPlotsLimit;                   ///< maximum number of created correlation plots
      CorrelationPlotsSelector correlationPlotsSelector;
      */

		// Histograms
		MonitorElement* h_test;
};

#endif

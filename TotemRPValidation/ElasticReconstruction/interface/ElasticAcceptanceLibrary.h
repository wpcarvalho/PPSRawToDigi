/****************************************************************************
*
* This is a part of the TOTEM offline software.
* Authors: 
*	Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $Id: ElasticAcceptanceLibrary.h 9977 2015-01-12 14:00:26Z tsodzawi $
* $Revision: 9977 $
* $Date: 2015-01-12 16:00:26 +0200 (pon, 12 sty 2015) $
*
****************************************************************************/


#ifndef _ElasticAcceptanceLibrary_h_
#define _ElasticAcceptanceLibrary_h_

namespace edm {
  class ParameterSet;
  class EventSetup;
  class Event;
}

class TH1D;
class TCanvas;

class TrackData;

#include <string>
#include <map>


/**
 *\brief 
**/
class ElasticAcceptanceLibrary
{
 public:
  ElasticAcceptanceLibrary(const edm::ParameterSet&);
  ~ElasticAcceptanceLibrary();

  void initialize(const edm::EventSetup&);
  void analyze(const edm::Event&, const edm::EventSetup&);
  void finalize();
  void writeHistogramsToFile();

  struct AcceptanceCollection {
    TH1D *t_s;      ///< linear t, smeared
    TH1D *t_o;      ///< linear t, original
    TH1D *tx_s;     ///< linear tx, smeared
    TH1D *tx_o;     ///< linear tx, original
    TH1D *ty_s;     ///< linear ty, smeared
    TH1D *ty_o;     ///< linear ty, original
    TH1D *logt_s;   ///< logarithmic t, smeared
    TH1D *logt_o;   ///< logarithmic t, original
    AcceptanceCollection() : t_s(NULL), t_o(NULL), logt_s(NULL), logt_o(NULL) {}
    AcceptanceCollection(const char *prefix, bool reco = false);
    void Fill(const TrackData &tr_o, const TrackData &tr_s);
    void Normalize(const AcceptanceCollection&);
    void Write();
  };

 private:
  unsigned char verbosity;														///< verbosity level
  std::string generatorLabel;                                                   ///< label of the HepMC product to be loaded as smeared event
  std::string originalLabel;                                                    ///< label of the HepMC product to be loaded as event before smaering
  std::string outputFile;														///< file to store output plots

  AcceptanceCollection total_l, total_r, total_reco;
  std::map<unsigned int, AcceptanceCollection> accepted_rp;
  std::map<unsigned int, AcceptanceCollection> accepted_unit;
  std::map<unsigned int, AcceptanceCollection> accepted_station;
  std::map<unsigned int, AcceptanceCollection> accepted_arm;
  AcceptanceCollection accepted_subsystem;
  edm::InputTag rpFittedTrackCollectionLabel;
  edm::InputTag rpRecoElasticEventLabel;
  /* interface to Leszek's validation kit */
 public:
  void ExportAllHistograms();
  std::vector<TCanvas*> getHistograms()
    { return elValHists; }

 private:
  std::vector<TCanvas*> elValHists;

  void ExportHistogram(TH1D *h);
};

#endif 


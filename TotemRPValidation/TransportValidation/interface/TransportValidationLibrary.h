/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $$RCSfile: TransportValidationLibrary.h,v $: $
* $Revision: 1.2 $
* $Date: 2008/12/17 13:34:49 $
*
****************************************************************************/

#include <vector>
#include <map>
#include <string>
#include <boost/shared_ptr.hpp>
#include "FWCore/Utilities/interface/InputTag.h"

class TFile;
class TH1D;
class LHCOpticsApproximator;
namespace edm {
  class ParameterSet;
  class EventSetup;
  class Event;
}

/**
 *\brief Library for a validation of the proton transport functionality. Called by TransportValidation analyzer.
 * Compares hit positions obtained by three ways:\\
 * a) optics -> directly by optics transport, using the same LHCOpticsApproximator object as which is used for reconstruction
 * b) sim hits -> PSimHit records produced by G4 simulation
 * c) reco -> local, one RP track fit, i.e. RPFittedTrack
**/
class TransportValidationLibrary
{
 public:
  /// set of histograms for each RP
  struct RPInfo {
    /// histograms of difference reco - optics
    boost::shared_ptr<TH1D> reco_diff_x, reco_diff_y;

    /// histograms of difference sim hits - optics
    boost::shared_ptr<TH1D> simhit_diff_x, simhit_diff_y;
    RPInfo();
  };

  /// simhit information converted to the global frame
  struct SimHitInfo
  {
    unsigned char det;
    double x, y;
    SimHitInfo(unsigned char _d = 0, double _x = 0., double _y = 0.) : det(_d), x(_x), y(_y) {}
  };

  TransportValidationLibrary(const edm::ParameterSet&);
  ~TransportValidationLibrary();

 private:
  edm::InputTag rPFittedTrackCollectionLabel;

  friend class TransportValidation;

  unsigned int verbosity;
  double forwardThLimit;
  std::string outputFile;
  
  std::map<unsigned int, RPInfo> histograms;

  void initialize(const edm::EventSetup&);
  void analyze(const edm::Event&, const edm::EventSetup&);
  void finalize();
  void writeHistogramsToFile();
};


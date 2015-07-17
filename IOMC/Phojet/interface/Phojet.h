/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *  Leszek Grzanka (Leszek.Grzanka@cern.ch)
 *
 ****************************************************************************/

#ifndef PHOJET_H_
#define PHOJET_H_

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CLHEP/Random/RandomEngine.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include "IOMC/Phojet/interface/HepMCFileReader.h"

#include <iostream>
#include "TFile.h"
#include "TClass.h"

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/convenience.hpp>

using namespace HepMC;
using namespace std;
using namespace edm;
using namespace boost::filesystem;


class Phojet : public edm::EDProducer
{
public:
  Phojet(const edm::ParameterSet &);
  virtual ~Phojet();
  virtual  void beginJob();
  virtual  void endJob();

  static const int PID = 2212;

protected:
  unsigned int verbosity;               ///< verbosity level

  virtual void produce(edm::Event&, const edm::EventSetup&);
  CLHEP::HepRandomEngine *rndEng;       ///< random engine
  HepMCFileReader *reader_;

private:
  void run_phojet( const path & workDirectory, const path & phojetExeFile, const path & phojetCfgFile, const bool printCrossSection);
  int generate_input_card( const path & file, double cms_energy_GeV, long number_of_events, string phojet_process_description);

  int buffer_size;
  double cms_energy;
  string phojet_process_description;
  string phojet_executable;
  path  tmpDirectoryPath;
  int phojetStatus;
};



#endif /* PHOJET_H_ */

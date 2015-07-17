/****************************************************************************
 *
 * This module is directly copied from IOMC/FlatProtonLogKsiLogTGun and
 * expanded to support the generation of multiple proton events.
 * Original Authors:
 *   Hubert Niewiadomski (Hubert.Niewiadomski@cern.ch)
 * Secondary Authors:
 *   Zhang Zhengkui (zhang.zhengkui.fin@gmail.com)
 *
 * $Id: FlatProtonLogKsiLogTGunMul.h 0001 2010-06-28 10:00:00Z zhangz $
 * $Revision: 0001 $
 * $Date: 2010-06-28 10:00:00 +0100 (Mon, 28 Jun 2010) $
 *
 ****************************************************************************/


#ifndef IOMC_FlatProtonLogKsiLogTGunMul_FlatProtonLogKsiLogTGunMul_h
#define IOMC_FlatProtonLogKsiLogTGunMul_FlatProtonLogKsiLogTGunMul_h

#include <string>
#include "HepPDT/defs.h"
#include "HepPDT/TableBuilder.hh"
#include "HepPDT/ParticleDataTable.hh"

#include "HepMC/GenEvent.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"

#include <memory>
#include "boost/shared_ptr.hpp"

using namespace edm;
namespace HepMC {
	class GenEvent;
}


class FlatProtonLogKsiLogTGunMul : public edm::EDProducer
{
public:
	FlatProtonLogKsiLogTGunMul(const edm::ParameterSet&);
	virtual ~FlatProtonLogKsiLogTGunMul();

private:
    virtual void beginRun(Run&, EventSetup const&);
    virtual void produce(edm::Event&, const edm::EventSetup&);
	bool ComputeTheta(double ksi, double t, double &theta);  //[-1..0], [GeV], [rad]
	bool GenerateProton(double &px, double &py, double &pz, double &E);  //[GeV x 4]

protected:
	// data members
	bool right_arm_;						///< generate proton on the right arm
	bool left_arm_;							///< generate proton on the left arm
	double nominalEnergy_;
	double min_t_;
	double max_t_;
	double min_ksi_;
	double max_ksi_;
	double min_phi_;
	double max_phi_;
	int max_iter_number_;
	int proton_number_;						///< the number of protons per arm per event

	bool log_ksi_distribution_;
	bool log_t_distribution_;
	bool z_symetric_mode_;

	double log10_t_min_;
	double log10_t_max_;
	double log10_ksi_min_;
	double log10_ksi_max_;

	HepMC::GenEvent* fEvt_;
	edm::ESHandle<HepPDT::ParticleDataTable> fPDGTable;
	int verbosity_;
	CLHEP::HepRandomEngine* fRandomEngine;
	CLHEP::RandFlat* fRandomGenerator;

	double m0_;
	double p1_;
	int PartID_;
};


#endif

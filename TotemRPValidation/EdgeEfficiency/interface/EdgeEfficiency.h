/****************************************************************************
*
* This is a part of the TOTEM offline software.
* Authors:
*   Frigyes Janos Nemes (frigyes.janos.nemes@cern.ch) 
*
* $$RCSfile: EdgeEfficiency.cc,v $: $
* $Revision: 1 $
* $Date: 2010-08-20 10:45:42       (Fri, 20 Aug 2010) $
*
****************************************************************************/

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"

#include "TH1D.h"

/**
 *\brief This module creates histograms about the efficiency near the edge of the detectors. The efficiency calculated by dividing the so called accepted and expected histograms. For a given detector the
 * expected histogram is defined by connecting points (from RPFittedTrackCollection) its two reference detectors wihere a detector is a reference one if overlap with the tested detector. For a given detector the accepted  
 * histogram contains such expected hits where exist at least one "close" detected hit (from DetSetVector<RPRecoHit>).
 **/

using namespace std ;
using namespace edm ;


class EdgeEfficiency : public edm::EDAnalyzer  
{
	public:
	EdgeEfficiency(const edm::ParameterSet&);
	~EdgeEfficiency();

	private:
	edm::InputTag rPFittedTrackCollectionLabel;
	edm::InputTag detSetVectorRPRecoHitLabel;
	virtual void beginRun(edm::Run const&, edm::EventSetup const&);
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();

	void test(const edm::Event&, unsigned int) ;
	void testTwoRPTestBeam(const edm::Event&) ;
	void testPackages(const edm::Event&) ;
	void testGivenPairOfRPs(int, int, const edm::Event&) ;

	std::string RootFileName ;	
	unsigned int Verbosity ;
	unsigned int selectedTest ; 
	double tolerance_in_v ; 		// [mm]  This is the tolerance in the w coordinate.
	double tolerance_in_angle ;		// [rad] This is the tolerance in the angle.
	double RP_angle_resolution ;
 	double RP_size_along_the_z_axis ;
	double resolution_in_w ;                // [mm]
	double stripZeroPosition ;

	edm::ESHandle<TotemRPGeometry> geom;	// Geometrial data of the RPs

	map<pair<unsigned int, unsigned int>, map<string, TH1D *> > histograms; // The key of the map is a pair which contain an RPId and a detector Id (from the given RP). The second will contain three histogram for the given detector.
	map< unsigned int, pair<unsigned int, unsigned int> > map_from_tested_RP_to_its_reference_RPs ; // The key of the map is the RPID of the tested detector. The value is a pair which contain the two reference detectors.
};



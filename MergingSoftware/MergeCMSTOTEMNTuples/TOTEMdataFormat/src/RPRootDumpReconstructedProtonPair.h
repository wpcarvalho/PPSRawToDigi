#ifndef _RPRootDumpReconstructedProtonPair_h_
#define _RPRootDumpReconstructedProtonPair_h_

#include <TObject.h>

/**
 * The reconstructed proton pair is stored in TotemNtuple as a separate branch. However, in common CMS-Totem
 * ntuple it will be stored as a part of RPEvent. For labels please look at RPEvent.h
 */
class RPRootDumpReconstructedProtonPair: public TObject {
public:
	RPRootDumpReconstructedProtonPair();

	bool valid;
	double thxr, thyr, xir, phir;
	double thxl, thyl, xil, phil;
	double x0, y0, z0, chi2, chindf;
	double tr, txr, tyr;
	double tl, txl, tyl;
	double t;

	ClassDef(RPRootDumpReconstructedProtonPair,2);
};

#endif

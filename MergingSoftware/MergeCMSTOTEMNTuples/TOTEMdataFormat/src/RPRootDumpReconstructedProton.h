#ifndef _RPRootDumpReconstructedProton_h_
#define _RPRootDumpReconstructedProton_h_

#include <TObject.h>

/**
 * The reconstructed proton is stored in TotemNtuple as a separate branch. However, in common CMS-Totem
 * ntuple it will be stored as a part of RPEvent. For labels please look at RPEvent.h
 */
class RPRootDumpReconstructedProton: public TObject {
public:
	RPRootDumpReconstructedProton();

	bool valid;
	double thx, thy;
	double phi;
	double t, tx, ty;
	double xi, x0, y0;
	double chi2, chindf;

	ClassDef(RPRootDumpReconstructedProton,2);
};

#endif

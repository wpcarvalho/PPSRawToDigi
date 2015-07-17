#ifndef _RPEvent_h_
#define _RPEvent_h_

#include <TObject.h>
#include "RPRootDumpReconstructedProton.h"
#include "RPRootDumpReconstructedProtonPair.h"

/**
 * This class contains information about reconstructed and simulated protons. In totem_ntuple there are
 * separate branches for that informations in main directory of ntuple. In common CMS-Totem ntuple we
 * decided to merge this data into one class.
 * In totem ntuple the branches for reconstructed protons are named:
 *   rec_prot_left - "rec_prot_left."
 *   rec_prot_right - "rec_prot_right."
 *   rec_prot_pair - "rec_prot_pair."
 */
class RPEvent: public TObject {
public:
	RPEvent();

	RPRootDumpReconstructedProton rec_prot_left, rec_prot_right;
	RPRootDumpReconstructedProtonPair rec_prot_pair;

	ClassDef(RPEvent, 1);
};

#endif

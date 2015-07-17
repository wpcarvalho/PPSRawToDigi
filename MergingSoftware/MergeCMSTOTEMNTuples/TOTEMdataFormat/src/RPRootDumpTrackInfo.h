#ifndef _RPRootDumpTrackInfo_h_
#define _RPRootDumpTrackInfo_h_

#include <vector>
#include <TObject.h>

/**
 * This class contains tracks information. For tracks, there is a separate branch for every pot
 * in the main directory of TotemNtuple.
 * The convention for naming branches is: track_rp_$number.
 * Where number is the number of pot: 0,1,2,3,4,5 20,21,22,23,24,25 100,101,102,103,104,105 120,121,122,123,124,125
 *
 * Morevover, in the main directory of TotemNtuple there are optionally branches called multi_track_rp_$number.
 * These branches contain vector<RPTrackInfo>
 */
class RPRootDumpTrackInfo: public TObject {
public:
	RPRootDumpTrackInfo();

	bool valid; // whether track fit is valid
	double x, y, z; // track fit interpolated to the middle of the RP
	double chi2; // fit chi square
	double chi2ndf; // fit chi square divided by the number of degrees of freedom
	unsigned int entries; // the number of contributing hits
	double res_x, res_y; // seem not used
	std::vector<int> u_sect, v_sect; // list of active trigger sectors calculated from (strip) data
	int u_sect_no, v_sect_no; // sizes of u_sect and v_sect vectors

	ClassDef(RPRootDumpTrackInfo,2);
};

#endif

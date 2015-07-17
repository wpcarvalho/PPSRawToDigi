#include "RPRootDumpTrackInfo.h"

RPRootDumpTrackInfo::RPRootDumpTrackInfo() {
	valid = false;
	chi2 = chi2ndf = x = y = z = res_x = res_y = 0.0;
	entries = 0;
	u_sect = std::vector<int>();
	v_sect = std::vector<int>();
	u_sect_no = v_sect_no = 0;
}

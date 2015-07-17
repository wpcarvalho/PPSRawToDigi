#include "TProfile3D.h"
#include "TMath.h"
#include "TotemRPValidation/ValidationTools/interface/ReconstructionProfile.h"

#include "vector"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TH2F.h"
#include "TH1.h"
#include "TStyle.h"
#include "TGraph.h"

class Inelastic1ArmPlots{

	public:

	explicit Inelastic1ArmPlots();
	~Inelastic1ArmPlots();
	std::vector<TCanvas*> getPlots();

	private:

	double tmin;
	double tmax;
	double interpoint_t_dist;
	double ksimin;
	double ksimax;
	double interpoint_ksi_dist;
	double event_cut;

	void SingleArmEntries(ReconstructionProfile *p, TPad *pad);
	void SingleArmSigma(ReconstructionProfile *p, const char *y_label, TPad *pad, bool logscale=false);
	void SingleArmMean(ReconstructionProfile *p, char *y_label, TPad *pad);
};

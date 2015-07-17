#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TSystem.h"
#include "TRint.h"
#include "TStyle.h"

#include "TotemRPValidation/ValidationPlots/interface/InelasticPlots.h"
#include "TotemRPValidation/ValidationPlots/interface/Inelastic1ArmPlots.h"

InelasticPlots::InelasticPlots(){
}

InelasticPlots::~InelasticPlots(){
}

void InelasticPlots::prepareResPlots( std::string filename ){
	TFile * fi = new TFile(filename.c_str());
	fi->cd("prot_rec_val/Arm_0001");
	Inelastic1ArmPlots * iep = new Inelastic1ArmPlots();
	resolutionPlots = iep->getPlots();
	fi->Close();
}

std::vector<TCanvas*> InelasticPlots::getResolutionPlots(){
	return resolutionPlots;
}


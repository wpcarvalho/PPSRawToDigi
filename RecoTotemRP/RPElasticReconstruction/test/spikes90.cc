#include "TRandom2.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"

//#define M_PI 3.141593

void RoundToStrip(double &y)
{
	double res = 66E-6 / sqrt(2.);
	if (gRandom->Rndm() < 1.6) y = floor(y / res + 0.5) * res;
}

void spikes90()
{
	gStyle->SetOptStat(111110);
	gStyle->SetCanvasColor(0);	
	gStyle->SetCanvasBorderMode(0);	
	gStyle->SetFrameBorderMode(0);	
	gStyle->SetTitleBorderSize(1);
	gStyle->SetStatBorderSize(1);	
	gStyle->SetFuncWidth(1);
	gStyle->SetFuncColor(4);
	gStyle->SetStatW(0.25);
	gStyle->SetStatH(0.25);

	TH1D *h = new TH1D("h", "h", 5000, -1E-5, 1E-5);
	TH1D *hits = new TH1D("hits", "hits", 200, -1E-4, 1E-4);

	double si_de = 2.4E-6;
	double B = 20.;
	double P = 7E3;

	//double L1 = 239., L2 = 264.6;
	double L1 = 264.6, L2 = 264.6;

	gRandom->SetSeed(0);

	for (int i = 0; i < 4000; i++) {
		// physics
		double p = gRandom->Rndm();
		double t = - 1./B * log(1. - p);
		double th = sqrt(t) / P;
		double phi = gRandom->Rndm() * 2 * M_PI; 
		//double th_x = th * cos(phi);
		double th_y = th * sin(phi);

		// beam divergence
		double th1 = fabs(gRandom->Gaus()) * si_de;
		double th2 = fabs(gRandom->Gaus()) * si_de;
		double ph1 = gRandom->Rndm() * 2. * M_PI; 
		double ph2 = gRandom->Rndm() * 2. * M_PI; 

		// the two protons
		double thr_y = th_y + th1 * sin(ph1);
		double thl_y = th_y + th2 * sin(ph2);

		// the four hits
		double y1r = L1 * thr_y;
		double y2r = L2 * thr_y;
		double y1l = -L1 * thl_y;
		double y2l = -L2 * thl_y;

		// rounding
		RoundToStrip(y1r);
		RoundToStrip(y2r);
		RoundToStrip(y1l);
		RoundToStrip(y2l);

		hits->Fill(y1r);

		// reconstruction
		double thpr_y = (y1r*L1 + y2r*L2) / (L1*L1 + L2*L2);
		double thpl_y = (- y1l*L1 - y2l*L2) / (L1*L1 + L2*L2);

		// histograms
		h->Fill(thpr_y - thpl_y);
	}

	TCanvas *can = new TCanvas("can", "peak sim", 500, 500);
	can->SetCrosshair(1);
	can->ToggleEventStatus();
	h->Draw("");

	//can = new TCanvas();
	//hits->Draw("");
}

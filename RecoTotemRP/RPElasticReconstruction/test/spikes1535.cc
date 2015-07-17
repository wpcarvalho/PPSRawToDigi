#include "TRandom2.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"

//#define M_PI 3.141593

double shiftSigma = 0E-6;
double shifts[12];


void RoundToStrip(double &y, signed int idx)
{
	double res = 66E-6 / sqrt(2.);
	if (gRandom->Rndm() < 1.7) y = floor(y / res + 0.5) * res;

	int index;
	if (idx < 0) index = 6; else index = 0;
	if (abs(idx) == 1 && (y < 0)) index += 0;
	if (abs(idx) == 1 && (y > 0)) index += 1;
	if (abs(idx) == 2) index += 2;
	if (abs(idx) == 3) index += 3;
	if (abs(idx) == 4 && (y < 0)) index += 4;
	if (abs(idx) == 4 && (y > 0)) index += 5;

	//printf("idx = %i => index = %i, shift = %.3E\n", idx, index, shifts[index]);

   	y += shifts[index];
}

void spikes1535()
{
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetCanvasColor(0);	
	gStyle->SetCanvasBorderMode(0);	
	gStyle->SetFrameBorderMode(0);	
	gStyle->SetTitleBorderSize(1);
	gStyle->SetStatBorderSize(1);	
	gStyle->SetFuncWidth(1);
	gStyle->SetFuncColor(4);
	gStyle->SetStatW(0.25);
	gStyle->SetStatH(0.25);
	gStyle->SetTitleYOffset(1.5);
	gStyle->SetTitleXOffset(1.2);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.08);
	

	TH1D *h = new TH1D("h", "h;#Delta_{R-L}#vartheta_{y}   (rad)", 50000, -1E-6, 1E-6);
	TH1D *hits = new TH1D("hits", "hits", 200, -1E-4, 1E-4);

	double si_de = 0.295E-6;
	double B = 20.;
	double P = 7E3;

	double L1, L2, L3, L4; L1 = L2 = L3 = L4 = 270.; // y
	//double L1 = 248.8, L2 = 250.8, L3 = 270.7, L4 = 272.7; // y
	//double L1 = 112.5, L2 = 111.4, L3 = 99.7, L4 = 98.6; // x

	gRandom->SetSeed(2);

	// initialize shifts
	for (int i = 0; i < 12; i++) {
		shifts[i] = gRandom->Gaus() * shiftSigma;
		printf("shifts[%i] = %.3E\n", i, shifts[i]);
	}

	for (int i = 0; i < 10000; i++) {
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
		double y3r = L3 * thr_y;
		double y4r = L4 * thr_y;
		double y1l = -L1 * thl_y;
		double y4l = -L4 * thl_y;

		// decisions
		bool config42 = false; //(gRandom->Rndm() > 0.4);

		// rounding
		if (true) {
			RoundToStrip(y1r, 1);
			RoundToStrip(y2r, 2);
			RoundToStrip(y3r, 3);
			RoundToStrip(y4r, 4);
			RoundToStrip(y1l, -1);
			RoundToStrip(y4l, -4);
		}

		//printf("hits: %.3E\t%.3E\t%.3E\t%.3E\n", y4l, y1l, y1r, y4r);

		hits->Fill(y1r);

		// reconstruction
		double thpr_y;
		if (config42) thpr_y = (y1r*L1 + y2r*L2 + y3r*L3 + y4r*L4) / (L1*L1 + L2*L2 + L3*L3 + L4*L4);
		else thpr_y = (y1r*L1 + y4r*L4) / (L1*L1 + L4*L4);
		double thpl_y = (- y1l*L1 - y4l*L4) / (L1*L1 + L4*L4);

		// histograms
		//printf("dth = %.3E\n\n", thpr_y - thpl_y);
		h->Fill(thpr_y - thpl_y);
	}

	TCanvas *can = new TCanvas("can", "peak sim", 300, 300);
	can->SetCrosshair(1);
	can->ToggleEventStatus();
	h->Draw("");

	//can = new TCanvas();
	//hits->Draw("");
}

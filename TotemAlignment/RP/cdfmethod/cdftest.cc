#include <stdio.h>

#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFile.h"
#include "TSystem.h"
#include "TRint.h"
#include "TRandom2.h"
#include "TMath.h"
#include "TStyle.h"

#include <vector>

//----------------------------------------------------------------------------------------------------

using namespace std;


vector<TGraph *> cdf;

char LoadCDF(char *fileName, double cdfCropMin = 1E-3, double cdfCropMax = 10.5)
{
	TFile cdfFile(fileName);
	if (cdfFile.IsZombie()) return 2;

	const int maxDNum = 6;
	char *tags[maxDNum] = {"islam", "ppp2", "ppp3", "bsw", "bh", "exp"};

	Double_t FullIntegral[maxDNum];
	int cdfPoints = 0;
	double cdfMin, cdfMax;

	for (int i = 0; i < maxDNum; i++) {
		/// try to load the graph
		TGraph *g = (TGraph *) cdfFile.Get(tags[i]);
		if (g) cdf.push_back(g);
		else continue;

		/// crop the cdf
		double x, y;
		g->GetPoint(0, x, y);
		while (x < cdfCropMin) { g->RemovePoint(0); g->GetPoint(0, x, y); }
		g->GetPoint(g->GetN() - 1, x, y);
		while (x > cdfCropMax) { g->RemovePoint(g->GetN() - 1); g->GetPoint(g->GetN() - 1, x, y); }

		/// normalize and invert cdf
		g->GetPoint(g->GetN() - 1, x, y);
		double intMin;
		g->GetPoint(0, x, intMin);
		FullIntegral[i] = y - intMin;
		cdfMin = x;
		for (int j = 0; j < g->GetN(); j++) {
			g->GetPoint(j, x, y);
			g->SetPoint(j, (y - intMin) / FullIntegral[i], x);
			cdfPoints++;
		}
		cdfMax = x;
	}

	return 0;
}


//----------------------------------------------------------------------------------------------------



int main(int argc, char *argv[])
{
	/// create fake application
	int fakeC = 2; char *fakeV[2] = {"prdlajs", "-l"};
	TRint app("application", &fakeC, fakeV, 0, true);
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
	gStyle->SetLegendBorderSize(1);

	// init random number generator
	gRandom->SetSeed(0);

	// color palette
	int colors[] = {1, 2, 4, 6, 8, 5, 7};
	
	// parameters
	double b = 20.;			// GeV^-2
	double t_min =   0.;	// GeV^2
	double t_max = 100.;	// ...
	int Ndet = 4;
	double Ly = 250.;		// m
	double Lx = 100.;		// m
	unsigned long N = (unsigned long)(5E7);
	double p = 7E3;			// GeV
	vector<double> sigmas;	// m
	int cdfIdx = 5;
	int bins = 500;

	// Load cdf file
	if (LoadCDF("../generator/sigma_int/sigma_int,KL,1E-3,1.5E0.root", t_min, t_max)) return 1;

	// list of sigmas
	sigmas.push_back(0.);

	sigmas.push_back(1E-6);
	sigmas.push_back(2E-6);
	sigmas.push_back(5E-6);
	sigmas.push_back(1E-5);
	sigmas.push_back(1.2E-5);
	sigmas.push_back(2E-5);
	sigmas.push_back(3E-5);
	sigmas.push_back(5E-5);
	sigmas.push_back(1E-4);
	sigmas.push_back(2E-4);
	sigmas.push_back(5E-4);
	sigmas.push_back(1E-3);

	sigmas.push_back(2E-3);
	sigmas.push_back(4E-3);
	sigmas.push_back(6E-3);
	sigmas.push_back(8E-3);
	sigmas.push_back(1E-2);

	sort(sigmas.begin(), sigmas.end());	

	// initilize histograms
	// generate angular errors
	vector<double> Dth_x, Dth_y;
	vector<TH1D *> tHist;
	for (unsigned int i = 0; i < sigmas.size(); i++) {
		Dth_x.push_back(sigmas[i]/sqrt(Ndet)/Lx);
		Dth_y.push_back(sigmas[i]/sqrt(Ndet)/Ly);
		//printf("sigma(De) = %.2E, sigma(th_x) = %.2E, Delta(th_x) = %.2E\n", sigmas[i], sigmas[i]/sqrt(Ndet)/Lx, Dth_x[i]);
		//printf("sigma(De) = %.2E, sigma(th_y) = %.2E, Delta(th_y) = %.2E\n\n", sigmas[i], sigmas[i]/sqrt(Ndet)/Ly, Dth_y[i]);

		char name[50], title[50];
		sprintf(name, "tHist%i", i);
		sprintf(title, "#sigma=%.2E m;|t|   (GeV^{2})", sigmas[i]);
		tHist.push_back(new TH1D(name, title, bins, 0., 1.5));
	}

	// simulation
	for (unsigned long i = 0; i < N; i++) {
		if (i % 1000000 == 0) printf("%i\n", i);
		double P = gRandom->Rndm();
		/*
		double I = exp(-b*t_min), M = exp(-b*t_max);
		double t = - log(I - P * (I - M)) / b;
		*/
		double t = cdf[cdfIdx]->Eval(P);
		double phi = gRandom->Rndm() * 2.*M_PI;
		double th_x = sqrt(t) / p * cos(phi);
		double th_y = sqrt(t) / p * sin(phi);

		for (unsigned int i = 0; i < sigmas.size(); i++) {
			double thp_x = th_x + Dth_x[i];
			double thp_y = th_y + Dth_y[i];
			double tp = p*p * (thp_x*thp_x + thp_y*thp_y);
			//double tp = sqrt(thp_x*thp_x + thp_y*thp_y);
			//double tp = fabs(thp_x);
			tHist[i]->Fill(tp);
		}

	}

	// fitting
	//vector<double> slopes;
	TF1 *ff = new TF1("ff", "[0] * exp(-[1]*x)");
	TGraphErrors *gSS = new TGraphErrors(); gSS->SetTitle(";log(#Delta / m);slope B   (GeV^{-2})"); gSS->SetMarkerStyle(2); gSS->SetMarkerSize(1);
	for (int i = 0; i < sigmas.size(); i++) {
		tHist[i]->Fit(ff, "QN", "", 0.002, 0.3);
		double slope = ff->GetParameter(1);
		double slopeErr = ff->GetParError(1);
		if (sigmas[i] > 0.) {
			gSS->SetPoint(gSS->GetN(), log10(sigmas[i]), slope);
			gSS->SetPointError(gSS->GetN() - 1, 0., slopeErr);
		}
		printf("%.2E\t%.2E\n", sigmas[i], slope);
	}

	// draw results
	TCanvas *c;
	c = new TCanvas("c", "", 400, 400);
	//c->Divide(3, 1);
	c->SetLogy(1);
	tHist[0]->Draw(); tHist[0]->SetLineColor(1);
	tHist[13]->Draw("same"); tHist[13]->SetLineColor(2);
	tHist[14]->Draw("same"); tHist[14]->SetLineColor(4);
	tHist[15]->Draw("same"); tHist[15]->SetLineColor(5);
	tHist[16]->Draw("same"); tHist[16]->SetLineColor(6);
	tHist[17]->Draw("same"); tHist[17]->SetLineColor(8);
	/*
	for (int i = 0; i < sigmas.size(); i++) {
		//tHist[i]->Draw((i) ? "same" : "");
		//tHist.back()->SetLineColor(colors[i]);
	}
	*/
	c->BuildLegend();


	c = new TCanvas("c2", "", 400, 400);
	//c->SetLogx(1);
	gSS->Draw("APL");

	/// run the fake application
	app.Run();
	return 0;
}

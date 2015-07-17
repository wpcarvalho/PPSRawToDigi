#include "TH2D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TRint.h"
#include "TStyle.h"
#include "TColor.h"

#include <math.h>
#include <vector>
#include <map>

using namespace std;

//----------------------------------------------------------------------------------------------------

struct optFun {
	double Lx, Ly, vx, vy;
	optFun(double _lx = 0., double _ly = 0., double _vx = 0., double _vy = 0.) : Lx(_lx), Ly(_ly), vx(_vx), vy(_vy) {}
};

map<unsigned short, optFun> optics;

void InitOptics()
{
	optics[0] = optFun(0., 0., 1., 1.);
	optics[20] = optFun(-2.867, -239., -2.177, 0.02153);
	optics[21] = optFun(-2.867, -239., -2.177, 0.02153);
	optics[22] = optFun(-2.627, -241.1, -2.151, 0.01998);
	optics[23] = optFun(-0.2413, -262.5, -1.9, 0.003149);
	optics[24] = optFun(-0.0012, -264.6, -1.875, 0.00142);
	optics[25] = optFun(-0.0012, -264.6, -1.875, 0.00142);
	optics[120] = optFun(2.867, 239., -2.177, 0.02153);
	optics[121] = optFun(2.867, 239., -2.177, 0.02153);
	optics[122] = optFun(2.627, 241.1, -2.151, 0.01998);
	optics[123] = optFun(0.2413, 262.5, -1.9, 0.003149);
	optics[124] = optFun(0.0012, 264.6, -1.875, 0.00142);
	optics[125] = optFun(0.0012, 264.6, -1.875, 0.00142);
}

//----------------------------------------------------------------------------------------------------


struct hit {
	unsigned short RPId;
	double x, y;
	/// converts from mm to m
	hit(unsigned short _id = 0, double _x = 0., double _y = 0.) : RPId(_id), x(_x * 1E-3), y(_y * 1E-3) {}
};

vector<hit> hits;


void FillHits(int event = 0)
{
	hits.clear();

	if (event == 0) {
		hits.push_back(hit(20, -0.140, 11.094));
		hits.push_back(hit(24, -8.010, 17.367));
		hits.push_back(hit(121, 0.152, -11.106));
		hits.push_back(hit(125, 0.011, -12.286));
	}

	if (event == 1) {
		hits.push_back(hit(20, -0.043, 7.924));
		hits.push_back(hit(24, -0.013, 8.784));
		hits.push_back(hit(121, 0.051, -7.917));
		hits.push_back(hit(122, 5.997, -7.733));
		hits.push_back(hit(123, 13.492, -12.865));
	}
	
	if (event == 2) {
		hits.push_back(hit(20, -0.007, 9.654));
		hits.push_back(hit(24, -0.000, 10.664));
		hits.push_back(hit(121, 0.000, -9.647));
		hits.push_back(hit(123, 19.645, -20.646));
	}
	
	if (event == 3) {
		hits.push_back(hit(20, 0.047, 10.161));
		hits.push_back(hit(24, -0.023, 11.255));
		hits.push_back(hit(121, -0.047, -10.161));
		hits.push_back(hit(122, 2.222, -5.699));
	}

	hits.push_back(hit(0, 0., 0.));
}

//----------------------------------------------------------------------------------------------------


void HistFillPoint(TH2D *h, unsigned int binx, unsigned int biny, double weight, bool swap)
{
	if (!swap) {
		//.printf("\t\t %i, %i\n", binx, biny);
		h->SetBinContent(binx, biny, h->GetBinContent(binx, biny) + weight);
	} else {
		//.printf("\t\t %i, %i\n", biny, binx);
		h->SetBinContent(biny, binx, h->GetBinContent(biny, binx) + weight);
	}
}

inline double fpart(double x)
{
	return x - floor(x);
}


inline double rfpart(double x)
{
	return 1. - fpart(x);
}


void HistFillLine(TH2D *h, double x1, double x2, double y1, double y2)
{
	printf("\tFilling line: x from %f to %f, y from %f to %f\n", x1, x2, y1, y2);

	double swapped = false;
	if (y2 - y1 > x2 - x1) {
		swapped = true;
		swap(x1, y1);
		swap(x2, y2);
	}

    double dx = x2 - x1;
    double dy = y2 - y1;
    double gradient = dy / dx;

    // handle first endpoint
    double xend = floor(0.5 + x1);
    double yend = double(y1) + gradient * (xend - x1);
    double xgap = rfpart(0.5 + x1);
    double xpxl1 = xend;
    double ypxl1 = floor(yend);
    HistFillPoint(h, xpxl1, ypxl1, rfpart(yend) * xgap, swapped);
    HistFillPoint(h, xpxl1, ypxl1 + 1, fpart(yend) * xgap, swapped);

    double intery = yend + gradient; // first y-intersection for the main loop

    // handle second endpoint
    xend = floor(0.5 + x2);
    yend = y2 + gradient * (xend - x2);
    xgap = fpart(0.5 + x2);
    double xpxl2 = xend;  
    double ypxl2 = floor(yend);
    HistFillPoint(h, xpxl2, ypxl2, rfpart(yend) * xgap, swapped);
    HistFillPoint(h, xpxl2, ypxl2 + 1, fpart(yend) * xgap, swapped);

    // main loop
    for (double x = xpxl1 + 1; x <= xpxl2 - 1; x++) {
            HistFillPoint(h, x, floor(intery), rfpart(intery), swapped);
            HistFillPoint(h, x, floor(intery) + 1, fpart(intery), swapped);
            intery = intery + gradient;
	}


	/*.
	double slope = (y2 - y1)/(x2 - x1);
	for (unsigned int x = x1; x <= x2; x++) {
		HistFillPoint(h, x, int(slope*(x - x1) + y1), 1., swapped);
	}
	*/
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// create fake application
	int fakeC = 2; char * fakeV[2] = {"prdlajs", "-l"};
	TRint app("application", &fakeC, fakeV, 0, true);
	gStyle->SetOptStat(1110);
	gStyle->SetCanvasColor(0);	
	gStyle->SetCanvasBorderMode(0);	
	gStyle->SetFrameBorderMode(0);	
	gStyle->SetTitleBorderSize(1);
	gStyle->SetStatBorderSize(1);	
	gStyle->SetFuncWidth(1);
	gStyle->SetFuncColor(4);
	gStyle->SetMarkerColor(2);
	gStyle->SetMarkerStyle(2);
	gStyle->SetMarkerSize(2.);
	double cLength[] = {0., 0.33, 0.66, 1.};
	double cRed[] = {1., 1., 0., 0.};
	double cGreen[] = {.9, 0., 0., 1.};
	double cBlue[] = {.9, 0., 1., 0.1};
	#if ROOT_VERSION_CODE >= ROOT_VERSION(5,16,0)
	TColor::CreateGradientColorTable(4, cLength, cRed, cGreen, cBlue, 12);
	#endif

	// init
	InitOptics();
	FillHits(1);

	// draw graphs
	TGraph *hitXvsLx = new TGraph(); hitXvsLx->SetTitle(";L_{x}   (m);x   (m)");
	TGraph *hitYvsLy = new TGraph(); hitYvsLy->SetTitle(";L_{y}   (m);y   (m)");
	TGraph *hitXvsvx = new TGraph(); hitXvsvx->SetTitle(";v_{x}   (m);x   (m)");
	TGraph *hitYvsvy = new TGraph(); hitYvsvy->SetTitle(";v_{y}   (m);y   (m)");

	for (unsigned int i = 0; i < hits.size(); i++) {
		hitXvsLx->SetPoint(hitXvsLx->GetN(), optics[hits[i].RPId].Lx, hits[i].x);
		hitYvsLy->SetPoint(hitYvsLy->GetN(), optics[hits[i].RPId].Ly, hits[i].y);

		hitXvsvx->SetPoint(hitXvsvx->GetN(), optics[hits[i].RPId].vx, hits[i].x);
		hitYvsvy->SetPoint(hitYvsvy->GetN(), optics[hits[i].RPId].vy, hits[i].y);
	}

	TCanvas *c = new TCanvas("c", "hits", 600, 600);
	c->Divide(2, 2);
	c->cd(1); hitXvsLx->Draw("AP");
	c->cd(2); hitYvsLy->Draw("AP");
	c->cd(3); hitXvsvx->Draw("AP");
	c->cd(4); hitYvsvy->Draw("AP");

	// do Hough transform
	int Nbins = 150;
	TH2D *accu = new TH2D("accu", ";#vartheta_{y};v^{*}_{y}", Nbins, -1E-4, 1E-4, Nbins, -1E1, 1E1); // theta [rad], y* [um]
	c = new TCanvas("c2", "bla", 1000, 200);
	c->Divide(5, 1);
	int count = 0;
	for (unsigned int i = 0; i < hits.size(); i++) {
		double &Ly = optics[hits[i].RPId].Ly;
		double &vy = optics[hits[i].RPId].vy;
		double &y = hits[i].y;
		printf("RPId = %i, y = %.2E, Ly = %.2E, vy = %.2E, th_y_guess = %.2E\n", hits[i].RPId, y, Ly, vy, y / Ly);

		// determine end-points
		unsigned int thBeg = 0, thEnd = 0, ystBeg = 0, ystEnd = 0;
		for (int bth = 1; bth <= accu->GetNbinsX(); bth++) {
			double th = accu->GetXaxis()->GetBinCenter(bth);
			double yst = (y - Ly * th) / vy;
			int byst = accu->GetYaxis()->FindBin(yst);

			if (byst > 0 && byst <= accu->GetNbinsY()) {
				if (thBeg == 0) { thBeg = bth; ystBeg = byst; }
				thEnd = bth; ystEnd = byst;		
				//.printf("\tbin(th) = %i, th = %.2E, y* = %.2E, bin(y*) = %i\n", bth, th, yst, byst);
			}
		}

		// fill line
		if (thBeg > 0) HistFillLine(accu, thBeg, thEnd, ystBeg, ystEnd);

		// draw progress
		TH2D *h = new TH2D(*accu);
		c->cd(++count);
		h->Draw("colz");	
	}

	new TCanvas();
	accu->Draw("colz");

	double max = -1., thMax, ystMax;
	// find maximum
	for (int i = 1; i <= accu->GetNbinsX(); i++)
		for (int j = 1; j <= accu->GetNbinsY(); j++)
			if (accu->GetBinContent(i, j) > max) {
				max = accu->GetBinContent(i, j);
				thMax = accu->GetXaxis()->GetBinCenter(i);
				ystMax = accu->GetYaxis()->GetBinCenter(j);
			}
	if (max > 0.) printf("track found at: th = %.2E, y* = %.2E\n", thMax, ystMax);


	
	// run the fake application
	app.Run();
	return 0;
}

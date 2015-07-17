TH1D *hlast;

void GetDraw(TFile *f, char *name, char *option, char *title, int color)
{
	TH1D *h = f->Get(name);
	h->SetLineColor(color);

	char buf[100];
	sprintf(buf, "%s, RMS = %.1E", title, h->GetRMS());
	h->SetTitle(buf);

	h->Draw(option);

	hlast = h;
}


void MakePlots()
{
	gStyle->SetOptStat(0);
	gStyle->SetCanvasColor(0);	
	gStyle->SetCanvasBorderMode(0);	
	gStyle->SetFrameBorderMode(0);	
	gStyle->SetOptTitle(0);
	gStyle->SetTitleBorderSize(1);
	gStyle->SetStatBorderSize(1);	
	gStyle->SetFuncWidth(1);
	gStyle->SetFuncColor(4);
	gStyle->SetStatW(0.25);
	gStyle->SetStatH(0.25);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetTitleYOffset(1.5);
	gStyle->SetTitleXOffset(1.2);
	gStyle->SetPadLeftMargin(0.10);
	gStyle->SetPadRightMargin(0.10);
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.08);
	gStyle->SetOptTitle(0);
	gStyle->ToggleEventStatus();
	

	TH2F *h2str;

	//----------------------------------------

	/*
	TFile *fNS = new TFile("../reconstructed/anal_roadsize_1535_1E4.root");

	
	//----------------------------------------
	TCanvas *c = new TCanvas("c1", "", 800, 400);
	c->Divide(2, 1);
	c->cd(1);
	//h2str = new TH2F("", "", 1, -4E-7, 4E-7, 1, 0., 2E3); h2str->Draw();
	GetDraw(fNS, "angles/thetaErr_x", "", "error in #vartheta_{x}", 1);
	GetDraw(fNS, "angles/thetaErr_y", "same", "error in #vartheta_{y}",  2);
	gPad->BuildLegend();

	c->cd(2);
	//h2str = new TH2F("", "", 1, -4E-7, 4E-7, 1, 0., 2E3); h2str->Draw();
	GetDraw(fNS, "vertex/vErr_x", "", "error in x*", 1);
	GetDraw(fNS, "vertex/vErr_y", "same", "error in y*",  2);
	gPad->BuildLegend();

	//----------------------------------------

	c = new TCanvas("c2", "", 400, 400);
	fNS->Get("t/tRes")->Draw("AP");

	//----------------------------------------

	gStyle->SetOptStat(1110);
	gStyle->SetStatW(0.35); gStyle->SetStatH(0.35);
	c = new TCanvas("c3", "", 800, 400);
	c->Divide(2, 1);
	c->cd(1); fNS->Get("statistics/yNdfHist")->Draw("");
	c->cd(2); fNS->Get("statistics/thYchi")->Draw("");

	c = new TCanvas("c4", "", 800, 400);
	c->Divide(2, 1);
	c->cd(1); fNS->Get("statistics/s2Min_y[2]")->Draw("");
	c->cd(2); fNS->Get("statistics/s2Min_y[4]")->Draw("");
	*/
	

	//----------------------------------------
	

	//TFile *fS = new TFile("../reconstructed/anal_roadsize_1535_1E4_fullSmear.root");
	TFile *fS = new TFile("../reconstructed/anal_reco_90_1E4_fullSmear.root");

	//----------------------------------------

	c = new TCanvas("c11", "", 800, 400);
	//gStyle->SetOptStat(0);
	c->Divide(2, 1);
	c->cd(1);
	//h2str = new TH2F("", "", 1, -4E-7, 4E-7, 1, 0., 2E3); h2str->Draw();
	GetDraw(fS, "angles/thetaErr_x", "", "error in #vartheta_{x}", 1);
	GetDraw(fS, "angles/thetaErr_y", "same", "error in #vartheta_{y}",  2);
	GetDraw(fS, "angles/thetaErr", "same", "error in #vartheta",  4);
	gPad->BuildLegend();

	c->cd(2);
	//h2str = new TH2F("", "", 1, -4E-7, 4E-7, 1, 0., 2E3); h2str->Draw();
	GetDraw(fS, "vertex/vErr_x", "", "error in x*", 1);
	GetDraw(fS, "vertex/vErr_y", "same", "error in y*",  2);
	gPad->BuildLegend();

	//----------------------------------------

	c = new TCanvas("c12", "", 400, 400);
	fS->Get("t/tResLog_y")->Draw("AP");

	//----------------------------------------

	c = new TCanvas("c13", "", 800, 400);
	c->Divide(2, 1);
	c->cd(1); fS->Get("right-left diff/dTh_y")->Draw("");
	c->cd(2); fS->Get("right-left diff/dVt_y")->Draw("");

	c = new TCanvas("c19", "", 300, 300);
	h = (TH1D *) fS->Get("one RP fit/120_vHist_y"); h->Draw(""); h->GetXaxis()->SetRangeUser(7., 7.5);
	
	
}

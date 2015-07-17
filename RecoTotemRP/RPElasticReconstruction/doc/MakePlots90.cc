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


void MakePlots90()
{
	gStyle->SetOptStat(1110);
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
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.08);
	//gStyle->SetTitleYOffset(1.3);
	//gStyle->SetTitleXOffset(1.3);
	

	TH2F *h2str;


	/*
	//----------------------------------------

	TFile *fNS = new TFile("../reconstructed/anal_reco_90_1E4.root");

	//----------------------------------------

	c = new TCanvas("c1", "", 800, 400);
	//gStyle->SetOptStat(0);
	c->Divide(2, 1);
	c->cd(1); h = (TH1D *)fNS->Get("angles/thetaErr_y"); h->Draw(""); h->GetXaxis()->SetTitleOffset(1.3);
	c->cd(2); h = (TH1D *)fNS->Get("angles/thetaErr_x"); h->Draw(""); h->GetXaxis()->SetTitleOffset(1.3);

	c = new TCanvas("c2", "", 800, 400);
	//gStyle->SetOptStat(0);
	c->Divide(2, 1);
	c->cd(1); h = (TH1D *)fNS->Get("vertex/vErr_y"); h->Draw(""); h->GetXaxis()->SetTitleOffset(1.3);
	c->cd(2); h = (TH1D *)fNS->Get("vertex/vErr_x"); h->Draw(""); h->GetXaxis()->SetTitleOffset(1.3);

	//----------------------------------------

	c = new TCanvas("c3", "", 400, 400);
	g = (TGraph *)fNS->Get("t/tRes_y"); g->Draw("AP"); g->GetXaxis()->SetTitleOffset(1.3); g->GetYaxis()->SetTitleOffset(1.5);

	//----------------------------------------

	//gStyle->SetOptStat(1110);
	//gStyle->SetStatW(0.35); gStyle->SetStatH(0.35);
	c = new TCanvas("c4", "", 800, 400);
	c->Divide(2, 1);
	c->cd(1); h = (TH1D *)fNS->Get("statistics/yNdfHist"); h->Draw(""); h->GetXaxis()->SetTitleOffset(1.3);
	c->cd(2); h = (TH1D *)fNS->Get("statistics/s2Min_y[2]"); h->Draw(""); h->GetXaxis()->SetTitleOffset(1.3);


	c = new TCanvas("c5", "", 400, 400);
	c->cd(3); h = (TH1D *)fNS->Get("statistics/thYchi"); h->Draw(""); h->GetXaxis()->SetTitleOffset(1.3);
	h->Fit("gaus");
	*/

	//----------------------------------------

	TFile *fS = new TFile("../reconstructed/anal_reco_90_1E4_fullSmear.root");

	//----------------------------------------

	c = new TCanvas("c11", "", 800, 400);
	//gStyle->SetOptStat(0);
	c->Divide(2, 1);
	c->cd(1); h = (TH1D *)fS->Get("angles/thetaErr_y"); h->Draw(""); h->GetXaxis()->SetTitleOffset(1.3);
	c->cd(2); h = (TH1D *)fS->Get("vertex/vErr_x"); h->Draw(""); h->GetXaxis()->SetTitleOffset(1.3);

	//----------------------------------------

	c = new TCanvas("c12", "", 800, 400);
	c->Divide(2, 1);
	c->cd(1); g = (TGraph *)fS->Get("t/tRes_y"); g->Draw("AP"); g->GetYaxis()->SetRangeUser(0., 45.); g->GetXaxis()->SetTitleOffset(1.3); g->GetYaxis()->SetTitleOffset(1.3);
	c->cd(2); g = (TGraph *)fS->Get("t/tRes_x"); g->Draw("AP"); g->GetYaxis()->SetRangeUser(0., 45.); g->GetXaxis()->SetTitleOffset(1.3); g->GetYaxis()->SetTitleOffset(1.3);

	c = new TCanvas("c12a", "", 400, 400);
	c->cd(3); g = (TGraph *)fS->Get("t/tRes"); g->Draw("AP"); g->GetYaxis()->SetRangeUser(0., 45.); g->GetXaxis()->SetTitleOffset(1.3); g->GetYaxis()->SetTitleOffset(1.3);

	//----------------------------------------

	c = new TCanvas("c13", "", 800, 400);
	c->Divide(2, 1);
	c->cd(1); h = (TH1D *)fS->Get("right-left diff/dTh_y"); h->Draw(""); h->GetXaxis()->SetTitleOffset(1.3);
	c->cd(2); h = (TH1D *)fS->Get("right-left diff/dVt_x"); h->Draw(""); h->GetXaxis()->SetTitleOffset(1.3);
}

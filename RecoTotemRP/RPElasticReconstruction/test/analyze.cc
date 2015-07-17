#include <stdio.h>
#include <iostream>

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPTrackCandidateCollection.h"
#include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDetTrigger.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "DataFormats/TotemRPDataTypes/interface/RPRecoHit.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecoElasticEvent.h"
#include "DataFormats/Provenance/interface/ParameterSetBlob.h"
#include "DataFormats/Provenance/interface/HashedTypes.h"
#include "DataFormats/Provenance/interface/Hash.h"

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


//----------------------------------------------------------------------------------------------------

void ShowTreeStructure(TTree *events)
{
	// brach and alias info
	printf("%i events\n", events->GetEntriesFast()); 
	printf("Branches: \n");
	TObjArray *ba = events->GetListOfBranches();
	for (unsigned int i = 0; i < ba->GetEntriesFast(); i++) {
		TBranch *b = (TBranch *) ba->At(i);
		printf("\t%i: '%s'\n", i, b->GetName());
	}

	printf("Aliases:\n");
	TList *al = (TList *) events->GetListOfAliases();
	for (unsigned int i = 0; i < al->GetEntries(); i++) {
		TObject *o = al->At(i);
		printf("\t%i: '%s' -> '%s'\n", i, o->GetName(), o->GetTitle());
	}
}

//----------------------------------------------------------------------------------------------------

void ExamineMetaData(TTree *md)
{
	TBranch *b = md->GetBranch("ParameterSetMap");

	typedef edm::Hash<edm::ParameterSetType> neco_type;
	typedef std::map< neco_type, edm::ParameterSetBlob > pset_type;
	pset_type parMap;
	b->SetAddress(&parMap);

	for (int i = 0; i < b->GetEntries(); i++) {
		b->GetEntry(i);
		std::cout << i << ": size = " << parMap.size() << std::endl;
		for (pset_type::iterator it = parMap.begin(); it != parMap.end(); ++it) {
			std::cout << "hash: " << it->first << std::endl;
			std::cout << "blob: " << it->second << std::endl;
		}
	}
}

//----------------------------------------------------------------------------------------------------

void AddFunction(TH1D* h, TF1 *f)
{
	TList *l = h->GetListOfFunctions();
	if (!l) {
		h->Fit("pol1", "Q");
		l = h->GetListOfFunctions();
		l->Clear();
	}
	l->Add(f);
}

//----------------------------------------------------------------------------------------------------

TGraphErrors* FitAndDrawErrProfile(TProfile *p, bool doFit = true, bool logScale = false)
{
	TGraphErrors* g = new TGraphErrors();

	double sumYF = 0., sumFF = 0.;
	for (int i = 1; i <= p->GetNbinsX() ; i++) {
		if (p->GetBinEntries(i) < 1) continue;
		// get and convert values
		double x = p->GetBinCenter(i);
		double y = sqrt(p->GetBinContent(i)) * 100.;
		double sy = p->GetBinError(i) / 2. / sqrt(p->GetBinContent(i)) * 100.;
		//printf("%.2E, %.2E +- %.2E\t(%f) \n", x, y, sy, p->GetBinEntries(i));

		// set values to the graph
		g->SetPoint(g->GetN(), x, y);
		g->SetPointError(g->GetN() - 1, 0., sy);

		if (p->GetBinEntries(i) < 3) continue;
		double xx = (logScale) ? pow(10., x) : x;
		sumYF += y /sqrt(xx) /sy/sy;
		sumFF += 1./xx /sy/sy;
	}

	double A = sumYF/sumFF;

	// set properties
	g->GetXaxis()->SetTitle(p->GetXaxis()->GetTitle());
	g->GetYaxis()->SetTitle(p->GetYaxis()->GetTitle());
	g->SetMarkerStyle(2); g->SetMarkerSize(1.);

	if (doFit) {
		// just a tric to initialize list of functions
		g->Fit("pol1", "Q");
		g->GetListOfFunctions()->Clear();

		// now the real stuff
		TF1 *fitFun;
		//TF1 *blaf = new TF1("blaf", "[0]/sqrt(pow(10, x))", -3., 1.);
		if (logScale) fitFun = new TF1("userLog", "[0]/sqrt(pow(10, x))", -3., -0.);
		else fitFun = new TF1("user", "[0]/sqrt(x)");
		fitFun->SetParameter(0, A);
		//fitFun->Draw("same");
		g->GetListOfFunctions()->Add(fitFun);
		printf("%s: A = %.2E\n", p->GetName(), A);
	}

	return g;
}

//----------------------------------------------------------------------------------------------------

struct optFun {
	double Lx, Ly, vx, vy;
	optFun(double _lx = 0., double _ly = 0., double _vx = 0., double _vy = 0.) : Lx(_lx), Ly(_ly), vx(_vx), vy(_vy) {}
};

std::map<unsigned short, optFun> optics;

void InitOptics_90()
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

void InitOptics_1535()
{
	optics[0] = optFun(0., 0., 1., 1.);
	optics[20] = optFun(-112.5, -248.8, 0.05947, 0.01992);
	optics[21] = optFun(-112.5, -248.8, 0.05947, 0.01992);
	optics[22] = optFun(-111.4, -250.8, 0.05485, 0.01827);
	optics[23] = optFun(-99.73, -270.7, 0.008973, 0.001909);
	optics[24] = optFun(-98.55, -272.7, 0.004355, 0.0001949);
	optics[25] = optFun(-98.55, -272.7, 0.004355, 0.0001949);
	optics[120] = optFun(112.5, 248.8, 0.05947, 0.01992);
	optics[121] = optFun(112.5, 248.8, 0.05947, 0.01992);
	optics[122] = optFun(111.4, 250.8, 0.05485, 0.01827);
	optics[123] = optFun(99.73, 270.7, 0.008973, 0.001909);
	optics[124] = optFun(98.55, 272.7, 0.004355, 0.0001949);
	optics[125] = optFun(98.55, 272.7, 0.004355, 0.0001949);
}


//----------------------------------------------------------------------------------------------------

struct OneRPFitResult {
	TH1D *thHist_x, *thHist_y;
	TH1D *vHist_x, *vHist_y;
	TH2D *hitHist;
	TGraph *thCorr_x, *thCorr_y;
	OneRPFitResult(int id = 0)
	{
		char buf[50];
		sprintf(buf, "%i_thCorr_x", id); thCorr_x = new TGraph(); thCorr_x->SetName(buf); thCorr_x->SetTitle(";#vartheta_{x}^{sim};#vartheta_{x}^{reco}");
		sprintf(buf, "%i_thCorr_y", id); thCorr_y = new TGraph(); thCorr_y->SetName(buf); thCorr_y->SetTitle(";#vartheta_{y}^{sim};#vartheta_{y}^{reco}");
		sprintf(buf, "%i_thHist_x", id); thHist_x = new TH1D(buf, ";#vartheta_{x}^{reco}", 1000, -5E-3, 5E-3);
		sprintf(buf, "%i_thHist_y", id); thHist_y = new TH1D(buf, ";#vartheta_{y}^{reco}", 1000, -5E-3, 5E-3);
		sprintf(buf, "%i_vHist_x", id); vHist_x = new TH1D(buf, ";x^{reco}   (mm)", 30000, -15E0, 15E0);
		sprintf(buf, "%i_vHist_y", id); vHist_y = new TH1D(buf, ";y^{reco}   (mm)", 30000, -15E0, 15E0);
		sprintf(buf, "%i_hitHist", id); hitHist = new TH2D(buf, ";y^{reco}   (mm)", 1000, -10E0, 10E0, 1000, -1E1, 1E1);
	}

	~OneRPFitResult()
	{
		delete thCorr_x; delete thCorr_y;
		delete thHist_x; delete thHist_y;
	}

	void Fill(RPFittedTrack ft, double th_x, double th_y)
	{
		thCorr_x->SetPoint(thCorr_x->GetN(), th_x, ft.GetTx());
		thCorr_y->SetPoint(thCorr_y->GetN(), th_y, ft.GetTy());
		thHist_x->Fill(ft.GetTx());
		thHist_y->Fill(ft.GetTy());
		vHist_x->Fill(ft.X0());
		vHist_y->Fill(ft.Y0());
		hitHist->Fill(ft.X0(), ft.Y0());
	}

	void Write()
	{
		thCorr_x->Write(); thCorr_y->Write();
		thHist_x->Write(); thHist_y->Write();
		vHist_x->Write(); vHist_y->Write();
		hitHist->Write();
	}
};

//----------------------------------------------------------------------------------------------------

TBranch* GetBranchByName(TTree *tree, const char *name, void *addr)
{
	TBranch *temp = tree->GetBranch(name); 
	if (!temp) { printf("ERROR: branch `%s' doesn't exit\n", name); exit(1); }
	temp->SetAddress(addr);
	return temp;
}

TBranch* GetBranchByAlias(TTree *tree, const char *alias, void *addr)
{
	const char *name = tree->GetAlias(alias);
	if (!name) { printf("ERROR: alias `%s' doesn't exit\n", alias); exit(1); }
	TBranch *temp = tree->GetBranch(name); 
	if (!temp) { printf("ERROR: branch `%s' doesn't exit\n", name); exit(1); }
	temp->SetAddress(addr);
	return temp;
}


//----------------------------------------------------------------------------------------------------

struct EffStruct {
	unsigned int tot, trig;
	unsigned int rOK, rNoRoad, rNoGoodRoad, rBadVertexX, rBadVertexY, rBadAngleX, rBadAngleY;

	TH1D *th_x, *th_y;

	EffStruct() : tot(0), trig(0), rOK(0), rNoRoad(0), rNoGoodRoad(0), rBadVertexX(0), rBadVertexY(0), rBadAngleX(0), rBadAngleY(0)
	{
		th_x = new TH1D("", ";#Delta#vartheta_{x}", 200, -3E-5, 3E-5);
		th_y = new TH1D("", ";#Delta#vartheta_{y}", 200, -3E-6, 3E-6);
	}

	void Update(bool triggered, RPRecoElasticEvent &ee, double dthx, double dthy)
	{
		tot++;
		if (triggered) trig++;
		if (ee.status == ee.sOK) rOK++;
		if (ee.status == ee.sNoRoad) rNoRoad++;
		if (ee.status == ee.sNoGoodRoad) rNoGoodRoad++;
		if (ee.status == ee.sBadVertexX) rBadVertexX++;
		if (ee.status == ee.sBadVertexY) rBadVertexY++;
		if (ee.status == ee.sBadAngleX) { rBadAngleX++; th_x->Fill(dthx); }
		if (ee.status == ee.sBadAngleY) { rBadAngleY++; th_y->Fill(dthy); }
	}
};

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: runAnalysis <options> <input file.root> <options>\n");
	printf("OPTIONS:\n");
	printf("\t-q\tTo quit ROOT in the end.\n");
	printf("\t-s\tTo show input file structure and quit.\n");
	printf("\t-v\tVerbose operation (medium level)\n");
	printf("\t-v2\tVerbose operation (full)\n");
	printf("\t-bad(vx|vy|ax|ay)\tVerbose operation for events with incompatibile vertex (v) or angular (a) fits. x and y for projections.\n");
	printf("\t-badroad\tVerbose operation for events with no good road.\n");
	exit(1);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int main(unsigned int argc, char *argv[])
{
	using namespace std;
	using namespace edm;

	// process command line parameters
	unsigned char inputIndex = 0;
	bool quitInTheEnd = false;
	bool printStructure = false, printParams = false;
	unsigned char verbGl = 0;
	bool fbvx = false, fbvy = false, fbax = false, fbay = false, fbroad = false;
	for (unsigned int i = 1; i < argc; i++) {
		if (!strncmp(argv[i], "-q", 2)) { quitInTheEnd = true; continue; }
		if (!strncmp(argv[i], "-s", 2)) { printStructure = true; continue; }
		if (!strcmp(argv[i], "-v")) { verbGl = 1; continue; }
		if (!strcmp(argv[i], "-v2")) { verbGl = 2; continue; }
		if (!strcmp(argv[i], "-badvx")) { fbvx = true; continue; }
		if (!strcmp(argv[i], "-badvy")) { fbvy = true; continue; }
		if (!strcmp(argv[i], "-badax")) { fbax = true; continue; }
		if (!strcmp(argv[i], "-baday")) { fbay = true; continue; }
		if (!strcmp(argv[i], "-badroad")) { fbroad = true; continue; }

		if (argv[i][0] != '-') { inputIndex = i; continue; }

		printf("Unrecognized parameter `%s'.\n", argv[i]);
		PrintUsage();
		return 5;
	}

	if (inputIndex == 0) {
		printf("You must specify input file\n");
		PrintUsage();
	}

	// create fake application
	int fakeC = 2; char *fakeV[2] = {"prdlajs", "-l"};
	TRint* app;
	if (!quitInTheEnd) app = new TRint("application", &fakeC, fakeV, 0, true);
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

	// enable automatic loading of dictionaries
	AutoLibraryLoader::enable();

	//InitOptics_1535();
	InitOptics_90();

	// open file and get event tree
	TFile *theFile = new TFile(argv[inputIndex]);
	if (theFile->IsZombie()) { printf("No such a file `%s'.\n", argv[inputIndex]); return 2; }
	TTree *events = (TTree*)theFile->Get("Events"); 
	if (!events) { printf("No CMSSW data in the file.\n"); return 3; }
	if (printStructure) { ShowTreeStructure(events); return 0; }
	TTree *metaData = (TTree*)theFile->Get("MetaData"); 
	if (printParams) { ExamineMetaData(metaData); }

	/// smeared or not smeared source
	bool smearedSource = false;

	// =================================== prepare branches
	edm::HepMCProduct srcPr; 
	TBranch *srcBr = events->GetBranch(events->GetAlias("source")); 
	srcBr->SetAddress(&srcPr); 

	edm::HepMCProduct srcSmearPr;
	TBranch *srcSmearBr;
	if (events->GetAlias("EnergyVertexSmeared")) {
		smearedSource = true;
		srcSmearBr = events->GetBranch(events->GetAlias("EnergyVertexSmeared")); 
		srcSmearBr->SetAddress(&srcSmearPr);
	}

	RPFittedTrackCollection fitTrCol;
	TBranch *fitTrColBr = events->GetBranch(events->GetAlias("RPSingleTrackCandCollFit")); 
	fitTrColBr->SetAddress(&fitTrCol);

	/*
	RPTrackCandidateCollection candTrCol;
	TBranch *candTrColBr = events->GetBranch(events->GetAlias("RPSinglTrackCandFind")); 
	candTrColBr->SetAddress(&candTrCol);

	DetSetVector<RPStripDigi> digiCol;
	TBranch *digiColBr = events->GetBranch(events->GetAlias("RPSiDetDigitizer")); 
	digiColBr->SetAddress(&digiCol);

	DetSetVector<RPRecoHit> recoCol;
	TBranch *recoColBr = events->GetBranch(events->GetAlias("RPHecoHitProd")); 
	recoColBr->SetAddress(&recoCol);
	*/

	DetSetVector<RPDetTrigger> trigCol;
	TBranch *trigColBr = GetBranchByName(events, "RPStripDigiedmDetSetVector_RPSiDetDigitizer__ElasticSimulation.obj", &trigCol); 

	RPRecoElasticEvent elasticReco;
	TBranch *elRecoBr = GetBranchByAlias(events, "ElasticReconstruction", &elasticReco);

	// =================================== prepare histograms and outputs
	TH1::AddDirectory(kFALSE);

	map<unsigned int, OneRPFitResult*> oneRPFitRes;
	oneRPFitRes[120] = new OneRPFitResult(120);
	oneRPFitRes[122] = new OneRPFitResult(122);
	oneRPFitRes[125] = new OneRPFitResult(125);

	TGraph *thetaYcorr = new TGraph(); thetaYcorr->SetTitle(";#theta_{y}^{sim};#theta_{y}^{reco}");
	TGraph *thetaXcorr = new TGraph(); thetaXcorr->SetTitle(";#theta_{x}^{sim};#theta_{x}^{reco}");

	TGraph *yStCorr = new TGraph(); yStCorr->SetTitle(";y*^{sim}   (m);y*^{reco}   (m)");
	TGraph *xStCorr = new TGraph(); xStCorr->SetTitle(";x*^{sim}   (m);x*^{reco}   (m)");

	TH1D *thHistAcc_y = new TH1D("thHistAcc_y", ";#vartheta_{y}^{reco}   (rad)", 20000, -1E-4, 1E-4);
	TH1D *thHistFull_y = new TH1D("thHistFull_y", ";#vartheta_{y}^{reco}   (rad)", 20000, -1E-4, 1E-4);

	TH1D *thetaYerr = new TH1D("thetaYerr", ";#vartheta_{y}^{reco} - #vartheta_{y}^{sim}   (rad)", 100, -3E-7, 3E-7);
	TH1D *thetaXerr = new TH1D("thetaXerr", ";#vartheta_{x}^{reco} - #vartheta_{x}^{sim}   (rad)", 100, -3E-5, 3E-5);
	TH1D *thetaErr = new TH1D("thetaErr", ";#vartheta^{reco} - #vartheta^{sim}   (rad)", 100, -1E-6, 1E-6);

	TH1D *phiHist_reco = new TH1D("phiHist_reco", ";#varphi", 100, -M_PI, +M_PI);
	TH1D *phiHist_sim = new TH1D("phiHist_sim", ";#varphi", 100, -M_PI, +M_PI);
	TH1D *thetaHist_reco = new TH1D("thetaHist_reco", ";#vartheta", 100, 0., 1E-4);
	TH1D *thetaHist_sim = new TH1D("thetaHist_sim", ";#vartheta", 100, 0., 1E-4);
	TH1D *tHist_reco = new TH1D("tHist_reco", ";t   (GeV^{2})", 100, 0., 0.5);
	TH1D *tHist_sim = new TH1D("tHist_sim", ";t   (GeV^{2})", 100, 0., 0.5);

	TH1D *xStErr = new TH1D("xStErr", ";x*^{reco} - x*^{sim}   (m)", 100, -1E-5, 1E-5);
	TH1D *yStErr = new TH1D("yStErr", ";y*^{reco} - y*^{sim}   (m)", 100, -1E-3, 1E-3);

	TProfile *tYrelErr = new TProfile("tYrelErr", ";|t_{y}^{sim}|   (GeV^{2});#sigma(#Deltat_{y}) / t_{y}^{}   (%)",	50, 0., 5E-1);
	TProfile *tXrelErr = new TProfile("tXrelErr", ";|t_{x}^{sim}|   (GeV^{2});#sigma(#Deltat_{x}) / t_{x}^{}   (%)",	50, 0., 5E-1);
	TProfile *tRelErr = new TProfile("tRelErr", ";|t^{sim}|   (GeV^{2});#sigma(#Deltat) / t^{}   (%)", 				50, 0., 5E-1);
	TProfile *tRelErrLog = new TProfile("tRelErrLog", ";log_{10}(|t^{sim}| / GeV^{2});#sigma(#Deltat) / t   (%)", 	50, -3., 1.);
	TProfile *tYRelErrLog = new TProfile("tYRelErrLog", ";log_{10}(|t_{y}^{sim}| / GeV^{2});#sigma(#Deltat_{y}) / t_{y}   (%)", 50, -3., 1.);

	TH1D* dTh_x = new TH1D("dTh_x", ";#Delta_{R-L}#vartheta_{x}   (rad)", 500, -2E-6, 2E-6);
	TH1D* dTh_y = new TH1D("dTh_y", ";#Delta_{R-L}#vartheta_{y}   (rad)", 5000, -10E-6, 10E-6);
	TH1D* dVt_x = new TH1D("dVt_x", ";#Delta_{R-L}x*   (m)", 100, -1E-4, 1E-4);
	TH1D* dVt_y = new TH1D("dVt_y", ";#Delta_{R-L}y*   (m)", 100, -1E-2, 1E-2);

	TH1D *thYchi = new TH1D("thYchi", "normalized error;(#vartheta_{y}^{reco} - #vartheta_{y}^{sim}) / #sigma(#vartheta_{y})^{reco}", 100, -10., 10.);
	TH1D *thXchi = new TH1D("thXchi", ";(#vartheta_{x}^{reco} - #vartheta_{x}^{sim}) / #sigma(#vartheta_{x})^{reco}", 100, -10., 10.);

	TH1D *yNdfHist = new TH1D("yNdfHist", "number of degrees of freedom;ndf_{y}", 10, -0.5, 9.5);
	TH1D *xNdfHist = new TH1D("xNdfHist", ";ndf_{x}", 10, -0.5, 9.5);

	TH1D *s2Min_x[5], *s2Min_y[5];
	for (int i = 0; i < 5; i++) {
		char buf[20];
		sprintf(buf, "s2Min_x[%i]", i); s2Min_x[i] = new TH1D(buf, ";S^{2}_{min}^{x}", 150, 0., 20.);
		sprintf(buf, "s2Min_y[%i]", i); s2Min_y[i] = new TH1D(buf, ";S^{2}_{min}^{y}", 150, 0., 20.);
	}

	TH1D *yResHist = new TH1D("yResHist", "normalized residua;(y - #hat{y}) / #sigma_{y}", 100, -1E1, 1E1);
	TH1D *xResHist = new TH1D("xResHist", ";(x - #hat{x}) / #sigma_{x}", 100, -1E1, 1E1);
	TH1D *yResHist0 = new TH1D("yResHist0", ";(y - #hat{y}) / #sigma_{y}", 500, -1E1, 1E1);

	TH1D *roadCountHist = new TH1D("roadCountHist", ";number of roads", 10, -0.5, 9.5);

	map<signed int, EffStruct> counters;

	TGraph *bla = new TGraph();
	TH1D *blaHist = new TH1D("blaHist", "", 200000, -1E3, 1E3);

	// =================================== loop over events

	for(unsigned int i = 0; i < events->GetEntries(); i++) { 
		// ------------------------------ GET EVENTS
		//.events->GetEntry(i); //. doesn't work !
		srcBr->GetEntry(i);
		if (smearedSource) srcSmearBr->GetEntry(i);
		fitTrColBr->GetEntry(i);
		//candTrColBr->GetEntry(i);
		//digiColBr->GetEntry(i);
		trigColBr->GetEntry(i);
		//recoColBr->GetEntry(i);
		elRecoBr->GetEntry(i);
		HepMC::GenEvent *srcEv = (HepMC::GenEvent *) srcPr.GetEvent(); 
		HepMC::GenEvent *smearedSrcEv =(HepMC::GenEvent *) srcSmearPr.GetEvent();

		// ------------------------------ VERBOSITY LEVEL
		unsigned char verbose = verbGl;
		if (fbvx && elasticReco.status == elasticReco.sBadVertexX) verbose = 2;
		if (fbvy && elasticReco.status == elasticReco.sBadVertexY) verbose = 2;
		if (fbax && elasticReco.status == elasticReco.sBadAngleX) verbose = 2;
		if (fbay && elasticReco.status == elasticReco.sBadAngleY) verbose = 2;
		if (fbroad && elasticReco.status == elasticReco.sNoGoodRoad) verbose = 2;

		// print info 
		if ((i % 1000) == 0 || verbose)
			printf("[0;32mEVENT[0m	MC=%i,   in order=%i\n", srcEv->event_number(), i);

		if (verbose) {
			cout << "[0;34mGENERATOR ";
			if (smearedSource) cout << "(SMEARED)[0m" << endl;
			else cout << "(NOT SMEARED)[0m" << endl;
		}

		// ------------------------------ VERTEX POSITION (only one vertex assumed)
		HepMC::GenEvent::vertex_iterator vI = (smearedSource) ? smearedSrcEv->vertices_begin() :  srcEv->vertices_begin();
		HepMC::FourVector vertex = (*vI)->position();		// in mm
		double x_st = vertex.x() * 1E-3;
		double y_st = vertex.y() * 1E-3;
		if (verbose) printf("vertex (in um): x = %.2E, y = %.2E    (z = %.2E mm)\n", x_st * 1E6, y_st * 1E6, vertex.z());


		// ------------------------------ PARTICLE INFORMATION
		HepMC::GenParticle *particle = NULL;
		if (srcEv->signal_process_id() == 91)		particle = srcEv->barcode_to_particle(3);	// elastic scattering
		if (srcEv->signal_process_id() == 123456)	particle = (* srcEv->particles_begin());	// background
		if (!particle) {
			cout << "No principal particle found. Event skipped." << endl;
			continue;
		}
		HepMC::FourVector p = particle->momentum();		// in GeV
		signed int ID = particle->pdg_id();
		if (verbose) cout << "ID = " << ID << endl;
		if (verbose) cout << "momentum = " << p << endl;
		double p_reco = 7E3;
		double theta_y = p.y() / p.z();
		double theta_x = p.x() / p.z();
		double theta = sqrt(theta_x*theta_x + theta_y*theta_y);
		double ty = theta_y*theta_y * p_reco*p_reco;
		double tx = theta_x*theta_x * p_reco*p_reco;
		double t = tx + ty;
		if (verbose) printf("t = %.2E, t_x = %.2E, t_y = %.2E\n", t, tx, ty);
		if (verbose) cout << "tan(theta) = " << theta << endl;
		if (verbose) cout << "tan(theta_x) = " << theta_x << endl;
		if (verbose) cout << "tan(theta_y) = " << theta_y << endl;
		phiHist_sim->Fill(atan2(theta_y, theta_x));
		thetaHist_sim->Fill(theta);
		tHist_sim->Fill(t);
		thHistFull_y->Fill(theta_y);


		// ------------------------------ TRIGGER ANALYZIS
		bool leftHit = false, rightHit = false;
		for (DetSetVector<RPDetTrigger>::const_iterator it = trigCol.begin(); it != trigCol.end(); ++it) {
			TotRPDetId detID(it->detId());
			unsigned int RPId = detID.DetectorDecId() / 10;
			if (it->size() && RPId < 100) leftHit = true;
			if (it->size() && RPId >= 100) rightHit = true;
		}
		bool triggered = leftHit && rightHit;	

		int nr = 0, nl = 0;
		double sr = 0, sl = 0;

		// ------------------------------ ONE RP FIT 
		if (verbose > 1) cout << "[0;34mLOCAL TRACK FITS[0m (mm)" << endl;
		for (std::map<RPId, RPFittedTrack>::const_iterator it = fitTrCol.begin(); it != fitTrCol.end(); ++it) {
			if (!it->second.IsValid()) continue;
			if (oneRPFitRes.find(it->first) != oneRPFitRes.end())
				oneRPFitRes[it->first]->Fill(it->second, theta_x, theta_y);
			if (verbose > 1)
				printf("\t%3i\tx=%8.3f\ty=%8.3f\t\tsi_x=%.2E, si_y=%.2E\n", it->first, it->second.X0(), it->second.Y0(), 
						it->second.X0Sigma(), it->second.Y0Sigma());


			double spatialSpacing = 66E-3 / sqrt(2.);

			if (it->first >= 100) {
				nr++;
				sr += (it->second.Y0() / spatialSpacing);
			} else {
				nl++;
				sl += (it->second.Y0() / spatialSpacing);
			}
		}

		blaHist->Fill(sr);

		// ------------------------------ ELASTIC RECONSTRUCTION REPORT
		if (verbose) cout << "[0;34mELASTIC RECONSTRUCTION[0m" << endl;
		if (verbose) cout << "status = " << elasticReco.status << endl;
		if (verbose > 100) {
			printf("roads size = %i\n", elasticReco.roads.size());
			for (unsigned int i = 0; i < elasticReco.roads.size(); i++) {
				printf("\t[%i] at x = %.2E, y = %.2E :  ", i, elasticReco.roads[i].centerX(), elasticReco.roads[i].centerY());
				for (int j = 0; j < elasticReco.roads[i].members.size(); j++)
					printf("%i, ", elasticReco.roads[i].members[j]);
				printf("\n");
			}
	
			printf("preferred road = %i\n", elasticReco.preferredRoad);

			if (elasticReco.preferredRoad >= 0) {
				printf("---------------------------------------------------------------------------------------------------------\n");
				printf("   fit | projection | theta (rad) | theta error |  vertex (m) |  vertex err |         ndf |  S2_min/ndf |\n");
				printf("---------------------------------------------------------------------------------------------------------\n");
				RPRecoElasticEvent::fit_type &f = elasticReco.leftFit;
				printf("  left |          x |%12.2E |%12.2E |%12.2E |%12.2E |%12i |%12.2E |\n", f.th_x, f.si_th_x, f.x, f.si_x, f.ndf_x, f.s2minPerDf_x());
				printf("                  y |%12.2E |%12.2E |%12.2E |%12.2E |%12i |%12.2E |\n", f.th_y, f.si_th_y, f.y, f.si_y, f.ndf_y, f.s2minPerDf_y());
				printf("---------------------------------------------------------------------------------------------------------\n");
				f = elasticReco.rightFit;
				printf(" right |          x |%12.2E |%12.2E |%12.2E |%12.2E |%12i |%12.2E |\n", f.th_x, f.si_th_x, f.x, f.si_x, f.ndf_x, f.s2minPerDf_x());
				printf("                  y |%12.2E |%12.2E |%12.2E |%12.2E |%12i |%12.2E |\n", f.th_y, f.si_th_y, f.y, f.si_y, f.ndf_y, f.s2minPerDf_y());
				printf("---------------------------------------------------------------------------------------------------------\n");
				f = elasticReco.globalFit;
				printf("global |          x |%12.2E |%12.2E |%12.2E |%12.2E |%12i |%12.2E |\n", f.th_x, f.si_th_x, f.x, f.si_x, f.ndf_x, f.s2minPerDf_x());
				printf("                  y |%12.2E |%12.2E |%12.2E |%12.2E |%12i |%12.2E |\n", f.th_y, f.si_th_y, f.y, f.si_y, f.ndf_y, f.s2minPerDf_y());
				printf("---------------------------------------------------------------------------------------------------------\n");
			}
		}

		double SpatialUncertaintyCorrection = 1.;
		// ------------------------------ ELASTIC RECONSTRUCTION PLOTS
		if (elasticReco.isValid()) {
			RPRecoElasticEvent::fit_type &f = elasticReco.result;
			double th_reco = sqrt(f.th_y*f.th_y + f.th_x*f.th_x);
			double tx_reco = f.th_x*f.th_x * p_reco*p_reco;
			double ty_reco = f.th_y*f.th_y * p_reco*p_reco;
			double t_reco = tx_reco + ty_reco;

			thetaYcorr->SetPoint(thetaYcorr->GetN(), theta_y, f.th_y);
			thetaXcorr->SetPoint(thetaXcorr->GetN(), theta_x, f.th_x);

			yStCorr->SetPoint(yStCorr->GetN(), y_st, f.y);
			xStCorr->SetPoint(xStCorr->GetN(), x_st, f.x);

			thHistAcc_y->Fill(theta_y);

			thetaYerr->Fill(f.th_y - theta_y);
			thetaXerr->Fill(f.th_x - theta_x);
			thetaErr->Fill(th_reco - theta);

			yStErr->Fill(f.y - y_st);
			xStErr->Fill(f.x - x_st);

			thXchi->Fill((f.th_x - theta_x) / (f.si_th_x * SpatialUncertaintyCorrection));
			thYchi->Fill((f.th_y - theta_y) / (f.si_th_y * SpatialUncertaintyCorrection));

			tYrelErr->Fill(fabs(ty), (ty_reco - ty) * (ty_reco - ty) / ty / ty);
			tXrelErr->Fill(fabs(tx), (tx_reco - tx) * (tx_reco - tx) / tx / tx);
			tRelErr->Fill(fabs(t), (t_reco - t) * (t_reco - t) / t / t);
			tYRelErrLog->Fill(log10(fabs(ty)), (ty_reco - ty) * (ty_reco - ty) / ty / ty);
			tRelErrLog->Fill(log10(fabs(t)), (t_reco - t) * (t_reco - t) / t / t);

			phiHist_reco->Fill(atan2(f.th_y, f.th_x));
			thetaHist_reco->Fill(th_reco);
			tHist_reco->Fill(t_reco);
		}

		double dthx = elasticReco.rightFit.th_x - elasticReco.leftFit.th_x;
		double dthy = elasticReco.rightFit.th_y - elasticReco.leftFit.th_y;
		double dvtx = elasticReco.rightFit.x - elasticReco.leftFit.x;
		double dvty = elasticReco.rightFit.y - elasticReco.leftFit.y;
		if (elasticReco.status != elasticReco.sNoGoodRoad && elasticReco.status != elasticReco.sNoRoad && elasticReco.result.ndf_x < 200) {
			dTh_x->Fill(dthx); dTh_y->Fill(dthy);
			dVt_x->Fill(dvtx); dVt_y->Fill(dvty);

			bla->SetPoint(bla->GetN(), dthy, (sr/nr + sl/nl));
		}

		// ------------------------------ ELASTIC RECONSTRUCTION RESIDUALS
		if (elasticReco.isValid()) {
			double s2min_y = 0, s2min_x = 0;
			int nHits = 0;
			RPRecoElasticEvent::fit_type &f = elasticReco.result;
			// RP residua
			for (map<RPId, RPFittedTrack>::const_iterator it = fitTrCol.begin(); it != fitTrCol.end(); ++it) {
				// skip invalid
				if (!it->second.IsValid()) continue;

				unsigned short RP;
				double x, y, sx, sy;

				RP = it->first;
				x = it->second.X0()*1E-3; y = it->second.Y0()*1E-3;
				sx = sy = 66E-6/sqrt(12) * SpatialUncertaintyCorrection;	//it->second.Y0Sigma()*1E-3;

				optFun &of = optics[RP];
				double y_res = (y - f.th_y * of.Ly - f.y * of.vy) / sy;
				double x_res = (x - f.th_x * of.Lx - f.x * of.vx) / sx;

				s2min_y += y_res * y_res;
				s2min_x += x_res * x_res;
				nHits++;

				yResHist->Fill(y_res);
				xResHist->Fill(x_res);
			}

			// degrees of freedom histograms
			xNdfHist->Fill(nHits - 2);
			yNdfHist->Fill(nHits - 2);

			// residual sum histograms
			if (nHits - 2 < 5) {
				s2Min_x[nHits-2]->Fill(s2min_x);
				s2Min_y[nHits-2]->Fill(s2min_y);
			}
		}

		// ------------------------------ INCREMENT COUNTERS
		counters[ID].Update(triggered, elasticReco, dthx, dthy);
		roadCountHist->Fill(elasticReco.roads.size());
	}


	// =================================== after processing
	TF1 *chiSq[5];
	chiSq[0] = NULL;
	chiSq[1] = new TF1("chiSq1", "exp(-x / 2) / sqrt(2 * pi * x)", 0., 20.); chiSq[1]->SetLineColor(8); chiSq[1]->SetTitle("#chi^2(#nu = 1)");
	chiSq[2] = new TF1("chiSq2", "exp(-x / 2) / 2", 0., 20.); chiSq[2]->SetLineColor(2);
	chiSq[3] = new TF1("chiSq3", "exp(-x / 2) * sqrt(x / 2 / pi)", 0., 20.); chiSq[3]->SetLineColor(6);
	chiSq[4] = new TF1("chiSq3", "exp(-x / 2) * x / 4", 0., 20.);
	for (int i = 0; i < 5; i++) {
		if (s2Min_x[i]->GetEntries() > 0) { s2Min_x[i]->Scale(1. / s2Min_x[i]->Integral("width")); if (chiSq[i]) AddFunction(s2Min_x[i], chiSq[i]); }
		if (s2Min_y[i]->GetEntries() > 0) { s2Min_y[i]->Scale(1. / s2Min_y[i]->Integral("width")); if (chiSq[i]) AddFunction(s2Min_y[i], chiSq[i]); }
	}

	TGraphErrors* tRes_y = FitAndDrawErrProfile(tYrelErr);
	TGraphErrors* tRes_x = FitAndDrawErrProfile(tXrelErr);
	TGraphErrors* tRes = FitAndDrawErrProfile(tRelErr);
	TGraphErrors* tResLog = FitAndDrawErrProfile(tRelErrLog, true, true);
	TGraphErrors* tResLog_y = FitAndDrawErrProfile(tYRelErrLog, true, true);


	// =================================== print summary
	map<signed int, string> pNames;
	pNames[2212] = "p+";
	pNames[2112] = "n";
	pNames[211] = "pi+";
	pNames[-211] = "pi-";
	pNames[11] = "e-";
	pNames[-11] = "e+";
	pNames[22] = "gamma";

	printf("\n\nSUMMARY:\n");
	printf("------------------------------------------------------------------------------------------------------------------------\n");
	printf("  particle |     total | triggered | OK recon. |   no road | no good r | incon. x* | incon. y* | inc. th_x | inc. th_y |\n");
	printf("------------------------------------------------------------------------------------------------------------------------\n");
	for (map<signed int, EffStruct>::iterator it = counters.begin(); it != counters.end(); ++it) {
		signed int ID = it->first;
		EffStruct& c = it->second;
		printf("%10s |%10i |%10i |%10i |%10i |%10i |%10i |%10i |%10i |%10i |\n", pNames[it->first].c_str(), c.tot, c.trig, c.rOK, c.rNoRoad, 
				c.rNoGoodRoad, c.rBadVertexX, c.rBadVertexY, c.rBadAngleX, c.rBadAngleY);
	}
	printf("------------------------------------------------------------------------------------------------------------------------\n");


	// =================================== draw
	int index = 0;
	TCanvas *can;

	/*
	can = new TCanvas("c1", "Angular results", 600, 600);
	can->Divide(2, 2); index = 0;
	can->cd(1); thetaYcorr->Draw("AP");
	can->cd(2); thetaYerr->Draw("");
	can->cd(3); thetaXcorr->Draw("AP");
	can->cd(4); thetaXerr->Draw("");
	
	can = new TCanvas("c2", "vertex results", 600, 600);
	can->Divide(2, 2);
	can->cd(1); yStCorr->Draw("AP");
	can->cd(2); yStErr->Draw();
	can->cd(3); xStCorr->Draw("AP");
	can->cd(4); xStErr->Draw();

	can = new TCanvas("c3", "t reco results", 900, 300);
	can->Divide(3, 1);
	can->cd(1); tRes_y->Draw("AP");
	can->cd(2); tRes_x->Draw("AP");
	can->cd(3); tRes->Draw("AP");
	*/
	
	/*
	can = new TCanvas("c4", "statistics", 1200, 600);
	can->Divide(4, 2); index = 0;
	can->cd(++index); yNdfHist->Draw("");
	//can->cd(++index); yResHist0->Draw(); yResHist->Fit("gaus", "q"); yResHist0->Draw("same"); yResHist0->SetLineColor(8);
	can->cd(++index); s2Min_y[2]->Draw("");
	can->cd(++index); //s2Min_y[4]->Draw("");
	can->cd(++index); thYchi->Draw(""); thYchi->Fit("gaus", "q");
	can->cd(++index); xNdfHist->Draw("");
	//can->cd(++index); xResHist->Draw(); xResHist->Fit("gaus", "q");
	can->cd(++index); s2Min_x[2]->Draw(""); 
	can->cd(++index); //s2Min_x[4]->Draw(""); 
	can->cd(++index); thXchi->Draw("");  thXchi->Fit("gaus", "q");
	*/

	/*
	can = new TCanvas("c5", "left-right differences", 1000, 1000);
	can->Divide(2, 2); index = 0;
	can->cd(1); dTh_x->Draw(); dTh_x->Fit("gaus", "q");
	can->cd(2); dTh_y->Draw(); dTh_y->Fit("gaus", "q");
	can->cd(3); dVt_x->Draw(); dVt_x->Fit("gaus", "q");
	can->cd(4); dVt_y->Draw(); dVt_y->Fit("gaus", "q");
	*/

	/*
	can = new TCanvas("c6", "physics", 1200, 400);
	can->Divide(3, 1);
	can->cd(1); phiHist_reco->Draw(); phiHist_sim->Draw("same"); phiHist_sim->SetLineColor(8);
	can->cd(2); thetaHist_sim->Draw(); thetaHist_reco->Draw("same"); thetaHist_sim->SetLineColor(8);
	can->cd(3); tHist_sim->Draw(); tHist_reco->Draw("same"); tHist_sim->SetLineColor(8);
	*/

	can = new TCanvas("c10", "", 500, 500);
	can->SetCrosshair(1);
	can->ToggleEventStatus();	
	//bla->Draw("AP");
	//blaHist->Draw();
	can->cd(2); dTh_y->Draw();



	// =================================== save
	string outFileName = argv[inputIndex];
	outFileName.insert(outFileName.rfind('/') + 1, "anal_");
	TFile *of = new TFile(outFileName.c_str(), "recreate");
	gDirectory = of->mkdir("angles");
	thetaYcorr->Write("thetaCorr_y");
	thetaYerr->Write("thetaErr_y");
	thetaXcorr->Write("thetaCorr_x");
	thetaXerr->Write("thetaErr_x");
	thetaErr->Write("thetaErr");

	gDirectory = of->mkdir("vertex");
	yStCorr->Write("vCorr_y");
	yStErr->Write("vErr_y");
	xStCorr->Write("vCorr_x");
	xStErr->Write("vErr_x");

	gDirectory = of->mkdir("t");
	tRes_y->Write("tRes_y");
	tRes_x->Write("tRes_x");
	tRes->Write("tRes");
	tResLog->Write("tResLog");
	tResLog_y->Write("tResLog_y");

	gDirectory = of->mkdir("statistics");
	yResHist0->Write("");
	yNdfHist->Write("");
	thYchi->Write("");
	xResHist->Write("");
	xNdfHist->Write("");
	thXchi->Write("");
	for (int i = 0; i < 5; i++) { s2Min_x[i]->Write(); s2Min_y[i]->Write(); }

	gDirectory = of->mkdir("one RP fit");
	for (map<unsigned int, OneRPFitResult*>::iterator it = oneRPFitRes.begin(); it != oneRPFitRes.end(); ++it)
		it->second->Write();

	gDirectory = of->mkdir("right-left diff");
	dTh_x->Write(); 
	dTh_y->Write();
	dVt_x->Write();
	dVt_y->Write();

	gDirectory = of->mkdir("physics");
	phiHist_reco->Write(); phiHist_sim->Write();
	thetaHist_reco->Write(); thetaHist_sim->Write();
	tHist_reco->Write(); tHist_sim->Write();
	thHistAcc_y->Write("thetaHistAccepted_y");
	thHistFull_y->Write("thetaHistFull_y");
	
	of->Close();

	// =================================== run root
	if (!quitInTheEnd) {
		app->Run();
		delete app;
	}
	return 0;
}

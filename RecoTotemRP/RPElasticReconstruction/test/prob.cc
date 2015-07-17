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

void PrintUsage()
{
	printf("USAGE: runAnalysis <options> <input file.root> <options>\n");
	printf("OPTIONS:\n");
	printf("\t-q\tTo quit ROOT in the end.\n");
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
	for (unsigned int i = 1; i < argc; i++) {
		if (!strncmp(argv[i], "-q", 2)) { quitInTheEnd = true; continue; }

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

	// open file and get event tree
	TFile *theFile = new TFile(argv[inputIndex]);
	if (theFile->IsZombie()) { printf("No such a file `%s'.\n", argv[inputIndex]); return 2; }
	TTree *events = (TTree*)theFile->Get("Events"); 
	if (!events) { printf("No CMSSW data in the file.\n"); return 3; }
	TTree *metaData = (TTree*)theFile->Get("MetaData"); 


	// =================================== prepare branches
	edm::HepMCProduct srcPr; 
	TBranch *srcBr = events->GetBranch(events->GetAlias("source")); 
	srcBr->SetAddress(&srcPr); 

	edm::HepMCProduct srcSmearPr;
	TBranch *srcSmearBr;
	if (events->GetAlias("EnergyVertexSmeared")) {
		srcSmearBr = events->GetBranch(events->GetAlias("EnergyVertexSmeared")); 
		srcSmearBr->SetAddress(&srcSmearPr);
	} else return 1;

	RPRecoElasticEvent elasticReco;
	TBranch *elRecoBr = GetBranchByAlias(events, "ElasticReconstruction", &elasticReco);

	// =================================== prepare histograms and outputs
	TH1::AddDirectory(kFALSE);


	TH1D *thetaYerr = new TH1D("thetaYerr", ";#vartheta_{y}^{reco} - #vartheta_{y}^{sim}   (rad)", 100, -1E-6, 1E-6);
	TH1D *thetaXerr = new TH1D("thetaXerr", ";#vartheta_{x}^{reco} - #vartheta_{x}^{sim}   (rad)", 100, -1E-6, 1E-6);
	TH1D *thetaErr = new TH1D("thetaErr", ";#vartheta^{reco} - #vartheta^{sim}   (rad)", 100, -1E-6, 1E-6);
	TGraph *bla = new TGraph();
	TGraph *bla2 = new TGraph();


	// =================================== loop over events

	for(unsigned int i = 0; i < events->GetEntries(); i++) { 
		if (i == 10) return 0;

		// ------------------------------ GET EVENTS
		//.events->GetEntry(i); //. doesn't work !
		srcBr->GetEntry(i);
		srcSmearBr->GetEntry(i);
		elRecoBr->GetEntry(i);
		HepMC::GenEvent *srcEv = (HepMC::GenEvent *) srcPr.GetEvent(); 
		HepMC::GenEvent *smearedSrcEv =(HepMC::GenEvent *) srcSmearPr.GetEvent();
		if ((i % 1000) == 0)
			printf("[0;32mEVENT[0m	MC=%i,   in order=%i\n", srcEv->event_number(), i);

		// ------------------------------ VERTEX POSITION (only one vertex assumed)
		HepMC::GenEvent::vertex_iterator vI = smearedSrcEv->vertices_begin();
		HepMC::FourVector vertex = (*vI)->position();		// in mm
		//.if (verbose) cout << "vertex = " << v << endl;
		double x_st = vertex.x() * 1E-3;
		double y_st = vertex.y() * 1E-3;

		// ------------------------------ UNSMEARED
		if (srcEv->signal_process_id() != 91) return 2;

		HepMC::FourVector p1, p2;
		p1 = srcEv->barcode_to_particle(3)->momentum();
		double th_y = p1.y() / p1.z();
		double th_x = p1.x() / p1.z();
		double th = sqrt(th_x*th_x + th_y*th_y);
		double phi = atan2(th_y, th_x);

		// ------------------------------ SMEARED
		p1 = smearedSrcEv->barcode_to_particle(3)->momentum();
		p2 = smearedSrcEv->barcode_to_particle(4)->momentum();
		double ths_x = ( p1.x() / p1.z()  +  p2.x() / p2.z() ) / 2.;
		double ths_y = ( p1.y() / p1.z()  +  p2.y() / p2.z() ) / 2.;
		double ths = sqrt(ths_x*ths_x + ths_y*ths_y);
		double phis = atan2(ths_y, ths_x);

		printf("%.2E %.2E\n", p1.x() / p1.z() - th_x, p2.x() / p2.z() - th_x);

		double SpatialUncertaintyCorrection = 1.;
		// ------------------------------ ELASTIC RECONSTRUCTION PLOTS
		if (elasticReco.isValid()) {
			RPRecoElasticEvent::fit_type &f = elasticReco.result;
			//double thp_x = f.th_x, thp_y = f.th_y;
			double thp_x = ths_x, thp_y = ths_y;

			double thp = sqrt(thp_y*thp_y + thp_x*thp_x);
			double phip = atan2(thp_y, thp_x);

			double dth_x = thp_x - th_x;
			double dth_y = thp_y - th_y;
			double dth = thp - th;
			double phid = atan2(dth_y, dth_x);

			thetaYerr->Fill(dth_x);
			thetaXerr->Fill(dth_y);
			thetaErr->Fill(dth);

			bla->SetPoint(bla->GetN(), phi, phip);
			bla2->SetPoint(bla2->GetN(), phi, phid);
		}
	}


	// =================================== draw
	int index = 0;
	TCanvas *can;

	can = new TCanvas("c7", "bla", 1200, 400);
	can->Divide(3, 1);
	can->cd(1); bla->Draw("AP");
	can->cd(2); bla2->Draw("AP");
	can->cd(3);
	thetaXerr->Draw(""); thetaXerr->SetLineColor(1);
	thetaYerr->Draw("same"); thetaYerr->SetLineColor(2);
	thetaErr->Draw("same"); thetaErr->SetLineColor(4);

	// =================================== run root
	if (!quitInTheEnd) {
		app->Run();
		delete app;
	}
	return 0;
}

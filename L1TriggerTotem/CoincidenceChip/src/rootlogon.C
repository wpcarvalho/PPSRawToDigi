#ifndef _L1TriggerTotemrootlogon_C_
#define _L1TriggerTotemrootlogon_C_

// system
#include <iostream>

// ROOT
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
// #include "FWCore/FWLite/interface/AutoLibraryLoader.h"

void rootlogon(){
	using namespace std;

	cout << endl << "Welcome to my rootlogon.C" << endl;
	cout << "For approved plots use: gROOT->SetStyle(\"totem\");";
	cout << endl << endl;

	//
	// gSystem->Load("libFWCoreFWLite.so");
	// AutoLibraryLoader::enable();

	//..totem style from RooLogon.C in workdir
	TStyle *totemStyle= new TStyle("totem","totem approved plots style");

	// use plain black on white colors
	// totemStyle->SetFillColor(0);
	// totemStyle->SetFillColor(kWhite);
	totemStyle->SetFillColor(10);
	totemStyle->SetFrameBorderMode(0);
	totemStyle->SetCanvasBorderMode(0);
	totemStyle->SetPadBorderMode(0);
	totemStyle->SetPadColor(0);
	totemStyle->SetCanvasColor(0);
	totemStyle->SetStatColor(0);
	totemStyle->SetTitleFillColor(0);
	totemStyle->SetHistFillColor(kBlue-7);
	totemStyle->SetFrameFillColor(0);


	// set the paper & margin sizes
	totemStyle->SetPaperSize(20,26);
	//totemStyle->SetPaperSize(TStyle::EPaperSize::kA4);
	// totemStyle->SetPadTopMargin(0.05);
	totemStyle->SetPadRightMargin(0.15);
	// totemStyle->SetPadBottomMargin(0.16);
	// totemStyle->SetPadLeftMargin(0.12);

	/*
	// use large Times-Roman fonts
	totemStyle->SetTextFont(132);
	totemStyle->SetTextSize(0.08);
	totemStyle->SetLabelFont(132,"x");
	totemStyle->SetLabelFont(132,"y");
	totemStyle->SetLabelFont(132,"z");
	totemStyle->SetLabelSize(0.05,"x");
	totemStyle->SetTitleSize(0.06,"x");
	totemStyle->SetLabelSize(0.05,"y");
	totemStyle->SetTitleSize(0.06,"y");
	totemStyle->SetLabelSize(0.05,"z");
	totemStyle->SetTitleSize(0.06,"z");*/

	// use bold lines and markers
	//totemStyle->SetMarkerStyle(20);
	//totemStyle->SetHistLineWidth(1.85);
	//totemStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

	// get rid of X error bars and y error bar caps
	//totemStyle->SetErrorX(0.001);

	// do not display any of the standard histogram decorations
	// totemStyle->SetOptTitle(0);
	totemStyle->SetOptTitle(1);

	// The type of information printed in the histogram statistics box
	//  can be selected via the parameter mode.
	//  The parameter mode can be = ksiourmen  (default = 000001111)
	//    k = 1;  kurtosis printed
	//    k = 2;  kurtosis and kurtosis error printed
	//    s = 1;  skewness printed
	//    s = 2;  skewness and skewness error printed
	//    i = 1;  integral of bins printed
	//    o = 1;  number of overflows printed
	//    u = 1;  number of underflows printed
	//    r = 1;  rms printed
	//    r = 2;  rms and rms error printed
	//    m = 1;  mean value printed
	//    m = 2;  mean and mean error values printed
	//    e = 1;  number of entries printed
	//    n = 1;  name of histogram is printed
	// all options On: totemStyle->SetOptStat(221112211);
	// totemStyle->SetOptStat(000110010);
	totemStyle->SetOptStat(110011);
	//totemStyle->SetOptFit(0);

	// put tick marks on top and RHS of plots
	//totemStyle->SetPadTickX(1);
	//totemStyle->SetPadTickY(1);

	totemStyle->SetPalette(1);
	//totemStyle->SetStatX(0);
	//totemStyle->SetStatY(0);

	//gROOT->SetStyle("Plain");

	gROOT->SetStyle("totem");
	// gROOT->ForceStyle(); 
	//gStyle->SetPadTickX(1);
	//gStyle->SetPadTickY(1);

	//gStyle->SetCanvasColor(0);
	//gStyle->SetFrameBorderMode(0);
	//gStyle->SetStatBorderSize(1);
	//gStyle->SetFrameFillColor(0);
	//gStyle->SetTitleFillColor(0);
	//gStyle->SetPalette(1);

	//gSystem->Load("/usr/lib/libpthread.so");
	//gSystem->Load("$ROOTSYS/lib/libThread.so");
	//gROOT->ProcessLine(".x /afs/cern.ch/user/j/jprochaz/private/config/rootlogon.C");
}

#endif

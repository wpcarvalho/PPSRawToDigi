/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $Id: SmearingValLibrary.cc 9977 2015-01-12 14:00:26Z tsodzawi $
* $Revision: 9977 $
* $Date: 2015-01-12 16:00:26 +0200 (pon, 12 sty 2015) $
*
****************************************************************************/

#include <stdio.h>
#include <iostream>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"

#include "TH1D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TF1.h"

#include "TSystem.h"

#include "TotemRPValidation/BeamSmearing/interface/SmearingValLibrary.h"

using namespace std;
using namespace edm;
using namespace HepMC;

//----------------------------------------------------------------------------------------------------

const char* SmearingValLibrary::regionNames[] = {"forward", "central", "backward"};

SmearingValLibrary::SmearingValLibrary(const edm::ParameterSet& conf) :
  verbosity(conf.getUntrackedParameter<unsigned int>("verbosity", 0)),
  generatorLabel(conf.getParameter<std::string>("generatorLabel")),
  originalLabel(conf.getParameter<std::string>("originalLabel")),
  outputFile(conf.getParameter<std::string>("outputFile")),
  thetaLim(conf.getParameter<double>("thetaLimit"))
{
}

//----------------------------------------------------------------------------------------------------


SmearingValLibrary::~SmearingValLibrary()
{
}

//----------------------------------------------------------------------------------------------------

SmearingValLibrary::SmearInfo::SmearInfo()
{
  xi = new TH1D("xi", "(E_{smeared} - E_{0}) / E_{0}", 500, -15E-4, 15E-4);
  th_x = new TH1D("th_x", "#theta_{x}: smeared - orig", 1000, -50E-6, 150E-6);
  th_y = new TH1D("th_y", "#theta_{y}: smeared - orig", 500, -50E-6, 50E-6);
  t = new TProfile("t", "t: (smeared - orig)/orig", 10, 0., 1.);
}

//----------------------------------------------------------------------------------------------------

void SmearingValLibrary::initialize()
{
  gSystem->Load("libMinuit.so");
  TH1::AddDirectory(kFALSE);

  vtx_x = new TH1D("vtx_x", "vertex X position;#mum", 500, -1000., 1000.);
  vtx_y = new TH1D("vtx_y", "vertex Y position;#mum", 500, -1000., 1000.);
  vtx_z = new TH1D("vtx_z", "vertex Z position;cm", 500, -20., 20.); 
}

//----------------------------------------------------------------------------------------------------

void SmearingValLibrary::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  edm::Handle< HepMCProduct > mcPr;

  event.getByLabel(generatorLabel, mcPr);
  HepMC::GenEvent *smear = (HepMC::GenEvent *) mcPr->GetEvent();

  event.getByLabel("SmearingGenerator", originalLabel, mcPr);
  HepMC::GenEvent *orig = (HepMC::GenEvent *) mcPr->GetEvent();

  // vertices
  GenEvent::vertex_const_iterator vold, vnew;
  for (vold = orig->vertices_begin(), vnew = smear->vertices_begin(); 
      vold != orig->vertices_end() && vnew != smear->vertices_end(); ++vold, ++vnew) {
    if (verbosity > 5) {
	  const FourVector &vo = (*vold)->position();
	  const FourVector &vn = (*vnew)->position();
	  cout << "Vertex\n\told: [" << vo.x() << ", " << vo.y() << ", " << vo.z() << ", " << vo.t() << "]\n"
		   << "\tnew: [" << vn.x() << ", " << vn.y() << ", " << vn.z() << ", " << vn.t() << "]\n";
    }
    
    vtx_x->Fill((*vnew)->position().x() * 1E3);    // from mm to um
    vtx_y->Fill((*vnew)->position().y() * 1E3);
    vtx_z->Fill((*vnew)->position().z() * 1E-1);  // from mm to cm
  }

  // particles
  GenEvent::particle_const_iterator pold, pnew;
  for (pold = orig->particles_begin(), pnew = smear->particles_begin();
      pold != orig->particles_end() && pnew != smear->particles_end(); ++pold, ++pnew) {
    FourVector o = (*pold)->momentum(), n = (*pnew)->momentum();
  
    // determine direction region
    double th = o.theta();
    region idx = rCent;
    if (th < thetaLim) idx = rForw;
    if (th > (M_PI - thetaLim)) idx = rBack;
  
      if (verbosity > 5)
        cout << "particle\n\told: [" << o.x() << ", " << o.y() << ", " << o.z() << ", " << o.t()
        << "]\n\tnew: [" << n.x() << ", " << n.y() << ", " << n.z() << ", " << n.t()
        << "]\n\tregion: " << idx << endl;
  
    // get pointer to histograms
    hists[idx].xi->Fill((n.e() - o.e()) / o.e());
    double othx = o.x() / o.rho(), othy = o.y() / o.rho();
    double nthx = n.x() / n.rho(), nthy = n.y() / n.rho();
    hists[idx].th_x->Fill(nthx - othx);
    hists[idx].th_y->Fill(nthy - othy);
    double ot = othx*othx + othy*othy;
    double nt = nthx*nthx + nthy*nthy;
    const double p2 = 7E3*7E3;
    if (ot != 0.)
      hists[idx].t->Fill(ot * p2, (nt - ot) / ot);
  }
}

//----------------------------------------------------------------------------------------------------

void SmearingValLibrary::finalize()
{
  if (verbosity) printf(">> SmearingValLibrary::finalize\n");
  TF1 ff("ff", "[0]/sqrt(x)");
  for (map<region, SmearInfo>::iterator it = hists.begin(); it != hists.end(); ++it) {
    it->second.t->Fit(&ff, "qw");
    if (verbosity)
      printf("\t%s: A = %E\n", regionNames[it->first], ff.GetParameter(0));
  }
}

//----------------------------------------------------------------------------------------------------

void SmearingValLibrary::writeHistogramsToFile()
{

  TFile *of = TFile::Open(outputFile.c_str(), "recreate");
  
  if(!of || !of->IsWritable()){
      std::cout<<"Output file not opened correctly!!"<<std::endl;
  }    

  gDirectory = of->mkdir("vertex");
  vtx_x->Write();
  vtx_y->Write();
  vtx_z->Write();

  for (map<region, SmearInfo>::iterator it = hists.begin(); it != hists.end(); ++it) {
    gDirectory = of->mkdir(regionNames[it->first]);
    it->second.xi->Write();
    it->second.th_x->Write();
    it->second.th_y->Write();
    it->second.t->Write();
  }

  of->Close();
}


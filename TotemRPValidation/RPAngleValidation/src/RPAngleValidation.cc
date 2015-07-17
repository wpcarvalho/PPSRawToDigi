#include <cmath>
#include "TH2D.h"
#include "TFile.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimG4Core/Notification/interface/SimG4Exception.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "TotemRPValidation/RPAngleValidation/interface/RPAngleValidation.h"

/*
 * TODO:
 * thetax, thetay wygenerowanych protonów
 * te same kąty ze zrekonstruowanej ścieżki z RP
 * dla różnych ksi (0-1%, 1-10%, 10-30%, >30%)
 * wykres theta_RP od theta_orig
 *
 * Algorithm:
 * 1. beginJob: get beam optics params, save theta_orig
 * 2. get two reconstructed hit positions
 * 3. calculate line equation
 * 4. calculate theta_RP
 * 5. store in histogram
 */

RPAngleValidation::RPAngleValidation(const edm::ParameterSet& conf)
 : conf_(conf)
{
	  OriginalHepMCProductLabel_ = conf.getParameter<std::string>("OriginalHepMCProductLabel");
	  OriginalHepMCModuleName_ = conf.getParameter<std::string>("OriginalHepMCModuleName");
	  SmearedHepMCProductLabel_ = conf.getParameter<std::string>("SmearedHepMCProductLabel");
	  SmearedHepMCModuleName_ = conf.getParameter<std::string>("SmearedHepMCModuleName");
	  verbosity_ = conf.getParameter<int>("Verbosity");
	  hist_file_name_ = conf.getParameter<std::string>("HistogramFileName");
	  rpFittedTrackCollectionLabel = conf.getParameter<edm::InputTag>("RPFittedTrackCollectionLabel");
	  //init rp ids for the stations
	  //order: from ip to outside, then top bottom
	  //top first vertical, bottom first vertical, first horizontal, second horizontal, last top vertical, last bottom vertical
	  for(int i=0; i<NPOTS; ++i)
	    station_ids_[rp_150_l].push_back(i);
	  for(int i=20; i<20+NPOTS; ++i)
	    station_ids_[rp_220_l].push_back(i);
	  for(int i=100; i<100+NPOTS; ++i)
	    station_ids_[rp_150_r].push_back(i);
	  for(int i=120; i<120+NPOTS; ++i)
	    station_ids_[rp_220_r].push_back(i);
}


RPAngleValidation::~RPAngleValidation()
{
}

std::string RPAngleValidation::StIdToName(int st_id)
{
	std::string name;
	switch (st_id) {
		case rp_150_l: name="rp_150_l"; break;
		case rp_150_r: name="rp_150_r"; break;
		case rp_220_l: name="rp_220_l"; break;
		case rp_220_r: name="rp_220_r"; break;
	}
	return name;
}

void RPAngleValidation::beginRun(edm::Run const&, edm::EventSetup const& es)
{
	  edm::ESHandle<BeamOpticsParams> BOParH;
	  es.get<BeamOpticsParamsRcd>().get(BOParH);
	  if(!BOParH.isValid())
		throw cms::Exception("RPAngleValidation") << " edm::ESHandle<BeamOpticsParams> is invalid";
	  BOPar_ = *BOParH;
	  for (int st_id=0;st_id<NSTATIONS;st_id++) {
		  std::string name="angle_dist_ip_vs_rp_theta_x_"+StIdToName(st_id);
		  angle_dists_[0][st_id]=std::auto_ptr<TH2D>(new TH2D(name.c_str(), name.c_str(), 1000, -1.0e-3, 1.0e-3, 1000, -1.0e-3, 1.0e-3));
		  name="angle_dist_ip_vs_rp_theta_y_"+StIdToName(st_id);
		  angle_dists_[1][st_id]=std::auto_ptr<TH2D>(new TH2D(name.c_str(), name.c_str(), 1000, -1.0e-3, 1.0e-3, 1000, -1.0e-3, 1.0e-3));
		  name="1d_angle_dist_ip_theta_x_"+StIdToName(st_id);
		  angle_dists_1d_[0][0][st_id]=std::auto_ptr<TH1D>(new TH1D(name.c_str(), name.c_str(), 1000, -5.0e-3, 5.0e-3));
		  name="1d_angle_dist_rp_theta_x_"+StIdToName(st_id);
		  angle_dists_1d_[1][0][st_id]=std::auto_ptr<TH1D>(new TH1D(name.c_str(), name.c_str(), 1000, -5.0e-3, 5.0e-3));
		  name="1d_angle_dist_ip_theta_y_"+StIdToName(st_id);
		  angle_dists_1d_[0][1][st_id]=std::auto_ptr<TH1D>(new TH1D(name.c_str(), name.c_str(), 1000, -5.0e-3, 5.0e-3));
		  name="1d_angle_dist_rp_theta_y_"+StIdToName(st_id);
		  angle_dists_1d_[1][1][st_id]=std::auto_ptr<TH1D>(new TH1D(name.c_str(), name.c_str(), 1000, -5.0e-3, 5.0e-3));
	  }
}

bool RPAngleValidation::FindProtons(const edm::Handle<edm::HepMCProduct> &HepMCEvt, int smeared)
{
  if(!HepMCEvt.isValid())
     throw SimG4Exception("RPReconstructedTracksValidation : Unable to find HepMCProduct(HepMC::GenEvent) in edm::Event  ");

  const HepMC::GenEvent *evt = HepMCEvt->GetEvent();

  int count[2] = {0,0};
  prot_[0][smeared].found = false;
  prot_[1][smeared].found = false;

  for(HepMC::GenEvent::particle_const_iterator it = evt->particles_begin();
      it != evt->particles_end();
      ++it )
  {
    HepMC::GenParticle * g = (*it);
    int g_status = g->status();
    int pdg_id = g->pdg_id();

    if(verbosity_ > 1)
    {
      g->print(std::cout);
// this is crashing for some reason
//      std::cout<<g->production_vertex()->position()<<std::endl;
    }

    // scanning only for undecayed particles, status == 1, see http://cepa.fnal.gov/psm/simulation/mcgen/lund/pythia_manual/pythia6.3/pythia6301/node39.html
    // protons have pdg_id == 2212, see http://www-pdg.lbl.gov/2009/reviews/rpp2009-rev-monte-carlo-numbering.pdf
    if (g_status == 1 && pdg_id == 2212)
    {
      const HepMC::FourVector &vtx = g->production_vertex()->position();
      const HepMC::FourVector &mom  = g->momentum();
      if(verbosity_ > 1) {
        std::cout<<"Setting the ";
        smeared?std::cout<<"smeared":std::cout<<"primary";
        std::cout<< " protons."<<std::endl;
      }
      if(mom.z())
      {
    	int arm=(mom.z()>0);
        ++count[arm];
        prot_[arm][smeared].vertex = vtx; //[mm]
        prot_[arm][smeared].momentum = mom;  //[GeV]
        prot_[arm][smeared].found = true;
        prot_[arm][smeared].thetaX = BOPar_.ComputeCrossingAngleCorrectedThetaX(mom);
        prot_[arm][smeared].thetaY = BOPar_.ComputeCrossingAngleCorrectedThetaY(mom);
      }
    }
  } // end loop on HepMC particles

  bool result = count[1]<=1 && count[0]<=1 && (count[1]>0 || count[0]>0);

  if(verbosity_ > 1)
    std::cout<<"right_count="<<count[1]<<" left_count="<<count[0]
        <<" result="<<result<<std::endl;
  return result;
}

//return number of pots through which the proton goes
int RPAngleValidation::SelectHits(const RPFittedTrackCollection &tracks,
    rec_tracks_collection & coll, const station_rp_ids_type& st_ids)
{
  RPFittedTrackCollection::const_iterator it;
  RPFittedTrackCollection::const_iterator low_it;
  RPFittedTrackCollection::const_iterator high_it;

  assert(st_ids.size() == 6);

  low_it = tracks.lower_bound(st_ids[0]);
  high_it = tracks.upper_bound(st_ids[5]);

  for(it = low_it; it!=high_it; ++it)
  {
    if(!it->second.IsValid())
      continue;
    TrackHit tmp;
    tmp.x=it->second.X0();
    tmp.y=it->second.Y0();
    tmp.z=it->second.Z0();
    coll[it->first]=tmp;
  }

  rec_tracks_collection::const_iterator end = coll.end();

  bool station_accepted =
    (coll.find(st_ids[near_top])!=end || coll.find(st_ids[near_bottom])!=end || coll.find(st_ids[near_horiz])!=end) &&
    (coll.find(st_ids[far_horiz])!=end || coll.find(st_ids[far_top])!=end || coll.find(st_ids[far_bottom])!=end) &&
    !(coll.find(st_ids[near_top])!=end && coll.find(st_ids[near_bottom])!=end) &&
    !(coll.find(st_ids[far_top])!=end && coll.find(st_ids[far_bottom])!=end);

  if(!station_accepted)
    coll.clear();

  return coll.size();
}

void RPAngleValidation::analyze(const edm::Event& e, const edm::EventSetup& es)
{
	  edm::Handle< RPFittedTrackCollection > input;
	  e.getByLabel(rpFittedTrackCollectionLabel, input);
	  if (!input.isValid())
	    throw cms::Exception("RPAngleValidation") << "edm::Handle< RPFittedTrackCollection > is invalid";

	  edm::Handle<edm::HepMCProduct> OriginalHepMCEvt;
	  e.getByLabel(OriginalHepMCModuleName_, OriginalHepMCProductLabel_, OriginalHepMCEvt );
	  if (!OriginalHepMCEvt.isValid())
	    throw cms::Exception("RPAngleValidation") << "edm::Handle<edm::HepMCProduct> OriginalHepMCEvt is invalid";

	  edm::Handle<edm::HepMCProduct> SmearedHepMCEvt;
	  e.getByLabel(SmearedHepMCModuleName_, SmearedHepMCProductLabel_, SmearedHepMCEvt );
	  if (!SmearedHepMCEvt.isValid())
	    throw cms::Exception("RPAngleValidation") << "edm::Handle<edm::HepMCProduct> SmearedHepMCEvt is invalid";

	  if(!FindProtons(OriginalHepMCEvt,0))
	  {
	    if(verbosity_ > 2)
	      std::cout<<"Primary protons not found!! Skipping the event."<<std::endl;
	    return;
	  }

	  if(!FindProtons(SmearedHepMCEvt,1))
	  {
	    if(verbosity_)
	      std::cout<<"Smeared protons not found!! Skipping the event."<<std::endl;
	    return;
	  }


	  for (int i=0;i<NSTATIONS;i++) {
		  rec_tracks_collection hits;
		  if (prot_[i/2][0].found && prot_[i/2][1].found && SelectHits(*input, hits, station_ids_[i])>1) {
			  std::vector<double> trackTheta;
			  if (calcTrackTheta(trackTheta, hits,i/2)<2) {
				  if (verbosity_) std::cout << "RPAngleValidation::analyze: calcTrackTheta returned no value." << std::endl;
				  continue;
			  }
			  double xi=BOPar_.IPSmearedProtonMomentumToXi(prot_[i%2][1].momentum.px(),prot_[i%2][1].momentum.py(),prot_[i%2][1].momentum.pz());
			  angle_dists_[0][i]->Fill(prot_[i/2][1].thetaX,trackTheta[0]);
			  angle_dists_[1][i]->Fill(prot_[i/2][1].thetaY,trackTheta[1]);
			  if (verbosity_) {
				  std::cout<<"RPAngleValidation::analyze: thetaX(IP)   thetaX(RP)   thetaY(IP)  thetaY(RP)   xi" << std::endl;
			      std::cout<<"                            " << prot_[i%2][1].thetaX << " " << trackTheta[0] << " " <<
			                                                   prot_[i%2][1].thetaY << " " << trackTheta[1] << " " <<
			                                                   xi << std::endl;
			  }
			  angle_dists_1d_[0][0][i]->Fill(prot_[i/2][1].thetaX);
			  angle_dists_1d_[0][1][i]->Fill(prot_[i/2][1].thetaY);
			  angle_dists_1d_[1][0][i]->Fill(trackTheta[0]);
			  angle_dists_1d_[1][1][i]->Fill(trackTheta[1]);
		  }
	  }
}

// calculates tangent of the angle of the fitted track in the ZX or ZY plane
int RPAngleValidation::calcTrackTheta(std::vector<double> &theta, const rec_tracks_collection &coll, int arm)
{
	rec_tracks_collection::const_iterator it;
	rec_tracks_collection vert, horiz;
	TrackHit th1, th2;
	double dz;
	if (verbosity_) {
		std::cout<<"RPAngleValidation::calcTrackTheta: arm "<< (arm?"right":"left")
		         << ", got "<<coll.size() << " tracks, RPs ";
		for (rec_tracks_collection::const_iterator i=coll.begin(); i!=coll.end(); i++) {
			std::cout << " " << i->first;
		}
		std::cout << std::endl;
	}
	if (coll.size()<2)
		return 0;

	for (it=coll.begin(); it!=coll.end(); it++) {
		int pot_id = it->first % 20;
		switch (pot_id) {
			case near_top:
			case near_bottom:
			case far_top:
			case far_bottom: vert[it->first]=it->second; break;
			case near_horiz:
			case far_horiz: horiz[it->first]=it->second; break;
		}
	}
	if (verbosity_) std::cout << "RPAngleValidation::calcTrackTheta: "
	                          << vert.size() << " vertical, "
	                          << horiz.size() << " horizontal pots... ";
	if (vert.size()>=2) {
		it=vert.begin();
		if (verbosity_) std::cout << "using vertical." << std::endl;
	} else if (horiz.size()>=2) {
		it=horiz.begin();
		if (verbosity_) std::cout << "using horizontal." << std::endl;
	} else
		return 0;
	th1=it->second;
	it++;
	th2=it->second;
	if (verbosity_) std::cout<<"RPAngleValidation::calcTrackTheta: hit1: "
	                         << th1.x << " " << th1.y << " " << th1.z << " hit2: "
	                         << th2.x << " " << th2.y << " " << th2.z << std::endl;
	dz=th2.z-th1.z;
	theta.reserve(2);
	theta[0]=(th2.x-th1.x)/dz;
	theta[1]=(th2.y-th1.y)/dz;
	if (verbosity_) std::cout << "RPAngleValidation::calcTrackTheta: thetaX: " << theta[0]
	                          << " thetaY: " << theta[1] << std::endl;
	return 2;
}

void RPAngleValidation::endJob()
{
	TFile *f = TFile::Open(hist_file_name_.c_str(), "recreate");
	if(f && f->IsWritable()) {
		if(verbosity_) std::cout<<"Writing histograms..."<<std::endl;
		  for (int st_id=0;st_id<NSTATIONS;st_id++) {
			  for (int arm=0; arm<2; arm++) {
				  angle_dists_[arm][st_id]->Write();
				  angle_dists_1d_[0][arm][st_id]->Write();
				  angle_dists_1d_[1][arm][st_id]->Write();
			  }
		  }
		if(verbosity_) std::cout<<"Writing histograms finished."<<std::endl;
	} else {
	    std::cout<<"Output file not opened correctly!!"<<std::endl;
	}
	f->Close();
}

DEFINE_FWK_MODULE(RPAngleValidation);

#include <TFile.h>

#include "MyEvtId.h"      // UA tree event desctiption
#include "MyL1Trig.h"  	 // UA tree event desctiption (arrays)
#include "MyL1TrigOld.h" // UA tree event desctiption (arrays)
#include "MyHLTrig.h"
#include "MyTracks.h"
#include "MyVertex.h"
#include "MyPFJet.h"
#include "MyCaloJet.h"
#include "MyMuon.h"
#include "MyElectron.h"
#include "MyPFCand.h"
#include "MyCaloTower.h"
#include "MyMet.h"
#include "MyFSCHit.h"
#include "MyFSCDigi.h"
#include "MyGenPart.h"

#include "MyCastorRecHit.h"

#include "utilities.h"
#include <map>
#include <limits>
#include <time.h>
#include <stdio.h>

using namespace std;

TBranch *trigDataBranch;

//it would be nice to move next three procedures (checkAndGetBranch, prepareRomanPotBranches, prepareTotemBranches)
//to a place common for both merge12 and combine programs but there are some memory problems that I
//encountered and currently I don't have time to investigate it

/**
 * This function checks if branch of a given name exists.
 * if branch exists - the pointer is returned
 * if not - the program is terminated with error message
 */
TBranch* checkAndGetBranch(TTree* tree, string branchName) {
	TBranch *branch = tree->GetBranch(branchName.c_str());
	if (!branch) {
		string dotBranchName = branchName + ".";
		branch = tree->GetBranch(dotBranchName.c_str());
	}
	if (!branch) {
		tree->Print();
		error(" No data branch " + branchName + " found in input file!");
	}
	return branch;
}

/**
 * This function prepares roman pot branches for merging.
 * First it connects roman pot variables with proper branches in totem ntuple.
 * Then it creates appropriate branches in output ntuple (with this variables).
 */
void prepareRomanPotBranches(TTree *outputTree, TTree *tree_totem) {
	checkAndGetBranch(tree_totem, "rec_prot_left")->SetAddress(&recProtLeft);
	checkAndGetBranch(tree_totem, "rec_prot_right")->SetAddress(&recProtRight);
	checkAndGetBranch(tree_totem, "rec_prot_pair")->SetAddress(&recProtonPair);

	outputTree->Branch("rec_prot_left.", &recProtLeft);
	outputTree->Branch("rec_prot_right.", &recProtRight);
	outputTree->Branch("rec_prot_pair.", &recProtonPair);

	char br_name[200];
	for (unsigned int a = 0; a < 2; ++a) { //the digis, patterns and tracks are stored separately
										   //for every RP, this loops iterate over RP ids
		int s = 2;
		for (unsigned int r = 0; r < 6; r++) {
			unsigned int id = 100 * a + 10 * s + r;
			if (includeDigi) {
				sprintf(br_name, "digi_rp_%u.", id);
				checkAndGetBranch(tree_totem, br_name)->SetAddress(
						&digi_info_[id]);
				outputTree->Branch(br_name, &digi_info_[id]);
			}
			if (includePatterns) {
				sprintf(br_name, "par_patterns_rp_%u.", id);
				checkAndGetBranch(tree_totem, br_name)->SetAddress(
						&par_patterns_info_[id]);
				outputTree->Branch(br_name, &par_patterns_info_[id]);

				sprintf(br_name, "nonpar_patterns_rp_%u.", id);
				checkAndGetBranch(tree_totem, br_name)->SetAddress(
						&nonpar_patterns_info_[id]);
				outputTree->Branch(br_name, &nonpar_patterns_info_[id]);
			}

			sprintf(br_name, "track_rp_%u.", id);
			checkAndGetBranch(tree_totem, br_name)->SetAddress(
					&track_info_[id]);
			outputTree->Branch(br_name, &track_info_[id]);
			if (includeMultiTracks) {
				sprintf(br_name, "multi_track_rp_%u", id);
				checkAndGetBranch(tree_totem, br_name)->SetAddress(
						&multi_track_info_[id]);
				outputTree->Branch(br_name, &multi_track_info_[id]);
			}
		}
	}
}

/**
 * This function prepares totem branches in output tree.
 */
void prepareTotemBranches(const string& totemFileName, TTree *outputTree,
		TTree *tree_totem) {
	checkAndGetBranch(tree_totem, "branchT2EV")->SetAddress(&t2Event);
	checkAndGetBranch(tree_totem, "branchT1EV")->SetAddress(&t1Event);
	trigDataBranch = checkAndGetBranch(tree_totem, "trigger_data");
	trigDataBranch->SetAddress(&trigData);
	checkAndGetBranch(tree_totem, "event_info")->SetAddress(&totemEventData);
	prepareRomanPotBranches(outputTree, tree_totem);
	outputTree->Branch("trigger_data.", &trigData);
	outputTree->Branch("event_info.", &totemEventData);
	outputTree->Branch("branchT1EV.", &t1Event);
	outputTree->Branch("branchT2EV.", &t2Event);
}

/**
 * This function merges the totem and cms ntuples. The result is saved in outputFileName.
 */
void mergeNtuples(const int dOrbit, const string& totemFileName,
		const string& cmsFileName, const string &outputFileName) {
	cout << "Opening Totem File: " << totemFileName << endl; //opening input ntuples
	TFile* totemFile = TFile::Open(totemFileName.c_str());
	if (!totemFile || totemFile->IsZombie())
		error("Error opening file");
	TTree* totemTree = (TTree*) totemFile->Get("TotemNtuple");
	if (!totemTree)
		error("No data tree found");
	cout << "Entries: " << totemTree->GetEntries() << endl;
	cout << "Opening CMS File: " << cmsFileName << endl;
	TFile* cmsFile = TFile::Open(cmsFileName.c_str());
	if (!cmsFile || cmsFile->IsZombie())
		error("Error opening file");
	TTree * cmsTree = (TTree*) cmsFile->Get("evt");
	unsigned int cmsSize = cmsTree->GetEntries();
	cout << "CMS entries: " << cmsTree->GetEntries() << endl;
	if(totemTree->GetEntries()!=cmsTree->GetEntries()){
		std::cerr << "Mismatched numbers of events!" << std::endl;
		return;
	}

	MyEvtId *evtcmsUA = 0; //assigns cms variables to input ntuple
	TBranch *evtId = checkAndGetBranch(cmsTree, "evtId");
	evtId->SetAddress(&evtcmsUA);
	MyL1TrigOld * trigcmsUA = 0;
	checkAndGetBranch(cmsTree, "L1TrigOld")->SetAddress(&trigcmsUA);
	//MyHLTrig * hltTrig_cmsUA = 0;
	//checkAndGetBranch(cmsTree, "HLTrig")->SetAddress(&hltTrig_cmsUA);
        // Tracks, vertices
	vector<MyTracks> *tracks_cmsUA = 0;
	checkAndGetBranch(cmsTree, "generalTracks")->SetAddress(&tracks_cmsUA);
	vector<MyVertex> *vertices_cmsUA = 0;
	checkAndGetBranch(cmsTree, "offlinePrimaryVertices")->SetAddress(&vertices_cmsUA);
            // Jets
	/*
	vector<MyPFJet> *ak5PFJets_cmsUA = 0;
	checkAndGetBranch(cmsTree, "ak5PFJets")->SetAddress(&ak5PFJets_cmsUA);
	vector<MyPFJet> *ak7PFJets_cmsUA = 0;
	checkAndGetBranch(cmsTree, "ak7PFJets")->SetAddress(&ak7PFJets_cmsUA);
	vector<MyCaloJet> *ak5CaloJets_cmsUA = 0;
	checkAndGetBranch(cmsTree, "ak5CaloJets")->SetAddress(&ak5CaloJets_cmsUA);
	vector<MyCaloJet> *ak7CaloJets_cmsUA = 0;
	checkAndGetBranch(cmsTree, "ak7CaloJets")->SetAddress(&ak7CaloJets_cmsUA);
        // Muons, electrons 
	vector<MyMuon> *muons_cmsUA = 0;
	checkAndGetBranch(cmsTree, "muons")->SetAddress(&muons_cmsUA);
	vector<MyElectron> *electrons_cmsUA = 0;
	checkAndGetBranch(cmsTree, "gsfElectrons")->SetAddress(&electrons_cmsUA);
        // Particle-flow, calo towers
	vector<MyPFCand> *pFlow_cmsUA = 0;
	checkAndGetBranch(cmsTree, "particleFlow")->SetAddress(&pFlow_cmsUA);
	vector<MyCaloTower> *caloTowers_cmsUA = 0;
	checkAndGetBranch(cmsTree, "caloTowers")->SetAddress(&caloTowers_cmsUA);
        // MET
	vector<MyMet> *met_cmsUA = 0;
	checkAndGetBranch(cmsTree, "met")->SetAddress(&met_cmsUA);
	vector<MyMet> *pfMet_cmsUA = 0;
	checkAndGetBranch(cmsTree, "pfMet")->SetAddress(&pfMet_cmsUA);
	vector<MyMet> *tcMet_cmsUA = 0;
	checkAndGetBranch(cmsTree, "tcMet")->SetAddress(&tcMet_cmsUA);
	*/
        // FSC
	//MyFSCInfo *fscInfo_cmsUA = 0;
	//checkAndGetBranch(cmsTree, "fscInfo")->SetAddress(&fscInfo_cmsUA);
	//vector<MyFSCHit> *fscHits_cmsUA = 0;
	//checkAndGetBranch(cmsTree, "fscHits")->SetAddress(&fscHits_cmsUA);
	//vector<MyFSCDigi> *fscDigis_cmsUA = 0;
	//checkAndGetBranch(cmsTree, "fscDigis")->SetAddress(&fscDigis_cmsUA);
		// Generator info
	vector<MyGenPart> *genPart_cmsUA = 0;
	checkAndGetBranch(cmsTree,"genPart")->SetAddress(&genPart_cmsUA);
	//Merijn: add the castor rechits
	vector<MyCastorRecHit> *castorRecHits_cmsUA = 0;
	checkAndGetBranch(cmsTree,"castorRecHits")->SetAddress(&castorRecHits_cmsUA);

	TFile* output = TFile::Open(outputFileName.c_str(), "RECREATE"); //output ntuple
	output->cd();
	TTree* outputTree = new TTree("cms_totem", "merged");
	outputTree->Branch("cmsEvtUA", &evtcmsUA); //preparing branches in output ntuple
	outputTree->Branch("cmsTrigUA", &trigcmsUA);
	//outputTree->Branch("cmsHLTTrigUA", &hltTrig_cmsUA);
        // Tracks, vertices
	outputTree->Branch("cmsTracksUA", &tracks_cmsUA);
	outputTree->Branch("cmsVerticesUA", &vertices_cmsUA);
        // Jets
	/*
	outputTree->Branch("cmsak5PFJetsUA", &ak5PFJets_cmsUA);
	outputTree->Branch("cmsak7PFJetsUA", &ak7PFJets_cmsUA);
	outputTree->Branch("cmsak5CaloJetsUA", &ak5CaloJets_cmsUA);
	outputTree->Branch("cmsak7CaloJetsUA", &ak7CaloJets_cmsUA);
        // Muons, electrons 
	outputTree->Branch("cmsMuonsUA", &muons_cmsUA);
	outputTree->Branch("cmsElectronsUA", &electrons_cmsUA);
        // Particle-flow, calo towers
	outputTree->Branch("cmsParticleFlowUA", &pFlow_cmsUA);
	outputTree->Branch("cmsCaloTowersUA", &caloTowers_cmsUA);
        // MET
	outputTree->Branch("cmsMETUA", &met_cmsUA);
	outputTree->Branch("cmsPFMETUA", &pfMet_cmsUA);
	outputTree->Branch("cmsTCMETUA", &tcMet_cmsUA);
	*/
        // FSC
	//outputTree->Branch("cmsFSCInfoUA", &fscInfo_cmsUA);
	//outputTree->Branch("cmsFSCHitsUA", &fscHits_cmsUA);
	//outputTree->Branch("cmsFSCDigisUA", &fscDigis_cmsUA);
	
		// Generator info
	outputTree->Branch("cmsGenPart", &genPart_cmsUA);
	//castor rechits:
	outputTree->Branch("castorRecHits", &castorRecHits_cmsUA);

        // TOTEM
	prepareTotemBranches(totemFileName, outputTree, totemTree);
/*
	vector<info> totem; //this vector contains only basic information about totem events
						//it will speed up merging procedure
	map<unsigned int, vector<int> > totemOrbits; //this map helps in searching for totem event with given orbit number
	//key - orbit number, value - vector of indexes of totem events that have this orbit number
	map<unsigned int, vector<int> >::iterator it;
*/
	for (unsigned int cms_i = 0; cms_i < cmsSize; cms_i++) { //iterate over events
		if(cms_i % 10000 == 0) std::cout << "Event: " << cms_i << std::endl;
		totemTree->GetEntry(cms_i);
		cmsTree->GetEntry(cms_i);
		if(evtcmsUA->Evt!=totemEventData->event_no){
			std::cerr << "Mismatched event numbers" << std::endl;
			return;
		}
		outputTree->Fill(); //fill the output ntuple
	} // end loop events

	outputTree->Write();
	output->Close();
}

/**
 * You can use this program by typing:
 * ./mergeNtuples totemNtuplePath cmsNtuplePath outputNtuple
 */
int main(int argc, char** argv) {
	if (argc < 4)
		error("Please specify TOTEM and CMS data files and output file.");
	const string totem(argv[1]);
	const string cms(argv[2]);
	const string outputFileName(argv[3]);
	mergeNtuples(0, totem, cms, outputFileName);
	return 0;
}

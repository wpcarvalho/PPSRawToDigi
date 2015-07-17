#include <TFile.h>
#include <TH1I.h>
#include <algorithm>
#include "MyEvtId.h"       // UA tree event desctiption
#include "utilities.h"
#include "MyL1TrigOld.h"

using namespace std;

const int maxDeltaBX = 10;
const int maxDeltaOrbit = 200;
unsigned int countTOTwrongBX = 0;
const bool lostEvtChck = false; // print-outs

/**
 * The basic cms info is stored in this structure. This structure is used to find the orbit offset.
 */
typedef struct {
	MyEvtId evt;
	MyL1TrigOld trigg;
} CmsInfo;

/**
 * This function compares two cms events. First, it compares orbit numbers and then (if equal)
 * bunch numbers.
 */
bool sortEvt(CmsInfo info1, CmsInfo info2) {
	if (info1.evt.Orbit != info2.evt.Orbit) {
		return (info1.evt.Orbit < info2.evt.Orbit);
	} else {
		return (info1.evt.Bunch < info2.evt.Bunch);
	}
	return true; //this is just to avoid compilator warnings
}

/**
 * This function returns sorted vector of cms events. The vector includes just that cms events that
 * correspond to totem events (according to orbit number + some margin).
 * The sorting is necessary before finding the orbit offset.
 *
 * totemOrbitMin - orbit number of the first totem event that will be used for offset searching
 * totemOrbitMax - orbit number of the last totem event that will be used for offset searching
 */
vector<CmsInfo> getSortedCmsInfo(TTree *tree_cms, unsigned int totemOrbitMin,
		unsigned int totemOrbitMax) {
	const unsigned int orbitMargin = 2000; //just for safety we take some more events
	unsigned int n_cms = tree_cms->GetEntries();
	MyEvtId *evt_cms = 0;
	MyL1TrigOld *trigg_cms = 0;
	tree_cms->SetBranchAddress("evtId", &evt_cms);
	tree_cms->SetBranchAddress("L1TrigOld", &trigg_cms);

	CmsInfo info_cms;
	vector<CmsInfo> vectorCmsInfo;
	int min = totemOrbitMin - orbitMargin;
	int max = totemOrbitMax + orbitMargin;

	for (unsigned int i = 0; i < n_cms; i++) { //choose cms events within given orbit range
		tree_cms->GetEntry(i);
		if (evt_cms->Orbit > min && evt_cms->Orbit < max) {
			info_cms.evt = *(evt_cms);
			info_cms.trigg = *(trigg_cms);
			vectorCmsInfo.push_back(info_cms);
		}
	}
	sort(vectorCmsInfo.begin(), vectorCmsInfo.end(), sortEvt);
	return vectorCmsInfo;
}

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
// connect totem tree, branches
	checkAndGetBranch(tree_totem, "branchT2EV")->SetAddress(&t2Event);
	checkAndGetBranch(tree_totem, "branchT1EV")->SetAddress(&t1Event);
	checkAndGetBranch(tree_totem, "trigger_data")->SetAddress(&trigData);
	checkAndGetBranch(tree_totem, "event_info")->SetAddress(&totemEventData);
	prepareRomanPotBranches(outputTree, tree_totem);
	outputTree->Branch("trigger_data.", &trigData);
	outputTree->Branch("event_info.", &totemEventData);
	outputTree->Branch("branchT1EV.", &t1Event);
	outputTree->Branch("branchT2EV.", &t2Event);
}

/**
 * The aim of this function is to find and return orbit offset. (cms.orbit - totem.orbit).
 * The cmsInfo vector must be sorted.
 * cmsInfo - vector of sorted cms events (these events correspond to totem events (orbit range)
 * totem - vector of totem events (not all events just a sample of a given number (you can set variable in the code "const int n_totem")
 */
unsigned int calculateOrbitOffset(vector<CmsInfo> &cmsInfo, vector<info> &totem,
		int *nmatchorb, int *nmatchall, TH1I *hScanOffset,
		TH1I *hScanOffset_match, TH1I *hBcnDelta) {
	const int n_cms = cmsInfo.size();
	CmsInfo info;
	const unsigned int nTotem = totem.size();
	cout << "scanning totem...\n";
	for (unsigned int iScanTotem = 0; iScanTotem < nTotem; ++iScanTotem) { //for every totem event from sample
		cout << "\t" << iScanTotem << endl;
		// checking totem BX :
		if ((totem[iScanTotem]).bcn < dTOTEMoRxBX)
			countTOTwrongBX++;
		int i = 0;
		info = cmsInfo[i];
		int minOrbit = (totem[iScanTotem]).orbit;
		if (info.evt.Orbit < minOrbit) { //this if is to evade the situation that after "while" statement we have i == -1
			// cms < totem - can't be => just pass
			while ((info.evt.Orbit < minOrbit) && (i < n_cms)) { //find first cms event such
				info = cmsInfo[i++]; //that cms.Orbit >= minOrbit
			}
			i--;
		}
		if (info.evt.Orbit < minOrbit) //if we haven't find cms event that cms.Orbit >= minOrbit go to next totem event
			continue;
		int deltaOrbit = info.evt.Orbit - (totem[iScanTotem]).orbit;
		while ((deltaOrbit < maxDeltaOrbit) && (i < n_cms)) { //iterate over all cms events such that cms.orbit belongs to range (totem.orbit; totem.orbit + maxOffset)
			info = cmsInfo[i++];
			deltaOrbit = info.evt.Orbit - (totem[iScanTotem]).orbit;
			// cms is supposed to be in increasing order
			if (deltaOrbit < maxDeltaOrbit) {
				nmatchorb[deltaOrbit]++; //increase the number of possibly matched events for orbit offset == deltaOrbit
				hScanOffset->Fill(deltaOrbit);
				if (info.evt.Bunch - (totem[iScanTotem]).bcn == 1) {
					// cms start bcn with 1, totem with 0
					nmatchall[deltaOrbit]++; //assumes that bunch crossing offset is always 1
					hScanOffset_match->Fill(deltaOrbit);
				}
				hBcnDelta->Fill(info.evt.Bunch - (totem[iScanTotem]).bcn);
			}
		}
	} // end loop totem events
	float maxNmatchall = 0;
	int dOrbit = -1;
	for (int i = 0; i < maxDeltaOrbit; i++) { //find the orbit offset that has the biggest number of possibly matched events
		if (nmatchall[i] > maxNmatchall) {
			maxNmatchall = nmatchall[i];
			dOrbit = i;
		}
	}
	if (dOrbit == -1)
		error("can't find any matched events; exiting...\n");

//print some info on STDOUT:
	cout << "\n\ndOrbit (cms - totem) = " << dOrbit << "; \"eff\" wrt TOTEM = "
			<< float(maxNmatchall) / float(nTotem) << " (" << maxNmatchall
			<< " matched events; " << nmatchorb[dOrbit]
			<< " matched with any bcn)" << endl;

	info = cmsInfo[0];
	cout << "start : (totem : cms)\n";
	cout << (totem[0]).orbit << "\t:\t" << info.evt.Orbit << std::endl;
	info = cmsInfo[n_cms - 1];
	cout << "stop  : (totem : cms)\n";
	cout << (totem[nTotem - 1]).orbit << "\t:\t" << info.evt.Orbit << std::endl;
	cout << "totem = " << nTotem << endl;
	return dOrbit;
}

/**
 * This function searches the orbit offset, creates some histograms presenting possible orbit
 * and bunch offsets and then merges the events that were used for searching.
 *
 * totemFileName - root file containing totem ntuple
 * cmsFileName - root file containing cms ntuple
 * outputFileName - root file where the histograms and sample of merged data will be stored
 */
void merge12(const string& totemFileName, const string& cmsFileName,
		const string &outputFileName) {
	const int verbose = 0; // 0: normal, 1: little,  3: extended
	cout << "Opening Totem File: " << totemFileName << endl;
	TFile* totemFile = TFile::Open(totemFileName.c_str()); //opening input files
	if (!totemFile || totemFile->IsZombie())
		error("Error opening file");
	TTree* tree_totem = (TTree*) totemFile->Get("TotemNtuple");
	if (!tree_totem)
		error("No data tree found");
	cout << "Entries: " << tree_totem->GetEntries() << endl;
	cout << "Opening CMS File: " << cmsFileName << endl;
	TFile* cmsFile = TFile::Open(cmsFileName.c_str());
	if (!cmsFile || cmsFile->IsZombie())
		error("Error opening file");
	TTree *tree_cms = (TTree*) cmsFile->Get("evt");
	cout << "Entries: " << tree_cms->GetEntries() << endl;

	MyEvtId *evtcmsUA = 0; //assigns variables with input files and prepares branches in output file
	MyL1TrigOld *trigcmsUA = 0;
	TFile* output = TFile::Open(outputFileName.c_str(), "RECREATE");
	TTree* outputTree = new TTree("cms_totem", "merged");
	outputTree->Branch("cmsEvtUA", &evtcmsUA);
	outputTree->Branch("cmsTrigUA", &trigcmsUA);
	prepareTotemBranches(totemFileName, outputTree, tree_totem);

	if (verbose > 2) {
		cout << " cms raw sequence: nEv=" << tree_cms->GetEntries() << endl;
		cout << " totem raw sequence: nEv=" << tree_totem->GetEntries() << endl;
	}

	vector<info> totem; //the basic totem data (trigger + bunch, orbit etc) will be stored in this vector
	const int n_totem = 50000; //number of totem events that will be used during orbit offset searching (these events are the first n_totem events chosen from the beggining of the ntuple)
	std::cout << "filling totem vector... Reduced statistics \n";
	for (int i = 0; i < n_totem; ++i) { //prepare the vector of totem events
		std::cout << "->" << i << std::endl;
		tree_totem->GetEntry(i);
		const info info_totem(i, trigData->event_num, trigData->orbit_num, trigData->bunch_num);
		totem.push_back(info_totem);
		if (verbose > 2)
			cout << trigData->event_num << " " << totemEventData->event_no
					<< "\t\t" << totemEventData->optoRx_BX.at(1) << " & "
					<< t2Event->ev_no << " & " << trigData->orbit_num
					<< " !!\n";
	}
	int nmatchorb[maxDeltaOrbit]; //prepare arrays and histograms with obrbit, bunch offsets
	memset(nmatchorb, 0, sizeof(nmatchorb));
	int nmatchall[maxDeltaOrbit];
	memset(nmatchall, 0, sizeof(nmatchall));
	TH1I *hScanOffset = new TH1I("hScanOffset", "hScanOffset",
			maxDeltaOrbit + 10, -10, maxDeltaOrbit);
	hScanOffset->SetXTitle("orbit_{cms}-orbit_{totem}");
	TH1I *hScanOffset_match = new TH1I("hScanOffset_match",
			"hScanOffset_match", maxDeltaOrbit + 10, -10, maxDeltaOrbit);
	hScanOffset_match->SetXTitle("orbit_{cms}-orbit_{totem}");
	hScanOffset_match->SetLineColor(2);
	TH1I *hBcnDelta = new TH1I("hBcnDelta", "hBcnDelta",
			maxDeltaBX * 2, -maxDeltaBX, maxDeltaBX);
	hBcnDelta->SetXTitle("bx_{cms}-bx_{totem}");

	vector<CmsInfo> cmsInfo = getSortedCmsInfo(tree_cms, totem[0].orbit,
			totem[n_totem - 1].orbit); //returns the sorted cms events vector (the events correspond to the orbit range of chosen totem sample)
	unsigned int dOrbit = calculateOrbitOffset(cmsInfo, totem, nmatchorb,
			nmatchall, hScanOffset, hScanOffset_match, hBcnDelta); //found orbit offset

	unsigned int nTotem = totem.size();
	unsigned int mcount = 0;
	unsigned int totem_i;
	unsigned int previous_totem_i = 0;
	unsigned int firstCMSorbit = (totem[0]).orbit + dOrbit;
	unsigned int lastCMSorbit = (totem[nTotem - 1]).orbit + dOrbit;
	unsigned int previousCMSorb = 0;
	unsigned int secondCMS = 0;
	unsigned int countCMS = 0;
	unsigned int possibleTOTEM = 0;
	unsigned int countCMSnomatch = 0;
	unsigned int secondCMSnomatch = 0; // previous of the same orb is found
	unsigned int nextCMSnomatch = 0; // next of the same orb is found
	const int n_cms = cmsInfo.size();
	CmsInfo info;

	for (int cms_i = 0; cms_i < n_cms; cms_i++) { //iterate over corresponding cms events
		info = cmsInfo[cms_i];
		if (cleanBX && !collidingBunch(info.evt.Bunch)) //if event bunch was not colliding go to next event
			continue;
		if (((unsigned int) info.evt.Orbit >= firstCMSorbit) //for cms events in this range we should find corresponding in totem sample
		&& ((unsigned int) info.evt.Orbit <= lastCMSorbit))
			countCMS++;
		unsigned int totem_orbit_corrected = info.evt.Orbit - dOrbit; //the corrected orbit (take into account orbit offset)
		bool match = false; //orbit and bunch
		bool oMatch; //just orbit
		totem_i = previous_totem_i; //assuming cms sorted
		while (!match && (totem_i < nTotem)
				&& ((totem[totem_i]).orbit <= totem_orbit_corrected)) { //looking for corresponding totem event
			oMatch = ((totem[totem_i]).orbit == totem_orbit_corrected); //if orbit matches
			match = (oMatch && ((info.evt.Bunch - (totem[totem_i]).bcn) == 1)); //if orbit matches and bunch numbers differ of excatly 1
			if (match) { // we found corresponding event
				// to take into account several events per orbit
				if (totem_i > maxEvtPerOrbit) {
					previous_totem_i = totem_i - maxEvtPerOrbit;
				} else {
					previous_totem_i = 0;
				}
				mcount++;
				evtcmsUA = &info.evt;
				trigcmsUA = &info.trigg;
				tree_totem->GetEntry(totem_i);
				outputTree->Fill(); //fill the merged (output) ntuple
				// the important part is above, bellow there are only some additional informations
			} else {
				// matched orbit, not matched bcn and totem.bcn<dTOTEMoRxBX :
				if (oMatch && ((totem[totem_i]).bcn < dTOTEMoRxBX)) {
					possibleTOTEM++;
				}
				if (lostEvtChck) {
					if (((unsigned int) info.evt.Orbit >= firstCMSorbit)
							&& ((unsigned int) info.evt.Orbit <= lastCMSorbit)) {
						unsigned int to_min = 0;
						int tb_min = 0;
						if (totem_i > 0) {
							to_min = (totem[totem_i - 1]).orbit;
							tb_min = (totem[totem_i - 1]).bcn;
						}
						unsigned int to_max = 0;
						int tb_max = 0;
						if (totem_i < (nTotem - 1)) {
							to_max = (totem[totem_i + 1]).orbit;
							tb_max = (totem[totem_i + 1]).bcn;
						}
						std::cout << totem_i << "\t\t" << to_min << "-> "
								<< (totem[totem_i]).orbit << " : "
								<< info.evt.Orbit - dOrbit << " <-" << to_max
								<< "\t\t" << tb_min << "-> "
								<< (totem[totem_i]).bcn << " : "
								<< info.evt.Orbit << " <-" << tb_max
								<< std::endl;
					}
				}
			} // end match check
			totem_i++;
		} // end loop totem events
		unsigned int orbit = info.evt.Orbit;
		if (!match && (orbit >= firstCMSorbit) && (orbit <= lastCMSorbit)) {
			countCMSnomatch++;
			if (previousCMSorb == orbit) {
				secondCMSnomatch++;
			}
		}
		if ((previousCMSorb == orbit) && (orbit >= firstCMSorbit)
				&& (orbit <= lastCMSorbit))
			secondCMS++;
		previousCMSorb = orbit;

		if (!match && (orbit >= firstCMSorbit) && (orbit <= lastCMSorbit)
				&& cms_i + 1 < n_cms) {
			info = cmsInfo[cms_i + 1];
			if (previousCMSorb == (unsigned int) info.evt.Orbit) {
				nextCMSnomatch++;
			}
		}
	} // end loop cms events
//print some information on STDOUT
	std::cout << "cleanBX (0/1)    = " << cleanBX << std::endl;
	std::cout << "match xcheck     = " << mcount << std::endl;
	std::cout << "countCMS         = " << countCMS << std::endl;
	std::cout << "eff vrt CMS      = " << float(mcount) / float(countCMS)
			<< std::endl;
	std::cout << "secondCMS        = " << secondCMS << std::endl;
	std::cout << "possibleTOTEM    = " << possibleTOTEM << std::endl;
	std::cout << "countCMSnomatch  = " << countCMSnomatch << std::endl;
	std::cout << "secondCMSnomatch = " << secondCMSnomatch << std::endl;
	std::cout << "nextCMSnomatch   = " << nextCMSnomatch << std::endl;
	std::cout << "countTOTwrongBX  = " << countTOTwrongBX << "(total TOTEM)"
			<< std::endl;
	output->Write();
	output->Flush();
}

/**
 * You can run program by typing:
 * ./merge12 totemNtuplePath cmsNtuplePath outputNtupleName
 */
int main(int argc, char** argv) {
	if (argc < 4)
		error("Please specify TOTEM and CMS data files and output file name.");
	const string totem(argv[1]);
	const string cms(argv[2]);
	const string output(argv[3]);
	merge12(totem, cms, output);
	return 0;
}

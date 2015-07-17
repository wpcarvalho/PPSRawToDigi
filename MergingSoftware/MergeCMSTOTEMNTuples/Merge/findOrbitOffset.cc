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
 * This function compares two cms events. First, it compares orbit numbers and then (if equal)
 * bunch numbers.
 */
bool sortEvt(MyEvtId event1, MyEvtId event2) {
	if (event1.Orbit != event2.Orbit) {
		return (event1.Orbit < event2.Orbit);
	} else {
		return (event1.Bunch < event2.Bunch);
	}
	return true; //this is just to avoid compilator warnings
}

/**
 * This function returns sorted vector of cms events. The vector includes just these cms events that
 * correspond to totem events (according to orbit number + some margin).
 * The sorting is necessary before finding the orbit offset.
 *
 * totemOrbitMin - orbit number of the first totem event that will be used for offset searching
 * totemOrbitMax - orbit number of the last totem event that will be used for offset searching
 */
vector<MyEvtId> getSortedCmsInfo(TTree *tree_cms, unsigned int totemOrbitMin,
		unsigned int totemOrbitMax) {
	const unsigned int orbitMargin = 2000; //just for safety we take some more events
	unsigned int cmsSize = tree_cms->GetEntries();
	MyEvtId *cmsEvent = 0;
	TBranch *cmsEventBranch = checkAndGetBranch(tree_cms, "evtId");
	cmsEventBranch->SetAddress(&cmsEvent);

	vector<MyEvtId> vectorCmsEvents;
	int min = totemOrbitMin - orbitMargin;
	int max = totemOrbitMax + orbitMargin;

	for (unsigned int i = 0; i < cmsSize; i++) { //choose cms events within given orbit range
		cmsEventBranch->GetEntry(i);
		if (cmsEvent->Orbit > min && cmsEvent->Orbit < max) {
			vectorCmsEvents.push_back(*(cmsEvent));
		}
	}
	sort(vectorCmsEvents.begin(), vectorCmsEvents.end(), sortEvt);
	return vectorCmsEvents;
}

/**
 * The aim of this function is to find and return orbit offset. (cms.orbit - totem.orbit).
 * The cmsInfo vector must be sorted.
 * cmsInfo - vector of sorted cms events (these events correspond to totem events (orbit range)
 * totem - vector of totem events (not all events just a sample of a given number (you can set variable in the code "const int n_totem")
 */
void calculateOrbitOffset(vector<MyEvtId> &cmsEvents, vector<info> &totem,
		int *nmatchorb, int *nmatchall, TH1I *hScanOffset,
		TH1I *hScanOffset_match, TH1I *hBcnDelta) {
	const int cmsSize = cmsEvents.size();
	MyEvtId cmsEvent;
	const unsigned int totemSize = totem.size();
	cout << "scanning totem...\n";
	for (unsigned int iScanTotem = 0; iScanTotem < totemSize; ++iScanTotem) { //for every totem event from sample
		cout << "\t" << iScanTotem << endl;
		// checking totem BX :
		if ((totem[iScanTotem]).bcn < dTOTEMoRxBX)
			countTOTwrongBX++;
		int i = 0;
		cmsEvent = cmsEvents[i];
		int minOrbit = (totem[iScanTotem]).orbit;
		if (cmsEvent.Orbit < minOrbit) { //this if is to evade the situation that after "while" statement we have i == -1
			// cms < totem - can't be => just pass
			while ((cmsEvent.Orbit < minOrbit) && (i < cmsSize)) { //find first cms event such
				cmsEvent = cmsEvents[i++]; //that cms.Orbit >= minOrbit
			}
			i--;
		}
		if (cmsEvent.Orbit < minOrbit) //if we haven't find cms event that cms.Orbit >= minOrbit go to next totem event
			continue;
		int deltaOrbit = cmsEvent.Orbit - (totem[iScanTotem]).orbit;
		while ((deltaOrbit < maxDeltaOrbit) && (i < cmsSize)) { //iterate over all cms events such that cms.orbit belongs to range (totem.orbit; totem.orbit + maxOffset)
			cmsEvent = cmsEvents[i++];
			deltaOrbit = cmsEvent.Orbit - (totem[iScanTotem]).orbit;
			// cms is supposed to be in increasing order
			if (deltaOrbit < maxDeltaOrbit) {
				nmatchorb[deltaOrbit]++; //increase the number of possibly matched events for orbit offset == deltaOrbit
				hScanOffset->Fill(deltaOrbit);
				if (cmsEvent.Bunch - (totem[iScanTotem]).bcn == 1) {
					// cms start bcn with 1, totem with 0
					nmatchall[deltaOrbit]++; //assumes that bunch crossing offset is always 1
					hScanOffset_match->Fill(deltaOrbit);
				}
				hBcnDelta->Fill(cmsEvent.Bunch - (totem[iScanTotem]).bcn);
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
			<< float(maxNmatchall) / float(totemSize) << " (" << maxNmatchall
			<< " matched events; " << nmatchorb[dOrbit]
			<< " matched with any bcn)" << endl;

	cmsEvent = cmsEvents[0];
	cout << "start : (totem : cms)\n";
	cout << (totem[0]).orbit << "\t:\t" << cmsEvent.Orbit << endl;
	cmsEvent = cmsEvents[cmsSize - 1];
	cout << "stop  : (totem : cms)\n";
	cout << (totem[totemSize - 1]).orbit << "\t:\t" << cmsEvent.Orbit << endl;
	cout << "totem = " << totemSize << endl;
}

/**
 * This function searches the orbit offset, creates some histograms presenting possible orbit
 * and bunch offsets and then merges the events that were used for searching.
 *
 * totemFileName - root file containing totem ntuple
 * cmsFileName - root file containing cms ntuple
 * outputFileName - root file where the histograms and sample of merged data will be stored
 */
void findOrbitOffset(const string& totemFileName, const string& cmsFileName,
		const string &outputFileName) {
	const int verbose = 0; // 0: normal, 1: little,  3: extended
	cout << "Opening Totem File: " << totemFileName << endl;
	TFile* totemFile = TFile::Open(totemFileName.c_str()); //opening input files
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
	TTree *cmsTree = (TTree*) cmsFile->Get("evt");
	cout << "Entries: " << cmsTree->GetEntries() << endl;

	TFile* output = TFile::Open(outputFileName.c_str(), "RECREATE");
	if (verbose > 2) {
		cout << " cms raw sequence: nEv=" << cmsTree->GetEntries() << endl;
		cout << " totem raw sequence: nEv=" << totemTree->GetEntries() << endl;
	}
	TBranch *totemTriggerData = checkAndGetBranch(totemTree, "trigger_data");
	totemTriggerData->SetAddress(&trigData);

	vector<info> totem; //the basic totem data (trigger + bunch, orbit etc) will be stored in this vector
	const int totemSize = totemTree->GetEntries(); //number of totem events that will be used during orbit offset searching (these events are the first n_totem events chosen from the beggining of the ntuple)
	cout << "filling totem vector... Reduced statistics \n";
	for (int i = 0; i < totemSize; ++i) { //prepare the vector of totem events
		cout << "->" << i << endl;
		totemTriggerData->GetEntry(i);
		const info info_totem(i, trigData->event_num, trigData->orbit_num,
				trigData->bunch_num);
		totem.push_back(info_totem);
		if (verbose > 2)
			cout << trigData->event_num << trigData->orbit_num << " !!\n";
	}
	int nmatchorb[maxDeltaOrbit]; //prepare arrays and histograms with obrbit, bunch offsets
	memset(nmatchorb, 0, sizeof(nmatchorb));
	int nmatchall[maxDeltaOrbit];
	memset(nmatchall, 0, sizeof(nmatchall));
	TH1I *hScanOffset = new TH1I("hScanOffset", "hScanOffset",
			maxDeltaOrbit + 10, -10, maxDeltaOrbit);
	hScanOffset->SetXTitle("orbit_{cms}-orbit_{totem}");
	TH1I *hScanOffset_match = new TH1I("hScanOffset_match", "hScanOffset_match",
			maxDeltaOrbit + 10, -10, maxDeltaOrbit);
	hScanOffset_match->SetXTitle("orbit_{cms}-orbit_{totem}");
	hScanOffset_match->SetLineColor(2);
	TH1I *hBcnDelta = new TH1I("hBcnDelta", "hBcnDelta", maxDeltaBX * 2,
			-maxDeltaBX, maxDeltaBX);
	hBcnDelta->SetXTitle("bx_{cms}-bx_{totem}");

	vector<MyEvtId> cmsEvents = getSortedCmsInfo(cmsTree, totem[0].orbit,
			totem[totemSize - 1].orbit); //returns the sorted cms events vector (the events correspond to the orbit range of chosen totem sample)
	calculateOrbitOffset(cmsEvents, totem, nmatchorb,	nmatchall, hScanOffset, hScanOffset_match, hBcnDelta); //found orbit offset

	output->Write();
	output->Flush();
}

/**
 * You can run program by typing:
 * ./findOrbitOffset totemNtuplePath cmsNtuplePath outputNtupleName
 */
int main(int argc, char** argv) {
	if (argc < 4)
		error("Please specify TOTEM and CMS data files and output file name.");
	const string totem(argv[1]);
	const string cms(argv[2]);
	const string output(argv[3]);
	findOrbitOffset(totem, cms, output);
	return 0;
}

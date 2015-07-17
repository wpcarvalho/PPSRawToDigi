#include <TFile.h>
#include <map>
#include "MyEvtId.h"      // UA tree event desctiption
#include "MyL1Trig.h"    // UA tree event desctiption (arrays)
#include "MyL1TrigOld.h" // UA tree event desctiption (arrays)
#include "MyHLTrig.h"
#include "MyTracks.h"
#include "MyVertex.h"
#include "MyPFJet.h"
#include "MyMuon.h"
#include "MyPFCand.h"
#include "MyCaloTower.h"
#include "utilities.h"

using namespace std;

//it would be nice to move next three procedures (checkAndGetBranch, prepareRomanPotBranches, prepareTotemBranches)
//to a place common for both merge12 and combine programs but there are some memory problems that I
//encountered and currently I don't have time to investigate it

/**
  * This function checks if branch of a given name exists.
  * if branch exists - the pointer is returned
  * if not - the program is terminated with error message
  */
 TBranch* checkAndGetBranch(TTree* tree, std::string branchName) {
         TBranch *branch = tree->GetBranch(branchName.c_str());
         if (!branch) {
                 std::string dotBranchName = branchName + ".";
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
void prepareTotemBranches(const std::string& totemFileName, TTree *outputTree,
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
  * This function merges the totem and cms ntuples. The result is saved in outputFileName.
  * dOrbit is previously found (using merge12 program) orbit offset (as far as we tested
  * it was always 0)
  */
void combineWithUA(const int dOrbit, const string& totemFileName,
                 const string& cmsFileName, const string &outputFileName) {
         cout << "Opening Totem File: " << totemFileName << endl; //opening input ntuples
         TFile* totemFile = TFile::Open(totemFileName.c_str());
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
         TTree * tree_cms = (TTree*) cmsFile->Get("evt");
         unsigned int n_cms = tree_cms->GetEntries();
         cout << "CMS entries: " << tree_cms->GetEntries() << endl;

         MyEvtId *evtcmsUA = 0; //assigns cms variables to input ntuple
         checkAndGetBranch(tree_cms, "evtId")->SetAddress(&evtcmsUA);
         MyL1TrigOld * trigcmsUA = 0;
         checkAndGetBranch(tree_cms, "L1TrigOld")->SetAddress(&trigcmsUA);
         MyHLTrig * hltTrig_cmsUA = 0;
         checkAndGetBranch(tree_cms, "HLTrig")->SetAddress(&hltTrig_cmsUA);
         std::vector<MyTracks> *tracks_cmsUA = 0;
         checkAndGetBranch(tree_cms, "generalTracks")->SetAddress(&tracks_cmsUA);
         std::vector<MyVertex> *vertices_cmsUA = 0;
         checkAndGetBranch(tree_cms, "offlinePrimaryVertices")->SetAddress(
                         &vertices_cmsUA);
         std::vector<MyPFJet> *pfJets_cmsUA = 0;
         checkAndGetBranch(tree_cms, "ak5PFJets")->SetAddress(&pfJets_cmsUA);
         std::vector<MyMuon> *muons_cmsUA = 0;
         checkAndGetBranch(tree_cms, "muons")->SetAddress(&muons_cmsUA);
         std::vector<MyPFCand> *pFlow_cmsUA = 0;
         checkAndGetBranch(tree_cms, "particleFlow")->SetAddress(&pFlow_cmsUA);
         std::vector<MyCaloTower> *caloTowers_cmsUA = 0;
         checkAndGetBranch(tree_cms, "caloTowers")->SetAddress(&caloTowers_cmsUA);

         TFile* output = TFile::Open(outputFileName.c_str(), "RECREATE"); //output ntuple
         TTree* outputTree = new TTree("cms_totem", "merged");
         outputTree->Branch("cmsEvtUA", &evtcmsUA); //preparing branches in output ntuple
         outputTree->Branch("cmsTrigUA", &trigcmsUA);
         outputTree->Branch("cmsHLTTrigUA", &hltTrig_cmsUA);
         outputTree->Branch("cmsTracksUA", &tracks_cmsUA);
         outputTree->Branch("cmsVerticesUA", &vertices_cmsUA);
         outputTree->Branch("cmsPFJetsUA", &pfJets_cmsUA);
         outputTree->Branch("cmsMuonsUA", &muons_cmsUA);
         outputTree->Branch("cmsParticleFlowUA", &pFlow_cmsUA);
         outputTree->Branch("cmsCaloTowersUA", &caloTowers_cmsUA);
         prepareTotemBranches(totemFileName, outputTree, tree_totem);

         vector<info> totem; //this vector contains only basic information about totem events
                                                 //it will speed up merging procedure
         map<unsigned int, vector<int> > orbits; //this map helps in searching for totem event with given orbit number
                                                                                         //key - orbit number, value - vector of indexes of totem events that have this orbit number
         map<unsigned int, vector<int> >::iterator it;
         const int n_totem = tree_totem->GetEntries(); //here you can set how many totem events you would like to merge
         for (int i = 0; i < n_totem; ++i) {
                 tree_totem->GetEntry(i);
                 const info info_totem(i, trigData->event_num, trigData->orbit_num, trigData->bunch_num);
                 totem.push_back(info_totem); //get event from totem tree and add it to the vector
                 unsigned int orbitNumber = trigData->orbit_num;
                 it = orbits.find(orbitNumber);
                 if (it != orbits.end()) { //if keys exist add to the existing vector
                         it->second.push_back(i);
                 } else { //if not create a vector for given orbit number
                         vector<int> vect;
                         vect.push_back(i);
                         orbits.insert(pair<int, vector<int> >(orbitNumber, vect));
                 }
         }
         // merge events and write them to the output ntuple
         unsigned int mcount = 0;
         unsigned int firstCMSorbit = (totem[0]).orbit + dOrbit;
         unsigned int lastCMSorbit = (totem[n_totem - 1]).orbit + dOrbit;
         unsigned int countCMS = 0;
         unsigned int countCMSnomatch = 0;
         int cms_start = 0;
         int cms_end = cms_start + n_cms;
         std::cout << "First, last CMS orbit: " << firstCMSorbit << ", "
                         << lastCMSorbit << std::endl;
         for (int cms_i = cms_start; cms_i < cms_end; cms_i++) { //iterate over cms events, cms ntuples is not sorted
                 tree_cms->GetEntry(cms_i);
                 std::cout << cms_i << "\t\t" << evtcmsUA->Orbit << "\t\t"
                                 << evtcmsUA->Bunch << std::endl;
                 if (cleanBX && !collidingBunch(evtcmsUA->Bunch)) { //if bunch is not colliding
                         std::cout << "Non colliding bunch.. skipping." << std::endl;
                         continue;
                 }
                 if (((unsigned int) evtcmsUA->Orbit >= firstCMSorbit)
                                 && ((unsigned int) evtcmsUA->Orbit <= lastCMSorbit)) { //in this we should be able to find corresponding event in totem vector
                         countCMS++;
                 } else {
                         continue;
                 }
                 unsigned int totem_orbit_corrected = evtcmsUA->Orbit - dOrbit;
                 bool match = false;
                 it = orbits.find(totem_orbit_corrected); //we are looking for totem events with totem_orbit_corrected orbit number
                 if (it != orbits.end()) {
                         vector<int> totemEvents = it->second;
                         for (unsigned int i = 0; (i < totemEvents.size()) && !match; i++) { //iterate over totem events with given orbit number and find one that has a proper bunch crossing
                                 int index = totemEvents[i];
                                 if ((evtcmsUA->Bunch - totem[index].bcn) == 1) {
                                         std::cout << (totem[index]).orbit << "\t\t"
                                                         << evtcmsUA->Orbit << "\t\t"
                                                         << evtcmsUA->Orbit - (totem[index]).orbit << "\t\t"
                                                         << (totem[index]).bcn << "\t\t" << evtcmsUA->Bunch
                                                         << "\t\t" << (evtcmsUA->Bunch - (totem[index]).bcn)
                                                         << std::endl;
                                         mcount++;
                                         tree_totem->GetEntry(index);
                                         outputTree->Fill(); //fill the output ntuple
                                         match = true;
                                 }
                         }
                 }
                 if (!match && ((unsigned int) evtcmsUA->Orbit >= firstCMSorbit)
                                 && ((unsigned int) evtcmsUA->Orbit <= lastCMSorbit)) {
                         countCMSnomatch++;
                 }
         } // end loop cms events
         outputTree->Write();
         output->Close();
         //print some basic information to STDOUT
         std::cout << "cleanBX (0/1)    = " << cleanBX << std::endl;
         std::cout << "match xcheck     = " << mcount << std::endl;
         std::cout << "countCMS         = " << countCMS << std::endl;
         std::cout << "countCMSnomatch  = " << countCMSnomatch << std::endl;
}

/**
  * You can use this program by typing:
  * ./combine totemNtuplePath cmsNtuplePath outputNtuple
  */
int main(int argc, char** argv) {
         if (argc < 4)
                 error("Please specify TOTEM and CMS data files and output file.");
         const string totem(argv[1]);
         const string cms(argv[2]);
         const string outputFileName(argv[3]);
         combineWithUA(0, totem, cms, outputFileName);
         return 0;
}

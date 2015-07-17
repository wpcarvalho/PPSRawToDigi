#include "L1TriggerTotem/CoincidenceChip/interface/RPTriggerAnalyzer.h"
#include "L1TriggerTotem/CoincidenceChip/interface/CoincidenceChip.h"
#include "L1TriggerTotem/CoincidenceChip/src/rootlogon.C"

#include "TotemCondFormats/DataRecord/interface/TotemDAQRecord.h"
#include "TotemCondFormats/DAQInformation/interface/DAQInformationRP.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"

#include "DataFormats/Common/interface/OwnVector.h"
// extern "C" void rootlogon() ;

using namespace std;

// ClassImpl(TriggerStat2);

RPTriggerAnalyzer::RPTriggerAnalyzer(const edm::ParameterSet& conf):
	verbosity(conf.getParameter<unsigned int>("verbosity")),
    fInfoCollector(),
    fPotCollection(),
	fConfig(conf)
   {

	 detTriggerLabel =  conf.getParameter<edm::InputTag>("DetTriggerLabel");
    rootlogon();
}

RPTriggerAnalyzer::~RPTriggerAnalyzer() { 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

// ------------ method called once each job just before starting event loop  ------------
void RPTriggerAnalyzer::beginRun(edm::Run const&, edm::EventSetup const& iSetup){
    if( verbosity > 0) std::cout << "RPTriggerAnalyzer: beginJob" << std::endl;

    std::vector<unsigned int> chosenPotsId = fConfig.getParameter<std::vector<unsigned int> >("chosenPotsId");
   // const edm::ParameterSet params = fConfig.getParameterSet("params");
   //bool useValidTrackCondition            = params.getParameter<bool>("useValidTrackCondition");

    edm::ESHandle< DAQInformationRP > daqInformationRP; 
    iSetup.get<TotemDAQRecord>().get(daqInformationRP);


    for(std::map<DAQInformationRP::RPID, VFATRegisters>::const_iterator it = daqInformationRP->coincidenceIdToRegisters.begin(); it!= daqInformationRP->coincidenceIdToRegisters.end();it++){
        unsigned int rpID=it->first;
        bool test=true;
        for(unsigned int i = 0; i < chosenPotsId.size(); i++){
          // cout << "1=" << chosenPotsId[i]  << "   2=" << rpID << endl;
          if(chosenPotsId[i] == rpID) test=false;
        }
        if(test and chosenPotsId.size()!=0) continue;

        if(verbosity >2) std::cout << " 1=" << it->first << std::endl;
       // TString potLabel="pot"; 
        // potLabel+=rpID; // potLabel+="_";
      //  potLabel = TotRPDetId::RPName(rpID);
        fPotCollection.push_back(new RawVsSimuPotComparator(rpID, verbosity,&fConfig));
        if(verbosity) cout << " id=" << fPotCollection.back()->GetDecRPIdFull() << "  label=" << fPotCollection.back()->GetName() << endl;

        fPotCollection.back()->beginJob();
    }
}

// ------------ method called to for each event  ------------
void RPTriggerAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& iSetup) {
  // const edm::TimeValue_t &timestamp = event.time().value();
  // printf(">>>>> %s\n", ctime((const time_t *) &timestamp));

    fInfoCollector.fNEvents++;

#if 0
    edm::Handle< edm::DetSetVector<RPStripDigi> >  stripDigiCollection;
    event.getByType( stripDigiCollection );

    for(edm::DetSetVector<RPStripDigi>::const_iterator stripDigiIt=stripDigiCollection->begin(); stripDigiIt != stripDigiCollection->end(); ++stripDigiIt) {
        for(edm::DetSet<RPStripDigi>::const_iterator it=stripDigiIt->begin(); it != stripDigiIt->end(); ++it) {
            TotRPDetId rpId(it->GetDetId());
            unsigned int plane = rpId.Detector()/2+1;
          //   cout <<  "rpId=" << rpId.DetectorDecId() << "  plane=" << plane << "  stripNo=" << it->GetStripNo() << endl;
        }
    }
#endif 
#if 0   
    bool test1 = (fInfoCollector.fNEvents == 26);
    bool test  = ((fInfoCollector.fNEvents == 26) or (fInfoCollector.fNEvents == 27));
   // cout << "#Event = " << fInfoCollector.fNEvents  << "   test=" << test << "  test1=" << test1 << endl;
    if(!((fInfoCollector.fNEvents == 26) or (fInfoCollector.fNEvents == 27))) {
       return;
    }else{
    for(edm::DetSetVector<RPStripDigi>::const_iterator stripDigiIt=stripDigiCollection->begin(); stripDigiIt != stripDigiCollection->end(); ++stripDigiIt) {
        for(edm::DetSet<RPStripDigi>::const_iterator it = stripDigiIt->begin(); it != stripDigiIt->end(); ++it) {
            TotRPDetId rpId(it->GetDetId());
            unsigned int plane = rpId.Detector()/2+1;
             cout << "plane=" << plane << "  stripNo=" << it->GetStripNo() << endl;

        }
    }
    }
#endif

//    cout << "ASdf" << endl;
    if(verbosity) std::cout << "### Event " << fInfoCollector.fNEvents    << " ################################ " << std::endl;
//    return;
    // if(verbosity) std::cout << " run = "    << event.run() << std::endl;
    // if(verbosity) std::cout << event.id()  << std::endl;
    
    // Check if event is empty (input to CC consists of zeros)
    edm::Handle<edm::DetSetVector<RPDetTrigger> > inputTrigger;
    event.getByLabel(detTriggerLabel, inputTrigger );
    if(inputTrigger->size()==0) {
        if(verbosity) std::cout << "empty" << endl;
        fInfoCollector.fNEmptyEvents++;
    }

    // create RP CC input bits
    fPotCollection.CreateNewRPCCInput(event, detTriggerLabel.label()); //label as string because InputTag requires c++11.

    for(unsigned int i=0; i < fPotCollection.size(); i++){
        if(verbosity) std::cout << "pot = " << TotRPDetId::RPName(fPotCollection[i]->GetDecRPIdFull())  << "   "  <<   fPotCollection[i]->GetDecRPIdFull() <<std::endl;
        fPotCollection[i]->analyze(event,iSetup);
    }

    if(verbosity) std::cout << endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void RPTriggerAnalyzer::endJob() { 

    const edm::ParameterSet& params = fConfig.getParameterSet("params");
    std::string runInfo             = params.getParameter<std::string>("runInfo");
    std::string runName             = params.getParameter<std::string>("runName");
    fInfoCollector.runName = runName;

	std::vector<std::string> sourceFileNames = params.getUntrackedParameter<std::vector<std::string> >("sourceFileNames");
    // planes = conf.getParameter<unsigned int>("planes");

    if(verbosity) std::cout << setw(70) << setfill('#') << "#" << endl;
    if(verbosity) std::cout << "Run name         :  " << runName << endl;
    if(verbosity) std::cout << "Source file names:  " << endl;
    for(unsigned int i=0; i < sourceFileNames.size(); i++){
        if(verbosity) std::cout << "   " << sourceFileNames[i] << endl;
    }

    if(verbosity) {
        CoincidenceChip coincidenceChip;
        // configure the coincidence chip
        coincidenceChip.configure(params.getParameterSet("coincidenceChipConfig")); 
        std::cout << coincidenceChip.configSummary();
    }
    if(verbosity) std::cout << endl;
    if(verbosity) std::cout << runInfo << endl;
    if(verbosity) std::cout << endl;
    if(verbosity) std::cout << "Stats global:" << endl;
    if(verbosity) std::cout << "NEvents                                   = " << fInfoCollector.fNEvents                                 << endl;
    if(verbosity) std::cout << "NEmptyEvents                              = " << fInfoCollector.fNEmptyEvents                            << endl;
    if(verbosity) std::cout << "NonEmptyEvents                            = " << fInfoCollector.fNEvents - fInfoCollector.fNEmptyEvents  << endl;
    if(verbosity) std::cout << endl;
    for(unsigned int i=0; i< fPotCollection.size(); i++){
        fPotCollection[i]->endJob();
    }

    for(unsigned int i=0; i< fPotCollection.size(); i++){
     //    fPotCollection[i].endJob2();
    }

    std::string outputDir   = fConfig.getParameter<std::string>("outputDir");
    std::string outputFile  = fConfig.getParameter<std::string>("outputFile");

    if(verbosity) std::cout << setw(70) << setfill('#') << "#" << endl;
    gSystem->Exec("mkdir --parents "+TString(outputDir.c_str()));
    TFile * theFile = TFile::Open(TString(outputDir.c_str())+"/"+outputFile.c_str(), "recreate");
    if(!theFile or !theFile->IsWritable()) { cout << "Output file not opened correctly!!" << endl;}
   // edm::OwnVector<RawVsSimuPotComparator> vec;
   // theFile->WriteObject(&vec, "testik");
    fInfoCollector.Write();
    fPotCollection.Write();
#if 0    
    if(verbosity){
        const edm::ParameterSet params = fConfig.getParameterSet("params");
        // std::string runInfo         = params.getParameter<std::string>("runInfo");
        std::string runName            = params.getParameter<std::string>("runName");

        // cout << " Run   POT   cond.   events     TM_raw  TM_simu  TM_notsame   TM_notsameEven TM_!=odd TM_!=EvenOdd" << endl;
        for(unsigned int i=0; i< fPotCollection.size(); i++){
            //fPotCollection[i].fConditionalTriggerStat.Write(fPotCollection[i].GetName());
           // fPotCollection[i]->Write();
            for(unsigned int j = 0; j < fPotCollection[i]->fConditionalTriggerStat.size(); j++){
                const TriggerStat2& trigStat = *(fPotCollection[i]->fConditionalTriggerStat[j]);
                 // trigStat.Write(trigStat.GetName());
#if 0
                TString name2 = fPotCollection[j].GetName()+"_";
                name2+=i;
                name2+="_";
                name2+=j;
                theFile->WriteObject(&trigStat, name2);
                const TriggerStat2* nic;
                theFile->GetObject(name2,nic);
#endif
                const TriggerStat& trigStatRaw  = *fPotCollection[i]->GetPot( kRaw).fConditionalTriggerStat.AssertFind(trigStat.GetName());
                const TriggerStat& trigStatSimu = *fPotCollection[i]->GetPot(kSimu).fConditionalTriggerStat.AssertFind(trigStat.GetName());
                // cout << "trigstat= " << trigStat.GetConditionNEvents() << endl;
                // cout << "raw     = " << trigStatRaw.GetConditionNEvents() << endl;
                // cout << "sim     = " << trigStatSimu.GetConditionNEvents() << endl;
                assert(trigStat.GetConditionNEvents() == trigStatRaw.GetConditionNEvents());
                assert(trigStat.GetConditionNEvents() == trigStatSimu.GetConditionNEvents());

                // cout << runName << " &  " << fPotCollection[i].GetDecRPIdFull() << " & " <<   trigStat.GetName()  << " & " << trigStat.GetConditionNEvents();
                // cout << runName << " &  " << TotRPDetId::RPName(fPotCollection[i].GetDecRPIdFull()) << " & " <<   trigStat.GetName()  << " & " << trigStat.GetConditionNEvents() << " & " ;
                TString forTex = "FORTEX  ";
                cout << forTex << "BEGIN_LABEL" << endl;
                cout << forTex << "runName  " << runName << endl;
                cout << forTex << "potName  " << TotRPDetId::RPName(fPotCollection[i]->GetDecRPIdFull()) << endl;
                cout << forTex << "condition  " << trigStat.GetName() << endl;
                cout << forTex << "nevents  " << trigStat.GetConditionNEvents() << endl;
                cout << forTex << "fConditionANDTriggerRaw  " << trigStatRaw.fConditionANDTrigger << endl;
                cout << forTex << "fConditionANDTriggerSimu  " << trigStatSimu.fConditionANDTrigger << endl;
                cout << forTex << "fConditionANDRawSimuCCoutputNotTheSame  " << trigStat.fConditionANDRawSimuCCoutputNotTheSame << endl;
                cout << forTex << "fConditionANDRawSimuCCoutputNotTheSameEven  " << trigStat.fConditionANDRawSimuCCoutputNotTheSameEven << endl;
                cout << forTex << "fConditionANDRawSimuCCoutputNotTheSameOdd  " << trigStat.fConditionANDRawSimuCCoutputNotTheSameOdd << endl;
                cout << forTex << "fConditionANDRawSimuCCoutputNotTheSameEvenOdd  " << trigStat.fConditionANDRawSimuCCoutputNotTheSameEvenOdd << endl;
                cout << forTex << "fConditionANDTriggerRaw0Simu0  " << trigStat.fConditionANDTriggerRaw0Simu0 << endl;
                cout << forTex << "fConditionANDTriggerRaw1Simu0  " << trigStat.fConditionANDTriggerRaw1Simu0 << endl;
                cout << forTex << "fConditionANDTriggerRaw0Simu1  " << trigStat.fConditionANDTriggerRaw0Simu1 << endl;
                cout << forTex << "fConditionANDTriggerRaw1Simu1  " << trigStat.fConditionANDTriggerRaw1Simu1 << endl;
#if 0
                double nevents = trigStat.GetConditionNEvents(); // keep nevents double or problems with converting int/int to int...
                if(nevents){
                    cout << " & " << fPotCollection[i].GetPot(kRaw).fConditionalTriggerStat.Find(trigStat.GetName())->fConditionANDTrigger/trigStat.GetConditionNEvents(); 
                    cout << " & " << fPotCollection[i].GetPot(kSimu).fConditionalTriggerStat.Find(trigStat.GetName())->fConditionANDTrigger/trigStat.GetConditionNEvents(); 
                    cout << " & " << trigStat.fConditionANDRawSimuCCoutputNotTheSame/trigStat.GetConditionNEvents() << "  \\\\ "<< endl;

//                cout << runName << " &  " << fPotCollection[i].GetDecRPIdFull() << " & " <<   trigStat.GetName()  << " & " << trigStat.GetConditionNEvents();

                PrintText(trigStatRaw.fConditionANDTrigger,"");
                PrintText(trigStatSimu.fConditionANDTrigger,"");
                PrintText(trigStat.fConditionANDRawSimuCCoutputNotTheSame,"");
                PrintText(trigStat.fConditionANDRawSimuCCoutputNotTheSameEven,"");
                PrintText(trigStat.fConditionANDRawSimuCCoutputNotTheSameOdd,"");
                PrintText(trigStat.fConditionANDRawSimuCCoutputNotTheSameEvenOdd,"");
                PrintText(trigStat.fConditionANDTriggerRaw0Simu0,"");
                PrintText(trigStat.fConditionANDTriggerRaw1Simu0,"");
                PrintText(trigStat.fConditionANDTriggerRaw0Simu1,"");
                PrintText(trigStat.fConditionANDTriggerRaw1Simu1,"");
                // statTree.Fill();
                // PrintText(trigStat.fConditionANDRawSimuCCoutputNotTheSame);
                // if(verbosity) std::cout << trigStat.fConditionANDRawSimuCCoutputNotTheSame << " \\\\" << endl;
                 if(verbosity) std::cout << endl;
                }else{
#if 0
                    cout << " & -" ; 
                    cout << " & -" ; 
                    cout << " & - \\\\" << endl;
#endif

//                cout << runName << " &  " << fPotCollection[i].GetDecRPIdFull() << " & " <<   trigStat.GetName()  << " & " << trigStat.GetConditionNEvents() << " & " ;
                PrintText(double(-1.));
                PrintText(double(-1.));
                PrintText(double(-1.),"  \\\\");
                if(verbosity) std::cout << endl;
                }
#endif                    

            }
        }
    }
#endif


//    boost::ptr_vector<TCanvas> can;

#if 0            
    for(unsigned int j = 0; j < fPotCollection.size(); j++){
//        gDirectory->WriteObject(&fPotCollection[j], fPotCollection[j].GetName());
//        theFile->WriteObject(&fPotCollection[j], fPotCollection[j].GetName()+"_2");

        for(unsigned int i = 0; i < fPotCollection[j].fHist2D.size(); i++){
         //   fPotCollection[j].fHist2D[i].Write();
#if 0
            for(unsigned int i=0; i< fPotLabels.size(); i++){
                TString hhname = fHist2D[i].GetName();
                if(hhname.BeginsWith(fPotLabels[i]))
                    cprefix = fPotLabels[i]; 
            }
            can.push_back(new TCanvas(cprefix+fHist2D[i].GetName(),cprefix+fHist2D[i].GetName()));
#endif      

            TString hname = fPotCollection[j].fHist2D[i].GetName();
           //     cout << hname << endl;
           // if(hname.BeginsWith("No")) {
           //     cout << hname << endl;
           //     assert(0);
              // continue;
           // }
            can.push_back(new TCanvas("can_"+hname,"can_"+hname));
            can.back().Divide(2,2);
            can.back().cd(1);
            fPotCollection[j].fHist2D[i].Draw("");
            fPotCollection[j].fHist2D[i].SetDrawOption("COLZTEXT");
            fPotCollection[j].fHist2D[i].GetXaxis()->SetNdivisions(0*10000+0*100+fPotCollection[j].GetNSectors());
            fPotCollection[j].fHist2D[i].GetYaxis()->SetNdivisions(0*10000+0*100+fPotCollection[j].GetRPNPlanesDividedByTwo());

            can.back().cd(2);
            TH1D* projYhist = fPotCollection[j].fHist2D[i].ProjectionY("ProjectionY");
            projYhist->Draw("");
            projYhist->GetYaxis()->SetRangeUser(0,projYhist->GetMaximum()*1.1);
            //fHist2D[i].Write();

            can.back().cd(3);
            TH1D* projXhist = fPotCollection[j].fHist2D[i].ProjectionX("ProjectionX");
            projXhist->Draw("");
            fPotCollection[j].fHist2D[i].Write();
        }
    }
#endif    

#if 0
    std::vector<TString> stripOrSector;
    TString stripName = "_StripsON"; 
    TString sectorName = "_SectorsON"; 
    stripOrSector.push_back(stripName);
    stripOrSector.push_back(sectorName);

    vector<TCanvas*> stripOrSectorCanvas;
    vector<unsigned int> stripOrSectorCanCounter;
    stripOrSectorCanCounter.push_back(1);
    stripOrSectorCanCounter.push_back(1);

    for(unsigned int potIndex=0; potIndex< fPotCollection.size(); potIndex++){
        for(unsigned int i=0; i < stripOrSector.size(); i++){
            for(unsigned int j=0; j < fPotCollection[potIndex].fConditionalTriggerStat.size(); j++){
                TString conditionName = fPotCollection[potIndex].fConditionalTriggerStat[j].GetName();

                can.push_back(new TCanvas(fPotCollection[potIndex].GetName()+"_"+conditionName+stripOrSector[i],""));
                can.back().SetTitle( can.back().GetName() );
                can.back().Divide(2,fPotCollection[potIndex].GetRPNPlanesDividedByTwo());
                stripOrSectorCanvas.push_back(&can.back());
                stripOrSectorCanCounter[i]=1;

                for(unsigned int hi = 0; hi < fPotCollection[potIndex].fHist1D.size(); hi++){
                    TString hname = fPotCollection[potIndex].fHist1D[hi].GetName();
                    // cout << "asdf= " << hname << "  " << stripName << endl;
                    if(!hname.BeginsWith(fPotCollection[potIndex].GetName())) continue;
                    if(!hname.Contains(conditionName)) continue;
                    if(hname.Contains(stripName) and !strcmp(stripOrSector[i], stripName)){
                        stripOrSectorCanvas.back()->cd(stripOrSectorCanCounter[i]++); 
                        fPotCollection[potIndex].fHist1D[hi].Draw("");
                        fPotCollection[potIndex].fHist1D[hi].GetXaxis()->SetNdivisions(0*10000+0*100+fPotCollection[potIndex].GetNStrips());

                    }else if(hname.Contains(sectorName) and !strcmp(stripOrSector[i], sectorName)){
                        // cout << "Asasdf sonda = " <<  stripOrSectorCanCounter[i] << endl;
                        stripOrSectorCanvas.back()->cd(stripOrSectorCanCounter[i]++); 
                        fPotCollection[potIndex].fHist1D[hi].Draw("");
                        fPotCollection[potIndex].fHist1D[hi].GetXaxis()->SetNdivisions(0*10000+0*100+fPotCollection[potIndex].GetNSectors());
                    }

                    //  }else if(hname.Contains("FailingSectors")){
                    //      failingSectorCanvas.cd(failingSectorCanCounter++); 
                    //      fHist1D[i].Draw("");
                    //      fHist1D[i].GetXaxis()->SetNdivisions(0*10000+0*100+fNSectors);
                    //  }
            }
         }
      }
   } 


  for(unsigned int i=0; i < can.size(); i++){
    can[i].Write();

    TString figname = TString(can[i].GetName())+".pdf";
    //can[i]->SaveAs(fname);
    //can[i]->Print(fname+".png");
    // can[i]->Print(figname+".ps");
    TString potName="";
    TString conditionName="";
    TString tmpfigname=figname;
    for(unsigned int i=0; i< fPotCollection.size(); i++){
        if(figname.BeginsWith( fPotCollection[i].GetName())){
            potName = fPotCollection[i].GetName(); 
            tmpfigname.ReplaceAll(potName+"_","");

            for(unsigned int j=0; j < fPotCollection[i].fConditionalTriggerStat.size(); j++){
                if(tmpfigname.BeginsWith(fPotCollection[i].fConditionalTriggerStat[j].GetName())){
                    conditionName = fPotCollection[i].fConditionalTriggerStat[j].GetName();
                }
            }
        }
    }
#endif

#if 0
    TString mydir = TString(outputDir.c_str())+"/"+potName+ "/"+conditionName;
     cout << "  mdir= " << mydir << "/" << figname << endl;
    //gPad->Print(figname+".ps");
    //cout << " asdf:= " << figname.Strip(TString::kTrailing, '_')  << endl;
    //cout << " asdf:= " << figname.Strip(TString::kLeading, '_')  << endl;
    //cout << " asdf:= " << figname.Strip(TString::kBoth, ' ')  << endl;
    //assert(0);
    gSystem->Exec("pdftk " + figname + " cat 1L output out.pdf; mkdir --parents "+mydir+"; mv out.pdf "+mydir+"/"+figname+"; rm "+figname);
    // assert(0);

    // gSystem->Exec("ps2eps --force --rotate=+ "+figname+".ps");
    // gSystem->Exec("epstopdf "+figname+".eps");
    //gSystem->Exec("ps2pdf "+figname);
  }
#endif    

  theFile->Close();
}

#if 0
void RPTriggerAnalyzer::PrintText(double val, const TString& separator) const{
    if(verbosity) {
        std::cout <<  setfill(' ') << setw(8);
        if(val == -1.) { std::cout << "NaN" << separator; return;}
        if(val < 1.1 and false){
            if(val >=0.994){ 
                std::cout << setprecision(6) <<  fixed << val << separator;
                // val = 0.994;
                // cout << setprecision(2) <<  fixed << val << separator;
            }else{
                // std::cout << setprecision(4) <<  fixed << val << separator;
                std::cout << setprecision(6) <<  fixed << val << separator;
            }
        }else{
            std::cout << val << separator;
        }
    }
}
#endif

// define this as a plug-in
DEFINE_FWK_MODULE(RPTriggerAnalyzer);

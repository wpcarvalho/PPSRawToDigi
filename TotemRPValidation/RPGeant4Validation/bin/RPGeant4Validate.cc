#include <iostream>
#include "stdlib.h"
#include "TotemRPValidation/BaseValidationClasses/interface/BaseHistogramManager.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "TotemRPValidation/RPGeant4Validation/interface/RPOscarSimReader.h"
#include "TotemRPValidation/RPGeant4Validation/interface/RPPSimHitDebugInfo.h"
//#include "TotemRPValidation/RPGeant4Validation/interface/RPEventIO.h"
//#include "TotemRPValidation/RPGeant4Validation/interface/RPEvent.h"
#include "TotemRPValidation/RPGeant4Validation/interface/RPGeant4Validation.h"

#include "TFile.h"
#include "TApplication.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TBrowser.h"
#include "TStyle.h"
#include <vector>


int main(int argc, char *args[])
{
//  BaseHistogramManager b("/test1/test2/");
//  TFile *f = TFile::Open("test.root", "recreate");
//  b.WriteRootFile(f);
  
  if(argc<3)
  {
    std::cout<<"RPGeant4Validate <input_simwatcher_file.root> <dest_histogram_output_file.root>"
      <<std::endl<<"needed!"<<std::endl;
    exit(0);
  }
  
  //  TApplication theApp("App", &argc, argv);
  gStyle->SetOptStat(1111111);
  
  RPGeant4Validation validation_analyser;
  std::cout<<"RPGeant4Validation has been created."<<std::endl;
 
  std::vector<PSimHit> input;
  
  std::cout<<"RPOscarSimReader is being created, root file:"<<args[1]<<std::endl;
  RPOscarSimReader reader = RPOscarSimReader(args[1]);
  std::cout<<"RPOscarSimReader has been created, root file:"<<args[1]<<std::endl;
  std::cout<<"Data processing started..."<<std::endl;

  for(int i=0; reader.RetrieveNewEvent(); i++)
  {
    validation_analyser.FillHistogramms(reader.GetRetrievedEvent(), reader.GetRetrievedEventDebugInfo());

    if(i%50==0)
      std::cout<<i<<std::endl;
  }
    
  std::cout<<"WriteHistograms..."<<std::endl;
  validation_analyser.WriteHistograms(args[2]);
  
//  TBrowser b;
//  cout<<"To terminate press <ctrl>+<c>"<<RPHepPDTWrapper::GetName(2212)<<endl;
//  theApp.Run();
  return 0;
}


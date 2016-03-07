#include "TotemRPValidation/BaseValidationClasses/interface/BaseHistogramManager.h"
#include "TFile.h"
#include <iostream>


BaseHistogramManager::BaseHistogramManager(const std::string &root_file_full_path)
{
  root_file_full_path_=root_file_full_path;
  std::string s;
  
  for(std::string::const_iterator it=root_file_full_path.begin(); it!=root_file_full_path.end();
    ++it)
  {
    if(*it!='/')
      s.push_back(*it);
    
    if(*it=='/' && s.size()>0)
    {
      root_file_full_path_segments_.push_back(s);
      s="";
      //std::cout<<s<<std::endl;
    }
  }
}

void BaseHistogramManager::Finalize()
{
}


void BaseHistogramManager::WriteRootFile(TFile *f)
{
  if(!f)
    return;

  //std::cout<<root_file_full_path_<<" being written"<<std::endl;
    
  Finalize();
    
  f->cd();
  gDirectory->cd("/");
  
  for(unsigned int i=0; i<root_file_full_path_segments_.size(); ++i)
  {
    try{
      gDirectory->cd(root_file_full_path_segments_[i].c_str());
    } catch(...){ 
      gDirectory->mkdir(root_file_full_path_segments_[i].c_str());
      gDirectory->cd(root_file_full_path_segments_[i].c_str());
    }
  }
  
  for(unsigned int i=0; i<reg_hists_.size(); ++i)
  {
    reg_hists_[i]->Write();
//    std::cout<<reg_hists_[i]->GetName()<<", "<<bytes<<" written"<<std::endl;
  }
}


void BaseHistogramManager::RegisterHistogram(TNamed &h)
{
  //h.SetDirectory(0);
  reg_hists_.push_back((TNamed*) &h);
}

#ifndef TotemRPValidation_BaseValidationClasses_BaseComponentClass_h
#define TotemRPValidation_BaseValidationClasses_BaseComponentClass_h

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "TFile.h"
#include "TH1.h"

class BaseHistogramManager
{
  public:
    BaseHistogramManager(const std::string &root_file_full_path);
    virtual ~BaseHistogramManager() {}
    void WriteRootFile(TFile *f);
    void RegisterHistogram(TNamed &h);
    virtual void Finalize();
  private:
    std::string root_file_full_path_;
    std::vector<std::string> root_file_full_path_segments_;
    std::vector<TNamed*> reg_hists_;
};

#endif

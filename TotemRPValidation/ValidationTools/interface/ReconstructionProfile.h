#ifndef TotemRPValidation_ValidationTools_ReconstructionProfile_h
#define TotemRPValidation_ValidationTools_ReconstructionProfile_h

//#define debug_ReconstructionProfile__

#include "TROOT.h"
#include "TAxis.h"
#include <vector>
#include "TH1.h"

struct FitParams
{
  double mean;
  double sigma;
  double chisq_over_N;
  int entries;
  bool valid;
  
  FitParams() : mean(0.0), sigma(0.0), chisq_over_N(0.0), entries(0), 
        valid(false) {}
};

class ReconstructionProfile : public TNamed
{
  public:
    ReconstructionProfile();
    ReconstructionProfile(const char* name, const char* title, Int_t nbinsx, 
          Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, 
          Int_t nbinsz, Double_t zlow, Double_t zup, Option_t* option = "");
    const ReconstructionProfile & operator=(const ReconstructionProfile &in);
    void Fill(Double_t x, Double_t y, Double_t z, Double_t t);
    int GetXBinId(double x) {return x_axis_.FindBin(x);}
    int GetYBinId(double y) {return y_axis_.FindBin(y);}
    int GetZBinId(double z) {return z_axis_.FindBin(z);}
    double Mean(double x1, double x2, double y1, double y2, 
          double z1, double z2, bool all=false);
    double Mean(bool all=false);
    double Sigma(double x1, double x2, double y1, double y2, 
          double z1, double z2, bool all=false);
    double Sigma(bool all=false);
    void FillHistogram(TH1 &h, double x1, double x2, double y1, 
          double y2, double z1, double z2, bool all=false);
    void FillHistogram(TH1 &h, bool all=false);
    FitParams FitStatisticalEstimators(double x1, double x2, double y1, 
          double y2, double z1, double z2, bool all=false);
    FitParams FitStatisticalEstimators(bool all=false);
    int Entries(double x1, double x2, double y1, double y2, 
          double z1, double z2, bool all=false);
    int Entries(bool all=false);
    const TAxis & GetXAxis() const {return x_axis_;}
    const TAxis & GetYAxis() const {return y_axis_;}
    const TAxis & GetZAxis() const {return z_axis_;}
    
  private:
    unsigned int ThreeDimToLinear(unsigned int x, unsigned int y, unsigned int z)
    {
      unsigned int xbins = x_axis_.GetNbins()+2;
      unsigned int ybins = y_axis_.GetNbins()+2;
      return x + xbins*y + xbins*ybins*z;
    }
    unsigned int LinearToThreeDimX(unsigned int lin)
    {
      unsigned int xbins = x_axis_.GetNbins()+2;
      unsigned int ybins = y_axis_.GetNbins()+2;
      return (lin%(xbins*ybins))%xbins;
    }
    unsigned int LinearToThreeDimY(unsigned int lin)
    {
      unsigned int xbins = x_axis_.GetNbins()+2;
      unsigned int ybins = y_axis_.GetNbins()+2;
      return (lin%(xbins*ybins))/xbins;
    }
    unsigned int LinearToThreeDimZ(unsigned int lin)
    {
      unsigned int xbins = x_axis_.GetNbins()+2;
      unsigned int ybins = y_axis_.GetNbins()+2;
      return lin/(xbins*ybins);
    }
    
    FitParams FitHistogram(TH1 &h);
    void SetHistogramRange(TH1 &h);
    
    typedef std::vector< std::vector< double > > contents_type;
    
    void AllocateContents();
    TAxis x_axis_;
    TAxis y_axis_;
    TAxis z_axis_;
    
    contents_type contents;
#ifdef debug_ReconstructionProfile__
    static int fitted_hist_id;
#endif
    
    ClassDef(ReconstructionProfile,1)
};



#endif /*TotemRPValidation_ValidationTools_ReconstructionProfile_h*/

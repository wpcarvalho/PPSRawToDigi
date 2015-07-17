#ifndef TotemRPValidation_ValidationTools_ReconstructionProfile6D_h
#define TotemRPValidation_ValidationTools_ReconstructionProfile6D_h

#include "TROOT.h"
#include "TAxis.h"
#include <vector>
#include "TH1.h"


struct entry_point
{
  double x0, x1, y0, y1, z0, z1;
};

struct entry_range
{
  entry_point down, up;
};


class ReconstructionProfile6D : public TNamed
{
  public:
    ReconstructionProfile6D();
    ReconstructionProfile6D(const char* name, const char* title, Int_t nbinsx, 
          Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, 
          Int_t nbinsz, Double_t zlow, Double_t zup, Option_t* option = "");
    const ReconstructionProfile6D & operator=(const ReconstructionProfile6D &in);
    void Fill(const entry_point &p, Double_t t);
    int GetXBinId(double x) {return x_axis_.FindBin(x);}
    int GetYBinId(double y) {return y_axis_.FindBin(y);}
    int GetZBinId(double z) {return z_axis_.FindBin(z);}
    double Mean(const entry_range &r);
    double Mean();
    double Sigma(const entry_range &r);
    double Sigma();
    void FillHistogram(TH1 &h, const entry_range &r);
    void FillHistogram(TH1 &h);
    int Entries(const entry_range &r);
    int Entries();
    
  private:
    entry_range GetFullRange() const;
    //[x1][x2][y1][y2][z1][z2]
    typedef std::vector< std::vector< std::vector< std::vector<std::vector< 
          std::vector< std::vector<double> > > > > > > contents_type;
    
    void AllocateContents();
    TAxis x_axis_;
    TAxis y_axis_;
    TAxis z_axis_;
    
    contents_type contents;
    
    //ClassDef(ReconstructionProfile6D,1)
};



#endif /*TotemRPValidation_ValidationTools_ReconstructionProfile6D_h*/

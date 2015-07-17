#include "TotemRPValidation/ValidationTools/interface/ReconstructionProfile.h"

#include "TROOT.h"
#include "TMath.h"
//#include "TH3D.h"
#include "TProfile3D.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>

ClassImp(ReconstructionProfile)

ReconstructionProfile::ReconstructionProfile()
 : TNamed()
{
}


const ReconstructionProfile & ReconstructionProfile::operator=(const ReconstructionProfile &in)
{
  TNamed::operator=(in);
  x_axis_.Set(in.x_axis_.GetNbins(), in.x_axis_.GetXmin(), in.x_axis_.GetXmax());
  y_axis_.Set(in.y_axis_.GetNbins(), in.y_axis_.GetXmin(), in.y_axis_.GetXmax());
  z_axis_.Set(in.z_axis_.GetNbins(), in.z_axis_.GetXmin(), in.z_axis_.GetXmax());
  contents = in.contents;
  return in;
}


ReconstructionProfile::ReconstructionProfile(const char* name, const char* title, Int_t nbinsx, 
      Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, 
      Int_t nbinsz, Double_t zlow, Double_t zup, Option_t* option)
 : TNamed(name, title), 
   x_axis_(nbinsx, xlow, xup),
   y_axis_(nbinsy, ylow, yup),
   z_axis_(nbinsz, zlow, zup)
{
  AllocateContents();
}


void ReconstructionProfile::AllocateContents()
{
  unsigned int xbins = x_axis_.GetNbins()+2;
  unsigned int ybins = y_axis_.GetNbins()+2;
  unsigned int zbins = z_axis_.GetNbins()+2;
  
  contents.resize(xbins*ybins*zbins);
}


void ReconstructionProfile::Fill(Double_t x, Double_t y, Double_t z, Double_t t)
{
  int binx,biny,binz;
  binx = x_axis_.FindBin(x);
  biny = y_axis_.FindBin(y);
  binz = z_axis_.FindBin(z);
  
  contents[ThreeDimToLinear(binx, biny, binz)].push_back(t);
}


void ReconstructionProfile::FillHistogram(TH1 &h, double x1, double x2, 
      double y1, double y2, double z1, double z2, bool all)
{
  int x1bin, x2bin, y1bin, y2bin, z1bin, z2bin;
  x1bin = x_axis_.FindBin(x1);
  x2bin = x_axis_.FindBin(x2);
  y1bin = y_axis_.FindBin(y1);
  y2bin = y_axis_.FindBin(y2);
  z1bin = z_axis_.FindBin(z1);
  z2bin = z_axis_.FindBin(z2);

  if(!all)
  {
    if(x1bin==0)
      x1bin++;
    if(y1bin==0)
      y1bin++;
    if(z1bin==0)
      z1bin++;
    if(x2bin==x_axis_.GetNbins()+1)
      x2bin--;
    if(y2bin==y_axis_.GetNbins()+1)
      y2bin--;
    if(z2bin==z_axis_.GetNbins()+1)
      z2bin--;
  }
  
  for(int x_id=x1bin; x_id<=x2bin; ++x_id)
  {
    for(int y_id=y1bin; y_id<=y2bin; ++y_id)
    {
      for(int z_id=z1bin; z_id<=z2bin; ++z_id)
      {
        for(unsigned int content=0; content<contents[ThreeDimToLinear(x_id, y_id, z_id)].size(); ++content)
        {
          h.Fill(contents[ThreeDimToLinear(x_id, y_id, z_id)][content]);
        }
      }
    }
  }
}


void ReconstructionProfile::FillHistogram(TH1 &h, bool all)
{
  double x1 = x_axis_.GetXmin();
  double x2 = x_axis_.GetXmax();
  double y1 = y_axis_.GetXmin();
  double y2 = y_axis_.GetXmax();
  double z1 = z_axis_.GetXmin();
  double z2 = z_axis_.GetXmax();
  
  if(all)
  {
    x1-=1;
    x2+=1;
    y1-=1;
    y2+=1;
    z1-=1;
    z2+=1;
  }
      
  FillHistogram(h, x1, x2, y1, y2, z1, z2, all);
}


double ReconstructionProfile::Mean(double x1, double x2, double y1, double y2, 
      double z1, double z2, bool all)
{
  int x1bin, x2bin, y1bin, y2bin, z1bin, z2bin;
  x1bin = x_axis_.FindBin(x1);
  x2bin = x_axis_.FindBin(x2);
  y1bin = y_axis_.FindBin(y1);
  y2bin = y_axis_.FindBin(y2);
  z1bin = z_axis_.FindBin(z1);
  z2bin = z_axis_.FindBin(z2);

  if(!all)
  {
    if(x1bin==0)
      x1bin++;
    if(y1bin==0)
      y1bin++;
    if(z1bin==0)
      z1bin++;
    if(x2bin==x_axis_.GetNbins()+1)
      x2bin--;
    if(y2bin==y_axis_.GetNbins()+1)
      y2bin--;
    if(z2bin==z_axis_.GetNbins()+1)
      z2bin--;
  }
  
  double sum = 0.0;
  int count = 0;
  
  for(int x_id=x1bin; x_id<=x2bin; ++x_id)
  {
    for(int y_id=y1bin; y_id<=y2bin; ++y_id)
    {
      for(int z_id=z1bin; z_id<=z2bin; ++z_id)
      {
        for(unsigned int content=0; content<contents[ThreeDimToLinear(x_id, y_id, z_id)].size(); ++content)
        {
          double val = contents[ThreeDimToLinear(x_id, y_id, z_id)][content];
          if(val>0.0 || val<=0.0)
          {
            sum += val;
            ++count;
//            std::cout<<x_id<<" "<<y_id<<" "<<z_id<<" "<<count<<" "<<sum<<" "<<
//            val<<std::endl;
          }
        }
      }
    }
  }
  return sum/count;
}


double ReconstructionProfile::Mean(bool all)
{
  double x1 = x_axis_.GetXmin();
  double x2 = x_axis_.GetXmax();
  double y1 = y_axis_.GetXmin();
  double y2 = y_axis_.GetXmax();
  double z1 = z_axis_.GetXmin();
  double z2 = z_axis_.GetXmax();

  if(all)
  {
    x1-=1;
    x2+=1;
    y1-=1;
    y2+=1;
    z1-=1;
    z2+=1;
  }

  return Mean(x1, x2, y1, y2, z1, z2, all);
}


double ReconstructionProfile::Sigma(double x1, double x2, double y1, double y2, 
      double z1, double z2, bool all)
{
  int x1bin, x2bin, y1bin, y2bin, z1bin, z2bin;
  x1bin = x_axis_.FindBin(x1);
  x2bin = x_axis_.FindBin(x2);
  y1bin = y_axis_.FindBin(y1);
  y2bin = y_axis_.FindBin(y2);
  z1bin = z_axis_.FindBin(z1);
  z2bin = z_axis_.FindBin(z2);

  if(!all)
  {
    if(x1bin==0)
      x1bin++;
    if(y1bin==0)
      y1bin++;
    if(z1bin==0)
      z1bin++;
    if(x2bin==x_axis_.GetNbins()+1)
      x2bin--;
    if(y2bin==y_axis_.GetNbins()+1)
      y2bin--;
    if(z2bin==z_axis_.GetNbins()+1)
      z2bin--;
  }
  
  double sum = 0.0;
  double sum2 = 0.0;
  int count = 0;
  
  for(int x_id=x1bin; x_id<=x2bin; ++x_id)
  {
    for(int y_id=y1bin; y_id<=y2bin; ++y_id)
    {
      for(int z_id=z1bin; z_id<=z2bin; ++z_id)
      {
        for(unsigned int content=0; content<contents[ThreeDimToLinear(x_id, y_id, z_id)].size(); ++content)
        {
          double val = contents[ThreeDimToLinear(x_id, y_id, z_id)][content];
          if(val>0.0 || val<=0.0)
          {
            sum += val;
            sum2 += val*val;
            ++count;
          }
        }
      }
    }
  }
  double mean_x = sum/count;
  double mean_x2 = sum2/count;
  return TMath::Sqrt(mean_x2 - mean_x*mean_x);
}


double ReconstructionProfile::Sigma(bool all)
{
  double x1 = x_axis_.GetXmin();
  double x2 = x_axis_.GetXmax();
  double y1 = y_axis_.GetXmin();
  double y2 = y_axis_.GetXmax();
  double z1 = z_axis_.GetXmin();
  double z2 = z_axis_.GetXmax();

  if(all)
  {
    x1-=1;
    x2+=1;
    y1-=1;
    y2+=1;
    z1-=1;
    z2+=1;
  }

  return Sigma(x1, x2, y1, y2, z1, z2, all);
}


int ReconstructionProfile::Entries(double x1, double x2, double y1, double y2, 
      double z1, double z2, bool all)
{
  int x1bin, x2bin, y1bin, y2bin, z1bin, z2bin;
  x1bin = x_axis_.FindBin(x1);
  x2bin = x_axis_.FindBin(x2);
  y1bin = y_axis_.FindBin(y1);
  y2bin = y_axis_.FindBin(y2);
  z1bin = z_axis_.FindBin(z1);
  z2bin = z_axis_.FindBin(z2);
  
  if(!all)
  {
    if(x1bin==0)
      x1bin++;
    if(y1bin==0)
      y1bin++;
    if(z1bin==0)
      z1bin++;
    if(x2bin==x_axis_.GetNbins()+1)
      x2bin--;
    if(y2bin==y_axis_.GetNbins()+1)
      y2bin--;
    if(z2bin==z_axis_.GetNbins()+1)
      z2bin--;
  }
  
  int count = 0;
  
  for(int x_id=x1bin; x_id<=x2bin; ++x_id)
  {
    for(int y_id=y1bin; y_id<=y2bin; ++y_id)
    {
      for(int z_id=z1bin; z_id<=z2bin; ++z_id)
      {
        for(unsigned int content=0; content<contents[ThreeDimToLinear(x_id, y_id, z_id)].size(); ++content)
        {
          double val = contents[ThreeDimToLinear(x_id, y_id, z_id)][content];
          if(val<0.0 || val>=0.0)
          {
            count++;
          }
        }
      }
    }
  }
  
  return count;
}


int ReconstructionProfile::Entries(bool all)
{
  double x1 = x_axis_.GetXmin();
  double x2 = x_axis_.GetXmax();
  double y1 = y_axis_.GetXmin();
  double y2 = y_axis_.GetXmax();
  double z1 = z_axis_.GetXmin();
  double z2 = z_axis_.GetXmax();

  if(all)
  {
    x1-=1;
    x2+=1;
    y1-=1;
    y2+=1;
    z1-=1;
    z2+=1;
  }

  return Entries(x1, x2, y1, y2, z1, z2, all);
}


FitParams ReconstructionProfile::FitStatisticalEstimators(double x1, 
    double x2, double y1, double y2, double z1, double z2, bool all)
{
  int expected_entries = Entries(x1, x2, y1, y2, z1, z2);
  if(expected_entries<50)
    return FitParams();
    
  double expected_mean = Mean(x1, x2, y1, y2, z1, z2);
  double expected_sigma = Sigma(x1, x2, y1, y2, z1, z2);
  TH1D temp("temp_rp__", "temp_rp__", 2000, expected_mean-6*expected_sigma, expected_mean+6*expected_sigma);
  //temp.SetBit(TH1::kCanRebin);
  FillHistogram(temp, x1, x2, y1, y2, z1, z2, all);
  return FitHistogram(temp);
}


FitParams ReconstructionProfile::FitStatisticalEstimators(bool all)
{
  double x1 = x_axis_.GetXmin();
  double x2 = x_axis_.GetXmax();
  double y1 = y_axis_.GetXmin();
  double y2 = y_axis_.GetXmax();
  double z1 = z_axis_.GetXmin();
  double z2 = z_axis_.GetXmax();

  if(all)
  {
    x1-=1;
    x2+=1;
    y1-=1;
    y2+=1;
    z1-=1;
    z2+=1;
  }

  return FitStatisticalEstimators(x1, x2, y1, y2, z1, z2, all);
}


FitParams ReconstructionProfile::FitHistogram(TH1 &h)
{
  FitParams par;
  
  if(h.GetEntries()<50)
    return par;
    
  int count = 0;
  while(h.GetMaximum()<0.05*h.GetEntries() && count<100)
  {
    h.Rebin();
    ++count;
  }
  if(h.GetMaximum()<5)
    return par;
  
  //SetHistogramRange(h);
  TF1 biased_gauss("biased_gauss", "[0]*exp(-(x-[1])*(x-[1])/[2]/[2])+[3]", -10, 10);
  
  double max = h.GetMaximum();
  double mean = h.GetXaxis()->GetBinCenter(h.GetMaximumBin());
  double spread = h.GetRMS()/10.0;
  double bias = 1.0;
  biased_gauss.SetParameters(max, mean, spread, bias);
  
  //h.Fit("gaus", "0Q");
  h.Fit(&biased_gauss, "0Q");
  
#ifdef debug_ReconstructionProfile__
  new TCanvas();
  char new_name[1000];
  sprintf(new_name, "fitted_%03d", fitted_hist_id);
  fitted_hist_id++;
  TH1 *h1 = (TH1 *) h.Clone(new_name); 
  h1->Draw();
#endif
  
//  par.mean = h.GetFunction("gaus")->GetParameter(1);
//  par.sigma = h.GetFunction("gaus")->GetParameter(2);
//  par.chisq_over_N = h.GetFunction("gaus")->GetChisquare() 
//        / h.GetFunction("gaus")->GetNDF();
//  par.entries = (int) h.GetEntries();
//  par.valid = true;

    par.mean = biased_gauss.GetParameter(1);
    par.sigma = biased_gauss.GetParameter(2);
    par.chisq_over_N = biased_gauss.GetChisquare() 
          / biased_gauss.GetNDF();
    par.entries = (int) h.GetEntries();
    par.valid = true;
  
  
  return par;
}

void ReconstructionProfile::SetHistogramRange(TH1 &h)
{
  double max = h.GetMaximum();
  double threshold = max*0.1;
  double wsumx = 0.0;
  double wsumx2 = 0.0;
  double sum = 0.0;
  for(int i=1; i<=h.GetNbinsX(); i++)
  {
    double bin_content = h.GetBinContent(i);
    if(bin_content>=threshold)
    {
      wsumx+=bin_content*h.GetBinCenter(i);
      wsumx2+=bin_content*h.GetBinCenter(i)*h.GetBinCenter(i);
      sum+=bin_content;
    }
  }
  double meanx = wsumx/sum;
  double meanx2 = wsumx2/sum;
  double varx = meanx2 - meanx*meanx;
  double range = 10*TMath::Sqrt(varx);
  int bin_beg = h.GetXaxis()->FindBin(meanx-range);
  int bin_end = h.GetXaxis()->FindBin(meanx+range);
  h.GetXaxis()->SetRange(bin_beg, bin_end);
}

#ifdef debug_ReconstructionProfile__
int ReconstructionProfile::fitted_hist_id = 1;
#endif

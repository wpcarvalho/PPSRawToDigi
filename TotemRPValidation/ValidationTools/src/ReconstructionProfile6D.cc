#include "TotemRPValidation/ValidationTools/interface/ReconstructionProfile6D.h"

#include "TROOT.h"
#include "TMath.h"
//#include "TH3D.h"
#include "TProfile3D.h"


//ClassImp(ReconstructionProfile6D);


ReconstructionProfile6D::ReconstructionProfile6D()
 : TNamed()
{
}


const ReconstructionProfile6D & ReconstructionProfile6D::operator=(const ReconstructionProfile6D &in)
{
  TNamed::operator=(in);
  x_axis_.Set(in.x_axis_.GetNbins(), in.x_axis_.GetXmin(), in.x_axis_.GetXmax());
  y_axis_.Set(in.y_axis_.GetNbins(), in.y_axis_.GetXmin(), in.y_axis_.GetXmax());
  z_axis_.Set(in.z_axis_.GetNbins(), in.z_axis_.GetXmin(), in.z_axis_.GetXmax());
  contents = in.contents;
  return in;
}


ReconstructionProfile6D::ReconstructionProfile6D(const char* name, const char* title, Int_t nbinsx, 
      Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup, 
      Int_t nbinsz, Double_t zlow, Double_t zup, Option_t* option)
 : TNamed(name, title), 
   x_axis_(nbinsx, xlow, xup),
   y_axis_(nbinsy, ylow, yup),
   z_axis_(nbinsz, zlow, zup)
{
  AllocateContents();
}


void ReconstructionProfile6D::AllocateContents()
{
  unsigned int xbins = x_axis_.GetNbins()+2;
  unsigned int ybins = y_axis_.GetNbins()+2;
  unsigned int zbins = z_axis_.GetNbins()+2;
  
  contents.resize(xbins);
  for(unsigned int x0=0; x0<xbins; ++x0)
  {
    contents[x0].resize(xbins);
    for(unsigned int x1=0; x1<xbins; ++x1)
    {
      contents[x0][x1].resize(ybins);
      for(unsigned int y0=0; y0<ybins; ++y0)
      {
        contents[x0][x1][y0].resize(ybins);
        for(unsigned int y1=0; y1<ybins; ++y1)
        {
          contents[x0][x1][y0][y1].resize(zbins);
          for(unsigned int z0=0; z0<zbins; ++z0)
          {
            contents[x0][x1][y0][y1][z0].resize(zbins);
          }
        }
      }
    }
  }
}


void ReconstructionProfile6D::Fill(const entry_point &p, Double_t t)
{
  int binx0,biny0,binz0;
  int binx1,biny1,binz1;
  binx0 = x_axis_.FindBin(p.x0);
  biny0 = y_axis_.FindBin(p.y0);
  binz0 = z_axis_.FindBin(p.z0);
  binx1 = x_axis_.FindBin(p.x1);
  biny1 = y_axis_.FindBin(p.y1);
  binz1 = z_axis_.FindBin(p.z1);
  
  contents[binx0][binx1][biny0][biny1][binz0][binz1].push_back(t);
}


void ReconstructionProfile6D::FillHistogram(TH1 &h, const entry_range &r)
{
  int x0bin_down, x0bin_up, x1bin_down, x1bin_up, y0bin_down, y0bin_up, 
      y1bin_down, y1bin_up, z0bin_down, z0bin_up, z1bin_down, z1bin_up;
  
  x0bin_down = x_axis_.FindBin(r.down.x0);
  x1bin_down = x_axis_.FindBin(r.down.x1);
  y0bin_down = y_axis_.FindBin(r.down.y0);
  y1bin_down = y_axis_.FindBin(r.down.y1);
  z0bin_down = z_axis_.FindBin(r.down.z0);
  z1bin_down = z_axis_.FindBin(r.down.z1);
  
  x0bin_up = x_axis_.FindBin(r.up.x0);
  x1bin_up = x_axis_.FindBin(r.up.x1);
  y0bin_up = y_axis_.FindBin(r.up.y0);
  y1bin_up = y_axis_.FindBin(r.up.y1);
  z0bin_up = z_axis_.FindBin(r.up.z0);
  z1bin_up = z_axis_.FindBin(r.up.z1);
  
  
  for(int x0_id=x0bin_down; x0_id<=x0bin_up; ++x0_id)
  {
    for(int x1_id=x1bin_down; x1_id<=x1bin_up; ++x1_id)
    {
      for(int y0_id=y0bin_down; y0_id<=y0bin_up; ++y0_id)
      {
        for(int y1_id=y1bin_down; y1_id<=y1bin_up; ++y1_id)
        {
          for(int z0_id=z0bin_down; z0_id<=z0bin_up; ++z0_id)
          {
            for(int z1_id=z1bin_down; z1_id<=z1bin_up; ++z1_id)
            {
              unsigned int size = contents[x0_id][x1_id][y0_id][y1_id][z0_id]
                    [z1_id].size();
              for(unsigned int content=0; content<size; ++content)
              {
                h.Fill(contents[x0_id][x1_id][y0_id][y1_id][z0_id][z1_id][content]);
              }
            }
          }
        }
      }
    }
  }
}


entry_range ReconstructionProfile6D::GetFullRange() const
{
  entry_range r;
  r.down.x0 = x_axis_.GetXmin();
  r.down.x1 = x_axis_.GetXmin();
  r.down.y0 = y_axis_.GetXmin();
  r.down.y1 = y_axis_.GetXmin();
  r.down.z0 = z_axis_.GetXmin();
  r.down.z1 = z_axis_.GetXmin();
  
  r.up.x0 = x_axis_.GetXmax();
  r.up.x1 = x_axis_.GetXmax();
  r.up.y0 = y_axis_.GetXmax();
  r.up.y1 = y_axis_.GetXmax();
  r.up.z0 = z_axis_.GetXmax();
  r.up.z1 = z_axis_.GetXmax();
  return r;
}


void ReconstructionProfile6D::FillHistogram(TH1 &h)
{
  FillHistogram(h, GetFullRange());
}


double ReconstructionProfile6D::Mean(const entry_range &r)
{
  int x0bin_down, x0bin_up, x1bin_down, x1bin_up, y0bin_down, y0bin_up, 
      y1bin_down, y1bin_up, z0bin_down, z0bin_up, z1bin_down, z1bin_up;
  
  x0bin_down = x_axis_.FindBin(r.down.x0);
  x1bin_down = x_axis_.FindBin(r.down.x1);
  y0bin_down = y_axis_.FindBin(r.down.y0);
  y1bin_down = y_axis_.FindBin(r.down.y1);
  z0bin_down = z_axis_.FindBin(r.down.z0);
  z1bin_down = z_axis_.FindBin(r.down.z1);
  
  x0bin_up = x_axis_.FindBin(r.up.x0);
  x1bin_up = x_axis_.FindBin(r.up.x1);
  y0bin_up = y_axis_.FindBin(r.up.y0);
  y1bin_up = y_axis_.FindBin(r.up.y1);
  z0bin_up = z_axis_.FindBin(r.up.z0);
  z1bin_up = z_axis_.FindBin(r.up.z1);
  
  double sum = 0.0;
  int count = 0;
  
  for(int x0_id=x0bin_down; x0_id<=x0bin_up; ++x0_id)
  {
    for(int x1_id=x1bin_down; x1_id<=x1bin_up; ++x1_id)
    {
      for(int y0_id=y0bin_down; y0_id<=y0bin_up; ++y0_id)
      {
        for(int y1_id=y1bin_down; y1_id<=y1bin_up; ++y1_id)
        {
          for(int z0_id=z0bin_down; z0_id<=z0bin_up; ++z0_id)
          {
            for(int z1_id=z1bin_down; z1_id<=z1bin_up; ++z1_id)
            {
              unsigned int size = contents[x0_id][x1_id][y0_id][y1_id][z0_id]
                    [z1_id].size();
              for(unsigned int content=0; content<size; ++content)
              {
                sum += contents[x0_id][x1_id][y0_id][y1_id][z0_id][z1_id][content];
                ++count;
              }
            }
          }
        }
      }
    }
  }
  return sum/count;
}


double ReconstructionProfile6D::Mean()
{
  return Mean(GetFullRange());
}


double ReconstructionProfile6D::Sigma(const entry_range &r)
{
  int x0bin_down, x0bin_up, x1bin_down, x1bin_up, y0bin_down, y0bin_up, 
      y1bin_down, y1bin_up, z0bin_down, z0bin_up, z1bin_down, z1bin_up;
  
  x0bin_down = x_axis_.FindBin(r.down.x0);
  x1bin_down = x_axis_.FindBin(r.down.x1);
  y0bin_down = y_axis_.FindBin(r.down.y0);
  y1bin_down = y_axis_.FindBin(r.down.y1);
  z0bin_down = z_axis_.FindBin(r.down.z0);
  z1bin_down = z_axis_.FindBin(r.down.z1);
  
  x0bin_up = x_axis_.FindBin(r.up.x0);
  x1bin_up = x_axis_.FindBin(r.up.x1);
  y0bin_up = y_axis_.FindBin(r.up.y0);
  y1bin_up = y_axis_.FindBin(r.up.y1);
  z0bin_up = z_axis_.FindBin(r.up.z0);
  z1bin_up = z_axis_.FindBin(r.up.z1);
  
  double sum = 0.0;
  double sum2 = 0.0;
  int count = 0;
  
  for(int x0_id=x0bin_down; x0_id<=x0bin_up; ++x0_id)
  {
    for(int x1_id=x1bin_down; x1_id<=x1bin_up; ++x1_id)
    {
      for(int y0_id=y0bin_down; y0_id<=y0bin_up; ++y0_id)
      {
        for(int y1_id=y1bin_down; y1_id<=y1bin_up; ++y1_id)
        {
          for(int z0_id=z0bin_down; z0_id<=z0bin_up; ++z0_id)
          {
            for(int z1_id=z1bin_down; z1_id<=z1bin_up; ++z1_id)
            {
              unsigned int size = contents[x0_id][x1_id][y0_id][y1_id][z0_id]
                    [z1_id].size();
              for(unsigned int content=0; content<size; ++content)
              {
                double val = contents[x0_id][x1_id][y0_id][y1_id][z0_id][z1_id][content];
                sum += val;
                sum2 += val*val;
                ++count;
              }
            }
          }
        }
      }
    }
  }
  double mean_x = sum/count;
  double mean_x2 = sum2/count;
  return TMath::Sqrt(mean_x2 - mean_x*mean_x);
}


double ReconstructionProfile6D::Sigma()
{
  return Sigma(GetFullRange());
}


int ReconstructionProfile6D::Entries(const entry_range &r)
{
  int x0bin_down, x0bin_up, x1bin_down, x1bin_up, y0bin_down, y0bin_up, 
      y1bin_down, y1bin_up, z0bin_down, z0bin_up, z1bin_down, z1bin_up;
  
  x0bin_down = x_axis_.FindBin(r.down.x0);
  x1bin_down = x_axis_.FindBin(r.down.x1);
  y0bin_down = y_axis_.FindBin(r.down.y0);
  y1bin_down = y_axis_.FindBin(r.down.y1);
  z0bin_down = z_axis_.FindBin(r.down.z0);
  z1bin_down = z_axis_.FindBin(r.down.z1);
  
  x0bin_up = x_axis_.FindBin(r.up.x0);
  x1bin_up = x_axis_.FindBin(r.up.x1);
  y0bin_up = y_axis_.FindBin(r.up.y0);
  y1bin_up = y_axis_.FindBin(r.up.y1);
  z0bin_up = z_axis_.FindBin(r.up.z0);
  z1bin_up = z_axis_.FindBin(r.up.z1);
  
  int count = 0;
  
  for(int x0_id=x0bin_down; x0_id<=x0bin_up; ++x0_id)
  {
    for(int x1_id=x1bin_down; x1_id<=x1bin_up; ++x1_id)
    {
      for(int y0_id=y0bin_down; y0_id<=y0bin_up; ++y0_id)
      {
        for(int y1_id=y1bin_down; y1_id<=y1bin_up; ++y1_id)
        {
          for(int z0_id=z0bin_down; z0_id<=z0bin_up; ++z0_id)
          {
            for(int z1_id=z1bin_down; z1_id<=z1bin_up; ++z1_id)
            {
              unsigned int size = contents[x0_id][x1_id][y0_id][y1_id][z0_id]
                    [z1_id].size();
              for(unsigned int content=0; content<size; ++content)
              {
                count+=contents[x0_id][x1_id][y0_id][y1_id][z0_id][z1_id].size();
              }
            }
          }
        }
      }
    }
  }
  
  return count;
}


int ReconstructionProfile6D::Entries()
{
  return Entries(GetFullRange());
}


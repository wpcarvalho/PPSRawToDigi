/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*  Jan Kašpar (jan.kaspar@gmail.com)
*
* $$RCSfile: ScoringPlaneDistributions.cc,v $: $
* $Revision: 9977 $
* $Date: 2015-01-12 15:00:26 +0100 (Mon, 12 Jan 2015) $
*
****************************************************************************/

#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"
#include "Geometry/TotemRPDetTopology/interface/RPTopology.h"

#include "TotemAlignment/RPRates/interface/ScoringPlaneDistributions.h"

#include "TH1D.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TDirectory.h"



ScoringPlaneDistributions::RP_profiles::RP_profiles(signed int Id, const char *label) : count(0)
{
  char buf[30];

  if (Id != 0) sprintf(buf, "%s: xProfile (%i)", label, Id);
  else sprintf(buf, "%s: xProfile", label);
  xProfile = new TH1D(buf, ";x   (mmm);hits/events total", 1000, -25., +25.);

  if (Id != 0) sprintf(buf, "%s: yProfile (%i)", label, Id);
  else sprintf(buf, "%s: yProfile", label);
  yProfile = new TH1D(buf, ";y   (mmm); hits/events total", 1000, -25., +25.);
}

//----------------------------------------------------------------------------------------------------

void ScoringPlaneDistributions::RP_profiles::Fill(double x, double y)
{
  count++;
  xProfile->Fill(x);
  yProfile->Fill(y);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

void ScoringPlaneDistributions::Init(const TotemRPGeometry *g, unsigned int N, double from, double to, signed int Id)
{
  geometry = g;

  eventCounter = 0;
  
  if (N == 0 || from >= to)
    throw cms::Exception("ScoringPlaneDistributions::Init") << "Wrong binning of offsets: N = " << N << ", from " << from << " to " << to;

  // build offset array
  double inc = (N > 1) ? (to - from) / (N - 1) : 0.;
  double o = from;
  offsets.clear();
  for (unsigned int i = 0; i < N; i++, o += inc)
    offsets.push_back(o); 

  // initialize members
  char buf[30];
  if (Id != 0) sprintf(buf, "hitDistribution (%i)", Id);
  else sprintf(buf, "hitDistribution");
  hitDistribution = new TGraph();
  hitDistribution->SetName(buf);
  hitDistribution->SetTitle("hit distribution;x   (mm);y   (mm)");

  // inialize RP profiles
  rpTop.clear();
  rpHor.clear();
  rpBot.clear();
  for (unsigned int i = 0; i < N; i++) {
    sprintf(buf, "rp_top, offset %.3f", offsets[i]);
    rpTop.push_back(RP_profiles(Id, buf)); 
    sprintf(buf, "rp_hor, offset %.3f", offsets[i]);
    rpHor.push_back(RP_profiles(Id, buf)); 
    sprintf(buf, "rp_bot, offset %.3f", offsets[i]);
    rpBot.push_back(RP_profiles(Id, buf)); 
  }

  // initialize overlap statistics
  for (unsigned int i = 0; i < N; i++) {
    vector<unsigned long> temp;
    for (unsigned int j = 0; j < N; j++)
      temp.push_back(0);
    overlapTopHor.push_back(temp); 
    overlapBotHor.push_back(temp); 
  }
}

//----------------------------------------------------------------------------------------------------

void ScoringPlaneDistributions::Fill(double x, double y)
{
  eventCounter++;

  hitDistribution->SetPoint(hitDistribution->GetN(), x, y);

  for (unsigned int i = 0; i < offsets.size(); i++) {
    double xp, yp, u, v;
    const double edgeToDetCenter = RPTopology::x_width_ / sqrt(2.) - RPTopology::phys_edge_lenght_ / 2.;

    // top
    xp = x; yp = y - edgeToDetCenter - offsets[i];
    u = (xp + yp) / sqrt(2); v = (yp - xp) / sqrt(2);
    bool hitTop = RPTopology::IsHit(u, v);
    if (hitTop) rpTop[i].Fill(xp, yp);

    // horizontal
    xp = x - edgeToDetCenter - offsets[i]; yp = y;
    u = (xp - yp) / sqrt(2); v = (xp + yp) / sqrt(2);
    if (RPTopology::IsHit(u, v)) rpHor[i].Fill(xp, yp);

    // bottom
    xp = x; yp = y + edgeToDetCenter + offsets[i];
    u = -(xp + yp) / sqrt(2); v = (xp - yp) / sqrt(2);
    bool hitBot = RPTopology::IsHit(u, v);
    if (hitBot) rpBot[i].Fill(xp, yp);

    for (unsigned int j = 0; j < offsets.size(); j++) {
      // retry horizontal with a different offset
      xp = x - edgeToDetCenter - offsets[j]; yp = y;
      u = (xp - yp) / sqrt(2); v = (xp + yp) / sqrt(2);
      bool hitHor = RPTopology::IsHit(u, v);

      if (hitTop && hitHor) overlapTopHor[i][j]++;
      if (hitTop && hitBot) overlapBotHor[i][j]++;
    }
  }
}

//----------------------------------------------------------------------------------------------------

void ScoringPlaneDistributions::WriteRPProfiles(const std::vector<RP_profiles> &prof)
{
  TGraphErrors *g = new TGraphErrors();
  g->SetTitle("relative rates;offset   (mm);hits / total events");
  for (unsigned int i = 0; i < prof.size(); i++) {
    int idx = g->GetN();
    g->SetPoint(idx, offsets[i], double(prof[i].count) / eventCounter);
    g->SetPointError(idx, 0., sqrt(prof[i].count) / eventCounter);
  }
  g->Write("rates");

  gDirectory = gDirectory->mkdir("offsets profiles");
  TDirectory *topDir = gDirectory;
  for (unsigned int i = 0; i < prof.size(); i++) {
    char buf[30];
    sprintf(buf, "offset %.3f", offsets[i]);
    gDirectory = topDir->mkdir(buf);
    prof[i].xProfile->Scale(1./eventCounter);
    prof[i].xProfile->Write();
    prof[i].yProfile->Scale(1./eventCounter);
    prof[i].yProfile->Write();
  }
}

//----------------------------------------------------------------------------------------------------

void ScoringPlaneDistributions::WriteOverlapProfiles(const overlap_statistics &os)
{
  TGraph2D *g = new TGraph2D();
  g->SetTitle("relative rates;horizontal offset   (mm);vertical offset   (mm);hits / total events");
  for (unsigned int i = 0; i < offsets.size(); i++) {
    for (unsigned int j = 0; j < offsets.size(); j++) {
      g->SetPoint(g->GetN(), offsets[j], offsets[i], double(os[i][j]) / eventCounter);
      //printf("* point: i = %i, j = %i, x = %f, y = %f, z = %f\n", i, j, offsets[j], offsets[i], double(os[i][j]) / eventCounter);
    }
  }
  g->Write("rates");
}

//----------------------------------------------------------------------------------------------------

void ScoringPlaneDistributions::Write()
{
  hitDistribution->Write();

  TDirectory *topDir = gDirectory;

  gDirectory = topDir->mkdir("rp_top");
  WriteRPProfiles(rpTop);

  gDirectory = topDir->mkdir("rp_horizontal");
  WriteRPProfiles(rpHor);

  gDirectory = topDir->mkdir("rp_bottom");
  WriteRPProfiles(rpBot);

  gDirectory = topDir->mkdir("overlap_top_horizontal");
  WriteOverlapProfiles(overlapTopHor);

  gDirectory = topDir->mkdir("overlap_bottom_horizontal");
  WriteOverlapProfiles(overlapBotHor);
}


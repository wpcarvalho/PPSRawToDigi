#include <iostream>
#include <cmath>
#include "TotemRPValidation/ValidationTools/interface/ReconstructionProfile.h"
#include "TPad.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"

double ksi_ratio_min = -5;
double ksi_ratio_max = 0;
double interpoint_ksi_ratio_dist = 0.5;

double Mmin = -2;
double Mmax = 5;
double interpoint_M_dist = 0.5;

double event_cut = 20;


double tmin = -5;
double tmax = 1;
double interpoint_t_dist = 0.5;

double ksimin = -4;
double ksimax = 0;
double interpoint_ksi_dist = 0.5;


void SingleArmEntries(ReconstructionProfile *p, TPad *pad)
{
  int t_points = (tmax-tmin)/interpoint_t_dist + 1;
  double int_t_delta = interpoint_t_dist/2.0;
  
  int ksi_points = (ksimax-ksimin)/interpoint_ksi_dist + 1;
  double int_ksi_delta = interpoint_ksi_dist/2.0;

  Double_t phi1 = 0;                                              // phi options
  Double_t phi2 = 2*TMath::Pi();

  TLegend *leg1 = new TLegend(.89,.98,.98,.65);
  
  pad->SetBorderMode(0);
  pad->SetFrameBorderMode(0);
  pad->SetLeftMargin(0.1);
  pad->SetRightMargin(0.05);
  
  TH1F *hr3 = pad->DrawFrame(ksimin-interpoint_ksi_dist,0,ksimax+interpoint_ksi_dist,5000);    
  hr3->SetXTitle("Log_{10}(-#xi)");
  hr3->SetYTitle("Events");
  hr3->GetXaxis()->SetTitleOffset(1);
  hr3->GetYaxis()->SetTitleOffset(1.3);
  
  TLegend *leg3 = new TLegend(.89,.98,.98,.65);
    
  
  int t_line_index = 0;
  
  double max_entry = -1e6;
  double min_entry = 1e6;
  
  double max_x_val = -1e6;
  double min_x_val = 1e6;
  
  for(int x=0; x<t_points; x++)
  {
    int entries_points = 0;
    double x_entries_points[1000];
    double y_entries_points[1000];
    
    double t_m = tmin+x*interpoint_t_dist;
    double t_low = t_m - int_t_delta;
    double t_hi = t_m + int_t_delta;
    for(int y=0; y<ksi_points; y++)
    {
      double ksi_m = ksimin + y*interpoint_ksi_dist;
      double ksi_low = ksi_m - int_ksi_delta;
      double ksi_hi = ksi_m + int_ksi_delta;
      x_entries_points[entries_points] = ksi_m;
      y_entries_points[entries_points] = 
        p->Entries(t_low, t_hi, ksi_low, ksi_hi, phi1, phi2);
      
      max_entry = (max_entry<y_entries_points[entries_points])?
          (y_entries_points[entries_points]):(max_entry);
      min_entry = (min_entry>y_entries_points[entries_points])?
          (y_entries_points[entries_points]):(min_entry);
            
      //FitParams par = p->FitStatisticalEstimators(t_low, t_hi, ksi_low, ksi_hi, phi1, phi2);
      int entries = y_entries_points[entries_points];
      if(entries>event_cut)
      {
        max_x_val = (max_x_val<x_entries_points[entries_points])?
            (x_entries_points[entries_points]):(max_x_val);
        min_x_val = (min_x_val>x_entries_points[entries_points])?
            (x_entries_points[entries_points]):(min_x_val);
        entries_points++;
      }
    }
      
    if(entries_points>0)
    {
      Char_t tvalue[6];                   // t value for legend
      sprintf(tvalue, "%4.1f", t_m);
      
      
      TGraph *gr3 = new TGraph(entries_points, x_entries_points, y_entries_points);
      leg1->AddEntry(gr3, tvalue, "p");     // legend
      leg1->SetHeader("Log(-t/GeV^{2})");
      gr3->SetMarkerColor(1);
      gr3->SetMarkerStyle(20+t_line_index);
      gr3->Draw("PL");
      
      Char_t tvalue1[6];                  // t value for legend
      sprintf(tvalue1, "%4.1f", t_m);
      leg1->Draw();
      t_line_index++;
    }
  }
  
  double d_ent = max_entry-min_entry;
  
  hr3->SetMaximum(max_entry+d_ent*0.1);
  hr3->SetMinimum(min_entry-d_ent*0.1);
  hr3->GetXaxis()->SetRangeUser(min_x_val-0.25, max_x_val+0.25);
}


void SingleArmSigma(ReconstructionProfile *p, char *y_label, TPad *pad, bool logscale=false)
{
  int t_points = (tmax-tmin)/interpoint_t_dist + 1;
  double int_t_delta = interpoint_t_dist/2.0;
  
  int ksi_points = (ksimax-ksimin)/interpoint_ksi_dist + 1;
  double int_ksi_delta = interpoint_ksi_dist/2.0;

  Double_t phi1 = 0;                                              // phi options
  Double_t phi2 = 2*TMath::Pi();
  
  pad->SetBorderMode(0);
  pad->SetFrameBorderMode(0);
  pad->SetLeftMargin(0.12);
  pad->SetRightMargin(0.05);
  TH1F *hr1 = pad->DrawFrame(ksimin-interpoint_ksi_dist,0,ksimax+interpoint_ksi_dist,0.02);    
  hr1->SetXTitle("Log_{10}(-#xi)");
  hr1->SetYTitle(y_label);
  hr1->GetXaxis()->SetTitleOffset(1);
  hr1->GetYaxis()->SetTitleOffset(1.6);
  TLegend *leg1 = new TLegend(.89,.98,.98,.65);  
  
  int t_line_index = 0;
  
  double max_y_val = -1e6;
  double min_y_val = 1e6;
  
  double max_x_val = -1e6;
  double min_x_val = 1e6;
  
  for(int x=0; x<t_points; x++)
  {
    int points = 0;
    double x_points[1000];
    double y_points[1000];
    
    double t_m = tmin+x*interpoint_t_dist;
    double t_low = t_m - int_t_delta;
    double t_hi = t_m + int_t_delta;
    for(int y=0; y<ksi_points; y++)
    {
      double ksi_m = ksimin + y*interpoint_ksi_dist;
      double ksi_low = ksi_m - int_ksi_delta;
      double ksi_hi = ksi_m + int_ksi_delta;

      int entries = p->Entries(t_low, t_hi, ksi_low, ksi_hi, phi1, phi2);
      if(entries>event_cut)
      {
        double val = p->Sigma(t_low, t_hi, ksi_low, ksi_hi, phi1, phi2);
        x_points[points] = ksi_m;
        y_points[points] = val;
        max_y_val = (max_y_val<val)?(val):(max_y_val);
        min_y_val = (min_y_val>val)?(val):(min_y_val);
        
        max_x_val = (max_x_val<x_points[points])?
            (x_points[points]):(max_x_val);
        min_x_val = (min_x_val>x_points[points])?
            (x_points[points]):(min_x_val);
        points++;
      }
    }
      
    if(points>0)
    {
      TGraph *gr1 = new TGraph(points, x_points, y_points);
      gr1->SetMarkerColor(1);
      gr1->SetMarkerStyle(20+t_line_index);
      gr1->Draw("PL");
      
      Char_t tvalue[6];                   // t value for legend
      sprintf(tvalue, "%4.1f", t_m);
      
      leg1->AddEntry(gr1, tvalue, "p");     // legend
      leg1->SetHeader("Log(-t/GeV^{2})");
      leg1->Draw();
      
      t_line_index++;
    }
  }
  
  double d_val = max_y_val-min_y_val;
  
  if(logscale)
  {
    hr1->SetMaximum(max_y_val*2);
    hr1->SetMinimum(min_y_val/2);
    pad->SetLogy();
  }
  else
  {
    hr1->SetMaximum(max_y_val+d_val*0.1);
    hr1->SetMinimum(min_y_val-d_val*0.1);
  }
  hr1->GetXaxis()->SetRangeUser(min_x_val-0.25, max_x_val+0.25);
}



void SingleArmMean(ReconstructionProfile *p, char *y_label, TPad *pad)
{
  int t_points = (tmax-tmin)/interpoint_t_dist + 1;
  double int_t_delta = interpoint_t_dist/2.0;
  
  int ksi_points = (ksimax-ksimin)/interpoint_ksi_dist + 1;
  double int_ksi_delta = interpoint_ksi_dist/2.0;

  Double_t phi1 = 0;                                              // phi options
  Double_t phi2 = 2*TMath::Pi();
  
  pad->SetBorderMode(0);
  pad->SetFrameBorderMode(0);
  pad->SetLeftMargin(0.12);
  pad->SetRightMargin(0.05);
  TH1F *hr1 = pad->DrawFrame(ksimin-interpoint_ksi_dist,0,ksimax+interpoint_ksi_dist,0.02);    
  hr1->SetXTitle("Log_{10}(-#xi)");
  hr1->SetYTitle(y_label);
  hr1->GetXaxis()->SetTitleOffset(1);
  hr1->GetYaxis()->SetTitleOffset(1.6);
  TLegend *leg1 = new TLegend(.89,.98,.98,.65);  
  
  int t_line_index = 0;
  
  double max_y_val = -1e6;
  double min_y_val = 1e6;
  
  double max_x_val = -1e6;
  double min_x_val = 1e6;
  
  for(int x=0; x<t_points; x++)
  {
    int points = 0;
    double x_points[1000];
    double y_points[1000];
    
    double t_m = tmin+x*interpoint_t_dist;
    double t_low = t_m - int_t_delta;
    double t_hi = t_m + int_t_delta;
    for(int y=0; y<ksi_points; y++)
    {
      double ksi_m = ksimin + y*interpoint_ksi_dist;
      double ksi_low = ksi_m - int_ksi_delta;
      double ksi_hi = ksi_m + int_ksi_delta;

      int entries = p->Entries(t_low, t_hi, ksi_low, ksi_hi, phi1, phi2);
      if(entries>event_cut)
      {
        double val = p->Mean(t_low, t_hi, ksi_low, ksi_hi, phi1, phi2);
        x_points[points] = ksi_m;
        y_points[points] = val;
        max_y_val = (max_y_val<val)?(val):(max_y_val);
        min_y_val = (min_y_val>val)?(val):(min_y_val);
        
        max_x_val = (max_x_val<x_points[points])?
            (x_points[points]):(max_x_val);
        min_x_val = (min_x_val>x_points[points])?
            (x_points[points]):(min_x_val);
        points++;
      }
    }
      
    if(points>0)
    {
      TGraph *gr1 = new TGraph(points, x_points, y_points);
      gr1->SetMarkerColor(1);
      gr1->SetMarkerStyle(20+t_line_index);
      gr1->Draw("PL");
      
      Char_t tvalue[6];                   // t value for legend
      sprintf(tvalue, "%4.1f", t_m);
      
      leg1->AddEntry(gr1, tvalue, "p");     // legend
      leg1->SetHeader("Log(-t/GeV^{2})");
      leg1->Draw();
      
      t_line_index++;
    }
  }
  
  double d_val = max_y_val-min_y_val;
  
  hr1->SetMaximum(max_y_val+d_val*0.1);
  hr1->SetMinimum(min_y_val-d_val*0.1);
  
  hr1->GetXaxis()->SetRangeUser(min_x_val-0.25, max_x_val+0.25);
}


void M_Entries(ReconstructionProfile *p, TPad *pad, bool logscale=false) 
{
  int ksi_ratio_points = (ksi_ratio_max-ksi_ratio_min)/interpoint_ksi_ratio_dist + 1;
  double int_ksi_ratio_delta = interpoint_ksi_ratio_dist/2.0;
  
  int M_points = (Mmax-Mmin)/interpoint_M_dist + 1;
  double int_M_delta = interpoint_M_dist/2.0;

  Double_t phi1 = 0;                                              // phi options
  Double_t phi2 = 2*TMath::Pi();
    
//  TCanvas *c3 = new TCanvas("c3","events (log #xi)",200,0,700,500);
  pad->SetBorderMode(0);
  pad->SetFrameBorderMode(0);
  pad->SetLeftMargin(0.1);
  pad->SetRightMargin(0.05);
  
  pad->cd();
  TH1F *hr3 = pad->DrawFrame(Mmin-interpoint_M_dist,0,Mmax+interpoint_M_dist,5000);    
  hr3->SetXTitle("Log_{10}(M/GeV)");
  hr3->SetYTitle("Events");
  hr3->GetXaxis()->SetTitleOffset(1);
  hr3->GetYaxis()->SetTitleOffset(1.3);
  
  TLegend *leg3 = new TLegend(.89,.98,.98,.65);
    
  int ksi_ratio_line_index = 0;
  
  double max_y_val = -1e6;
  double max_entry = -1e6;
  double min_y_val = 1e6;
  double min_entry = 1e6;
  
  double max_x_val = -1e6;
  double min_x_val = 1e6;
  
  for(int x=0; x<ksi_ratio_points; x++)
  {
    int points = 0;
    double x_points[1000];
    double y_points[1000];
    
    int entries_points = 0;
    double x_entries_points[1000];
    double y_entries_points[1000];
    
    double ksi_ratio_m = ksi_ratio_min+x*interpoint_ksi_ratio_dist;
    double ksi_ratio_low = ksi_ratio_m - int_ksi_ratio_delta;
    double ksi_ratio_hi = ksi_ratio_m + int_ksi_ratio_delta;
    for(int y=0; y<M_points; y++)
    {
      double M_m = Mmin + y*interpoint_M_dist;
      double M_low = M_m - int_M_delta;
      double M_hi = M_m + int_M_delta;
      
//      std::cout<<"M_low="<<M_low<<" M_hi="<<M_hi
//          <<" ksi_ratio_low="<<ksi_ratio_low<<" ksi_ratio_hi="<<ksi_ratio_hi<<std::endl;
//      std::cout<<"Entries="<<p->Entries(M_low, M_hi, ksi_ratio_low, ksi_ratio_hi, -10, 10)<<std::endl;
      
      x_entries_points[entries_points] = M_m;
      y_entries_points[entries_points] = 
        p->Entries(M_low, M_hi, ksi_ratio_low, ksi_ratio_hi, -10, 10);

      max_entry = (max_entry<y_entries_points[entries_points])?
          (y_entries_points[entries_points]):(max_entry);
      min_entry = (min_entry>y_entries_points[entries_points])?
          (y_entries_points[entries_points]):(min_entry);
          

      if(y_entries_points[entries_points]>event_cut)
      {
        max_x_val = (max_x_val<x_entries_points[entries_points])?
            (x_entries_points[entries_points]):(max_x_val);
        min_x_val = (min_x_val>x_entries_points[entries_points])?
            (x_entries_points[entries_points]):(min_x_val);
        entries_points++;
      }
    }
      
    if(entries_points>0)
    {
      Char_t tvalue[6];                   // t value for legend
      sprintf(tvalue, "%4.1f", ksi_ratio_m);
      
      TGraph *gr3 = new TGraph(entries_points, x_entries_points, y_entries_points);
      leg3->AddEntry(gr3, tvalue, "p");     // legend
      leg3->SetHeader("Log_{10}(#xi_{lo}/#xi_{hi})");
      gr3->SetMarkerColor(1);
      gr3->SetMarkerStyle(20+ksi_ratio_line_index);
      gr3->Draw("PL");
      
      Char_t tvalue1[6];                  // t value for legend
      sprintf(tvalue1, "%4.1f", ksi_ratio_m);
      leg3->Draw();
      ksi_ratio_line_index++;
    }
  }
  
  double d_ent = max_entry-min_entry;
  
  hr3->SetMaximum(max_entry+d_ent*0.0);
  hr3->SetMinimum(min_entry-d_ent*0.0);

  hr3->GetXaxis()->SetRangeUser(min_x_val-0.25, max_x_val+0.25);
}



void M_plot(ReconstructionProfile *p, char *y_label, TPad *pad, bool logscale=false) 
{
  int ksi_ratio_points = (ksi_ratio_max-ksi_ratio_min)/interpoint_ksi_ratio_dist + 1;
  double int_ksi_ratio_delta = interpoint_ksi_ratio_dist/2.0;
  
  int M_points = (Mmax-Mmin)/interpoint_M_dist + 1;
  double int_M_delta = interpoint_M_dist/2.0;

  Double_t phi1 = 0;                                              // phi options
  Double_t phi2 = 2*TMath::Pi();
  
  //TCanvas *c1 = new TCanvas("c1","#sigma(#xi)(log #xi)",200,0,700,500);
  pad->SetBorderMode(0);
  pad->SetFrameBorderMode(0);
  pad->SetLeftMargin(0.12);
  pad->SetRightMargin(0.05);
  
  pad->cd();
  TH1F *hr1 = pad->DrawFrame(Mmin-interpoint_M_dist,0,Mmax+interpoint_M_dist,0.02);    
  hr1->SetXTitle("Log_{10}(M/GeV)");
  hr1->SetYTitle(y_label);
  hr1->GetXaxis()->SetTitleOffset(1);
  hr1->GetYaxis()->SetTitleOffset(1.6);
  TLegend *leg1 = new TLegend(.89,.98,.98,.65);
    
  int ksi_ratio_line_index = 0;
  
  double max_y_val = -1e6;
  double max_entry = -1e6;
  double min_y_val = 1e6;
  double min_entry = 1e6;
  
  double max_x_val = -1e6;
  double min_x_val = 1e6;
  
  for(int x=0; x<ksi_ratio_points; x++)
  {
    int points = 0;
    double x_points[1000];
    double y_points[1000];
    
    double ksi_ratio_m = ksi_ratio_min+x*interpoint_ksi_ratio_dist;
    double ksi_ratio_low = ksi_ratio_m - int_ksi_ratio_delta;
    double ksi_ratio_hi = ksi_ratio_m + int_ksi_ratio_delta;
    for(int y=0; y<M_points; y++)
    {
      double M_m = Mmin + y*interpoint_M_dist;
      double M_low = M_m - int_M_delta;
      double M_hi = M_m + int_M_delta;

      if(p->Entries(M_low, M_hi, ksi_ratio_low, ksi_ratio_hi, -10, 10)>event_cut)
      {
        double val = p->Sigma(M_low, M_hi, ksi_ratio_low, ksi_ratio_hi, -10, 10);
        x_points[points] = M_m;
        y_points[points] = val;
        max_y_val = (max_y_val<val)?(val):(max_y_val);
        min_y_val = (min_y_val>val)?(val):(min_y_val);
        
        max_x_val = (max_x_val<x_points[points])?
            (x_points[points]):(max_x_val);
        min_x_val = (min_x_val>x_points[points])?
            (x_points[points]):(min_x_val);
        points++;
      }
    }
      
    if(points>0)
    {
      TGraph *gr1 = new TGraph(points, x_points, y_points);
      gr1->SetMarkerColor(1);
      gr1->SetMarkerStyle(20+ksi_ratio_line_index);
      gr1->Draw("PL");
      
      Char_t tvalue[6];                   // t value for legend
      sprintf(tvalue, "%4.1f", ksi_ratio_m);
      
      leg1->AddEntry(gr1, tvalue, "p");     // legend
      leg1->SetHeader("Log_{10}(#xi_{lo}/#xi_{hi})");
      leg1->Draw();
      
      ksi_ratio_line_index++;
    }
  }
  
  double d_val = max_y_val-min_y_val;

  if(logscale)
  {
    hr1->SetMaximum(max_y_val*2);
    hr1->SetMinimum(min_y_val/2);
    pad->SetLogy();
  }
  else
  {
    hr1->SetMaximum(max_y_val+d_val*0.1);
    hr1->SetMinimum(min_y_val-d_val*0.1);
  }
  
  hr1->GetXaxis()->SetRangeUser(min_x_val-0.25, max_x_val+0.25);
}


void M_Mean_plot(ReconstructionProfile *p, char *y_label, TPad *pad) 
{
  int ksi_ratio_points = (ksi_ratio_max-ksi_ratio_min)/interpoint_ksi_ratio_dist + 1;
  double int_ksi_ratio_delta = interpoint_ksi_ratio_dist/2.0;
  
  int M_points = (Mmax-Mmin)/interpoint_M_dist + 1;
  double int_M_delta = interpoint_M_dist/2.0;

  Double_t phi1 = 0;                                              // phi options
  Double_t phi2 = 2*TMath::Pi();
  
  //TCanvas *c1 = new TCanvas("c1","#sigma(#xi)(log #xi)",200,0,700,500);
  pad->SetBorderMode(0);
  pad->SetFrameBorderMode(0);
  pad->SetLeftMargin(0.12);
  pad->SetRightMargin(0.05);
  
  pad->cd();
  TH1F *hr1 = pad->DrawFrame(Mmin-interpoint_M_dist,0,Mmax+interpoint_M_dist,0.02);    
  hr1->SetXTitle("Log_{10}(M/GeV)");
  hr1->SetYTitle(y_label);
  hr1->GetXaxis()->SetTitleOffset(1);
  hr1->GetYaxis()->SetTitleOffset(1.6);
  TLegend *leg1 = new TLegend(.89,.98,.98,.65);
    
  int ksi_ratio_line_index = 0;
  
  double max_y_val = -1e6;
  double max_entry = -1e6;
  double min_y_val = 1e6;
  double min_entry = 1e6;
  
  double max_x_val = -1e6;
  double min_x_val = 1e6;
  
  for(int x=0; x<ksi_ratio_points; x++)
  {
    int points = 0;
    double x_points[1000];
    double y_points[1000];
    
    double ksi_ratio_m = ksi_ratio_min+x*interpoint_ksi_ratio_dist;
    double ksi_ratio_low = ksi_ratio_m - int_ksi_ratio_delta;
    double ksi_ratio_hi = ksi_ratio_m + int_ksi_ratio_delta;
    for(int y=0; y<M_points; y++)
    {
      double M_m = Mmin + y*interpoint_M_dist;
      double M_low = M_m - int_M_delta;
      double M_hi = M_m + int_M_delta;

      if(p->Entries(M_low, M_hi, ksi_ratio_low, ksi_ratio_hi, -10, 10)>event_cut)
      {
        double val = p->Mean(M_low, M_hi, ksi_ratio_low, ksi_ratio_hi, -10, 10);
        x_points[points] = M_m;
        y_points[points] = val;
        max_y_val = (max_y_val<val)?(val):(max_y_val);
        min_y_val = (min_y_val>val)?(val):(min_y_val);
        
        max_x_val = (max_x_val<x_points[points])?
            (x_points[points]):(max_x_val);
        min_x_val = (min_x_val>x_points[points])?
            (x_points[points]):(min_x_val);
        points++;
      }
    }
      
    if(points>0)
    {
      TGraph *gr1 = new TGraph(points, x_points, y_points);
      gr1->SetMarkerColor(1);
      gr1->SetMarkerStyle(20+ksi_ratio_line_index);
      gr1->Draw("PL");
      
      Char_t tvalue[6];                   // t value for legend
      sprintf(tvalue, "%4.1f", ksi_ratio_m);
      
      leg1->AddEntry(gr1, tvalue, "p");     // legend
      leg1->SetHeader("Log_{10}(#xi_{lo}/#xi_{hi})");
      leg1->Draw();
      
      ksi_ratio_line_index++;
    }
  }
  
  double d_val = max_y_val-min_y_val;
  
  hr1->SetMaximum(max_y_val+d_val*0.1);
  hr1->SetMinimum(min_y_val-d_val*0.1);
  hr1->GetXaxis()->SetRangeUser(min_x_val-0.25, max_x_val+0.25);
}


void DrawMPlots(void)
{
  ReconstructionProfile *ksi_error_right = (ReconstructionProfile *) gDirectory->Get("ksi_error_vs_log_t_log_ksi_phi_right_");
  TCanvas *c2_right = new TCanvas("Sigma(xi) right", "Sigma(xi)");
  SingleArmSigma(ksi_error_right, "#sigma(#xi)", c2_right);
  c2_right->SaveAs("01_xi_xi.gif");

  ReconstructionProfile *ksi_rel_error_right = (ReconstructionProfile *) gDirectory->Get("ksi_rel_error_vs_log_t_log_ksi_phi_right_");
  TCanvas *c2a_right = new TCanvas("Sigma(xi)/xi right", "Sigma(xi)/xi");
  SingleArmSigma(ksi_rel_error_right, "#sigma(#xi)/#xi", c2a_right, true);
  c2a_right->SaveAs("02_xi_rel_xi.gif");
  
  ReconstructionProfile *t_x_error_right = (ReconstructionProfile *) gDirectory->Get("t_x_error_vs_log_t_log_ksi_phi_right_");
  TCanvas *c3a_right = new TCanvas("Sigma(t_x) right", "Sigma(t_x)");
  SingleArmSigma(t_x_error_right, "#sigma(t_{x})", c3a_right, true);
  c3a_right->SaveAs("03_tx_xi.gif");
  
  ReconstructionProfile *t_x_rel_error_right = (ReconstructionProfile *) gDirectory->Get("t_x_rel_error_vs_log_t_log_ksi_phi_right_");
  TCanvas *c3b_right = new TCanvas("Sigma(t_x)/t_x right", "Sigma(t_x)/t_x");
  SingleArmSigma(t_x_rel_error_right, "#sigma(t_{x})/t_{x}", c3b_right, true);
  c3b_right->SaveAs("04_tx_rel_xi.gif");
  
  ReconstructionProfile *t_y_error_right = (ReconstructionProfile *) gDirectory->Get("t_y_error_vs_log_t_log_ksi_phi_right_");
  TCanvas *c3c_right = new TCanvas("Sigma(t_y) right", "Sigma(t_y)");
  SingleArmSigma(t_y_error_right, "#sigma(t_{y})", c3c_right);
  c3c_right->SaveAs("05_ty_xi.gif");
  
  ReconstructionProfile *t_y_rel_error_right = (ReconstructionProfile *) gDirectory->Get("t_y_rel_error_vs_log_t_log_ksi_phi_right_");
  TCanvas *c3d_right = new TCanvas("Sigma(t_y)/t_y right", "Sigma(t_y)/t_y");
  SingleArmSigma(t_y_rel_error_right, "#sigma(t_{y})/t_{y}", c3d_right, true);
  c3d_right->SaveAs("06_ty_rel_xi.gif");
  
  ReconstructionProfile *t_error_right = (ReconstructionProfile *) gDirectory->Get("t_error_vs_log_t_log_ksi_phi_right_");
  TCanvas *c3_right = new TCanvas("Sigma(t) right", "Sigma(t)");
  SingleArmSigma(t_error_right, "#sigma(t)", c3_right, true);
  c3_right->SaveAs("07_t_xi.gif");
  
  ReconstructionProfile *t_rel_error_right = (ReconstructionProfile *) gDirectory->Get("t_rel_error_vs_log_t_log_ksi_phi_right_");
  TCanvas *c4_right = new TCanvas("Sigma(t)/t right", "Sigma(t)/t");
  SingleArmSigma(t_rel_error_right, "#sigma(t)/t", c4_right);
  c4_right->SaveAs("08_t_rel_xi.gif");
  
  ReconstructionProfile *phi_error_right = (ReconstructionProfile *) gDirectory->Get("phi_error_vs_log_t_log_ksi_phi_right_");
  TCanvas *c5_right = new TCanvas("Sigma(phi) right", "Sigma(phi)");
  SingleArmSigma(phi_error_right, "#sigma(#phi) [rad]", c5_right, true);
  c5_right->SaveAs("09_phi_xi.gif");
  
  ReconstructionProfile *thetax_right = (ReconstructionProfile *) gDirectory->Get("thetax_error_vs_log_t_log_ksi_phi_right_");
  TCanvas *c8a_right = new TCanvas("thetax_error_vs_log_t_log_ksi_phi_1 right", "thetax_error_vs_log_t_log_ksi_phi_1");
  SingleArmSigma(thetax_right, "#sigma(#theta_{x})", c8a_right);
  c8a_right->SaveAs("10_thetax_xi.gif");
  
  ReconstructionProfile *thetay_right = (ReconstructionProfile *) gDirectory->Get("thetay_error_vs_log_t_log_ksi_phi_right_");
  TCanvas *c8b_right = new TCanvas("thetay_error_vs_log_t_log_ksi_phi_1 right", "thetay_error_vs_log_t_log_ksi_phi_1");
  SingleArmSigma(thetay_right, "#sigma(#theta_{y})", c8b_right);
  c8b_right->SaveAs("11_thetay_xi.gif");
  
  
  
  ReconstructionProfile *m_error = (ReconstructionProfile *) gDirectory->Get("M_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  TCanvas *c1 = new TCanvas("Entries", "Entries");
  M_Entries(m_error, c1);
  c1->SaveAs("12_entries.gif");
  
  TCanvas *c2 = new TCanvas("Sigma(M)", "Sigma(M)");
  M_plot(m_error, "#sigma(M) [GeV]", c2);
  c2->SaveAs("13_M.gif");

  ReconstructionProfile *m_rel_error = (ReconstructionProfile *) gDirectory->Get("M_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");  
  TCanvas *c2ab = new TCanvas("Sigma(M)/M", "Sigma(M)/M");
  M_plot(m_rel_error, "#sigma(M)/M", c2ab, true);
  c2ab->SaveAs("14_M_rel.gif");

  ReconstructionProfile *ksi_hi_error = (ReconstructionProfile *) gDirectory->Get("ksi_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");  
  TCanvas *c2a = new TCanvas("Sigma(xi_hi)", "Sigma(xi_hi)");
  M_plot(ksi_hi_error, "#sigma(#xi_{hi}) [GeV^{2}]", c2a);
  c2a->SaveAs("15_xi_hi.gif");

  ReconstructionProfile *ksi_hi_rel_error = (ReconstructionProfile *) gDirectory->Get("ksi_higher_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");  
  TCanvas *c2ac = new TCanvas("Sigma(xi_hi)/xi_hi", "Sigma(xi_hi)/xi_hi");
  M_plot(ksi_hi_rel_error, "#sigma(#xi_{hi})/#xi_{hi}", c2ac, true);
  c2ac->SaveAs("16_xi_hi_rel.gif");

  ReconstructionProfile *ksi_lo_error = (ReconstructionProfile *) gDirectory->Get("ksi_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");  
  TCanvas *c2b = new TCanvas("Sigma(xi_lo)", "Sigma(xi_lo)");
  M_plot(ksi_lo_error, "#sigma(#xi_{lo}) [GeV^{2}]", c2b);
  c2b->SaveAs("17_xi_lo.gif");

  ReconstructionProfile *ksi_lo_rel_error = (ReconstructionProfile *) gDirectory->Get("ksi_lower_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");  
  TCanvas *c2ba = new TCanvas("Sigma(xi_lo)/xi_lo", "Sigma(xi_lo)/xi_lo");
  M_plot(ksi_lo_rel_error, "#sigma(#xi_{lo})/#xi_{lo}", c2ba, true);
  c2ba->SaveAs("18_xi_lo_rel.gif");
  
  ReconstructionProfile *t_higher_error = (ReconstructionProfile *) gDirectory->Get("t_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");  
  TCanvas *c3a = new TCanvas("Sigma(t_hi)", "Sigma(t_hi)");
  M_plot(t_higher_error, "#sigma(t_{hi}) [GeV^{2}]", c3a);
  c3a->SaveAs("19_t_hi.gif");
  
  ReconstructionProfile *t_lower_error = (ReconstructionProfile *) gDirectory->Get("t_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");  
  TCanvas *c3 = new TCanvas("Sigma(t_lo)", "Sigma(t_lo)");
  M_plot(t_lower_error, "#sigma(t_{lo}) [GeV^{2}]", c3);
  c3->SaveAs("20_t_lo.gif");

  ReconstructionProfile *t_higher_rel_error = (ReconstructionProfile *) gDirectory->Get("t_higher_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");  
  TCanvas *c4 = new TCanvas("Sigma(t_hi)/t_hi", "Sigma(t_hi)/t_hi");
  M_plot(t_higher_rel_error, "#sigma(t_{hi})/t_{hi}", c4);
  c4->SaveAs("21_t_hi_rel.gif");

  ReconstructionProfile *t_lower_rel_error = (ReconstructionProfile *) gDirectory->Get("t_lower_rel_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");  
  TCanvas *c4a = new TCanvas("Sigma(t_lo)/t_lo", "Sigma(t_lo)/t_lo");
  M_plot(t_lower_rel_error, "#sigma(t_{lo})/t_{lo}", c4a);
  c4a->SaveAs("22_t_lo_rel.gif");

  ReconstructionProfile *phi_higher_error = (ReconstructionProfile *) gDirectory->Get("phi_higher_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");  
  TCanvas *c4b = new TCanvas("Sigma(phi_hi)", "Sigma(phi_hi)");
  M_plot(phi_higher_error, "#sigma(#phi_{hi}) [rad]", c4b);
  c4b->SaveAs("23_phi_hi.gif");

  ReconstructionProfile *phi_lower_error = (ReconstructionProfile *) gDirectory->Get("phi_lower_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");  
  TCanvas *c4c = new TCanvas("Sigma(phi_lo)", "Sigma(phi_lo)");
  M_plot(phi_lower_error, "#sigma(#phi_{lo}) [rad]", c4c);
  c4c->SaveAs("24_phi_lo.gif");

  ReconstructionProfile *x_error = (ReconstructionProfile *) gDirectory->Get("x_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");  
  TCanvas *c5 = new TCanvas("Sigma(x*)", "Sigma(x*)");
  M_plot(x_error, "#sigma(x*) [mm]", c5);
  c5->SaveAs("25_x.gif");

  ReconstructionProfile *y_error = (ReconstructionProfile *) gDirectory->Get("y_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");  
  TCanvas *c6 = new TCanvas("Sigma(y*)", "Sigma(y*)");
  M_plot(y_error, "#sigma(y*) [mm]", c6);
  c6->SaveAs("26_y.gif");

  ReconstructionProfile *z_error = (ReconstructionProfile *) gDirectory->Get("z_error_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");  
  TCanvas *c7 = new TCanvas("Sigma(z*)", "Sigma(z*)");
  M_plot(z_error, "#sigma(z*) [mm]", c7);
  c7->SaveAs("27_z.gif");

  ReconstructionProfile *chi = (ReconstructionProfile *) gDirectory->Get("chsqndf_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  TCanvas *c8 = new TCanvas("Sigma(chi2/ndf)", "Sigma(chi2/ndf)");
  M_Mean_plot(chi, "#chi^{2}/ndf", c8);
  c8->SaveAs("28_chi2_mean.gif");
  
  TH1D *prob = new TH1D("prob", "prob", 200, 0, 1);
  ReconstructionProfile *probe = (ReconstructionProfile *) gDirectory->Get("probe_vs_log_M_log_ksi_lo_div_ksi_hi_log_ksi_higher_");
  probe->FillHistogram(*prob);
  TCanvas *c9 = new TCanvas("Prob", "Prob");
  prob->Draw();
  c9->SaveAs("29_prob.gif");
  
  TH1D *chi2 = new TH1D("chi2/n", "chi2/n", 2000, 0, 10);
  ReconstructionProfile *chi2_prof = (ReconstructionProfile *) gDirectory->Get("chi2_error_overn_N_vs_log_t_log_ksi_phi_right_");
  chi2_prof->FillHistogram(*chi2);
  TCanvas *c10 = new TCanvas("chi2/n", "chi2/n");
  chi2->Draw();
  c10->SaveAs("30_chi2.gif");
}


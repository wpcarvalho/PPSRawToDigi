#include "TotemRPValidation/ValidationPlots/interface/Inelastic1ArmPlots.h"

Inelastic1ArmPlots::Inelastic1ArmPlots(){
  tmin = -5;
  tmax = 1;
  interpoint_t_dist = 0.5;

  ksimin = -4;
  ksimax = 0;
  interpoint_ksi_dist = 0.5;

  event_cut = 20;
}

Inelastic1ArmPlots::~Inelastic1ArmPlots(){
}

void Inelastic1ArmPlots::SingleArmEntries(ReconstructionProfile *p, TPad *pad)
{
  int t_points = (int)((tmax-tmin)/interpoint_t_dist + 1);
  double int_t_delta = interpoint_t_dist/2.0;
  
  int ksi_points = (int)((ksimax-ksimin)/interpoint_ksi_dist + 1);
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
  
//  TLegend *leg3 = new TLegend(.89,.98,.98,.65);
    
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
      int entries = (int)(y_entries_points[entries_points]);
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


void Inelastic1ArmPlots::SingleArmSigma(ReconstructionProfile *p, const char *y_label, TPad *pad, bool logscale)
{
  int t_points = (int)((tmax-tmin)/interpoint_t_dist + 1);
  double int_t_delta = interpoint_t_dist/2.0;
  
  int ksi_points = (int)((ksimax-ksimin)/interpoint_ksi_dist + 1);
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

void Inelastic1ArmPlots::SingleArmMean(ReconstructionProfile *p, char *y_label, TPad *pad)
{
  int t_points = (int)((tmax-tmin)/interpoint_t_dist + 1);
  double int_t_delta = interpoint_t_dist/2.0;
  
  int ksi_points = (int)((ksimax-ksimin)/interpoint_ksi_dist + 1);
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



std::vector< TCanvas* > Inelastic1ArmPlots::getPlots()
{
  std::vector< TCanvas* > l;

  ReconstructionProfile *ksi_error = (ReconstructionProfile *) gDirectory->Get("ksi_error_vs_log_t_log_ksi_phi_1");
  TCanvas *cc1 = new TCanvas("Entries", "Entries");
  SingleArmEntries(ksi_error, cc1);
  l.push_back(cc1);
  
  TCanvas *cc2 = new TCanvas("Sigma(xi)", "Sigma(xi)");
  SingleArmSigma(ksi_error, "#sigma(#xi)", cc2, true);
  l.push_back(cc2);

  ReconstructionProfile *ksi_ref_error = (ReconstructionProfile *) gDirectory->Get("ksi_rel_error_vs_log_t_log_ksi_phi_1");
  TCanvas *cc2a = new TCanvas("Sigma(xi)/xi", "Sigma(xi)/xi");
  SingleArmSigma(ksi_ref_error, "#sigma(#xi)/#xi", cc2a, true);
  l.push_back(cc2a);
  
  ReconstructionProfile *t_error = (ReconstructionProfile *) gDirectory->Get("t_error_vs_log_t_log_ksi_phi_1");
  TCanvas *cc3 = new TCanvas("Sigma(t)", "Sigma(t)");
  SingleArmSigma(t_error, "#sigma(t)", cc3, true);
  l.push_back(cc3);
  
  ReconstructionProfile *t_rel_error = (ReconstructionProfile *) gDirectory->Get("t_rel_error_vs_log_t_log_ksi_phi_1");
  TCanvas *cc4 = new TCanvas("Sigma(t)/t", "Sigma(t)/t");
  SingleArmSigma(t_rel_error, "#sigma(t)/t", cc4, false);
  l.push_back(cc4);
  
  ReconstructionProfile *t_x_error = (ReconstructionProfile *) gDirectory->Get("t_x_error_vs_log_t_log_ksi_phi_1");
  TCanvas *cc4a = new TCanvas("Sigma(t_x)", "Sigma(t_x)");
  SingleArmSigma(t_x_error, "#sigma(t_{x})", cc4a, true);
  l.push_back(cc4a);
  
  ReconstructionProfile *t_x_rel_error = (ReconstructionProfile *) gDirectory->Get("t_x_rel_error_vs_log_t_log_ksi_phi_1");
  TCanvas *cc4b = new TCanvas("Sigma(t_x)/t_x", "Sigma(t_x)/t_x");
  SingleArmSigma(t_x_rel_error, "#sigma(t_{x})/t_{x}", cc4b, true);
  l.push_back(cc4b);
  
  ReconstructionProfile *t_y_error = (ReconstructionProfile *) gDirectory->Get("t_y_error_vs_log_t_log_ksi_phi_1");
  TCanvas *cc4c = new TCanvas("Sigma(t_y)", "Sigma(t_y)");
  SingleArmSigma(t_y_error, "#sigma(t_{y})", cc4c);
  l.push_back(cc4c);
  
  ReconstructionProfile *t_y_rel_error = (ReconstructionProfile *) gDirectory->Get("t_y_rel_error_vs_log_t_log_ksi_phi_1");
  TCanvas *cc4d = new TCanvas("Sigma(t_y)/t_y", "Sigma(t_y)/t_y");
  SingleArmSigma(t_y_rel_error, "#sigma(t_{y})/t_{y}", cc4d, true);
  l.push_back(cc4d);
  
  ReconstructionProfile *phi_error = (ReconstructionProfile *) gDirectory->Get("phi_error_vs_log_t_log_ksi_phi_1");
  TCanvas *cc5 = new TCanvas("Sigma(phi)", "Sigma(phi)");
  SingleArmSigma(phi_error, "#sigma(#phi) [rad]", cc5, true);
  l.push_back(cc5);
  
  ReconstructionProfile *vx_error = (ReconstructionProfile *) gDirectory->Get("vx_error_vs_log_t_log_ksi_phi_1");
  TCanvas *cc6 = new TCanvas("Sigma(vx)", "Sigma(vx)");
  SingleArmSigma(vx_error, "#sigma(vx)", cc6, false);
  l.push_back(cc6);
  
  ReconstructionProfile *vy_error = (ReconstructionProfile *) gDirectory->Get("vy_error_vs_log_t_log_ksi_phi_1");
  TCanvas *cc7 = new TCanvas("Sigma(vy)", "Sigma(vy)");
  SingleArmSigma(vy_error, "#sigma(vy)", cc7, false);
  l.push_back(cc7);
  
  ReconstructionProfile *chi = (ReconstructionProfile *) gDirectory->Get("chi2_error_overn_N_vs_log_t_log_ksi_phi_1");
  TCanvas *cc8 = new TCanvas("chi2/ndf", "chi2/ndf");
  SingleArmSigma(chi, "#chi^{2}/ndf", cc8, false);
  l.push_back(cc8);
  
  ReconstructionProfile *thetax = (ReconstructionProfile *) gDirectory->Get("thetax_error_vs_log_t_log_ksi_phi_1");
  TCanvas *cc8a = new TCanvas("thetax_error_vs_log_t_log_ksi_phi_1", "thetax_error_vs_log_t_log_ksi_phi_1");
  SingleArmSigma(thetax, "#sigma(#theta_{x})", cc8a, true);
  l.push_back(cc8a);
  
  ReconstructionProfile *thetay = (ReconstructionProfile *) gDirectory->Get("thetay_error_vs_log_t_log_ksi_phi_1");
  TCanvas *cc8b = new TCanvas("thetay_error_vs_log_t_log_ksi_phi_1", "thetay_error_vs_log_t_log_ksi_phi_1");
  SingleArmSigma(thetay, "#sigma(#theta_{y})", cc8b, false);
  l.push_back(cc8b);
  
  TH1D *prob = new TH1D("prob", "prob", 200, 0, 1);
  ReconstructionProfile *probe = (ReconstructionProfile *) gDirectory->Get("prob_function_vs_log_t_log_ksi_phi_1");
  probe->FillHistogram(*prob);
  TCanvas *cc9 = new TCanvas("Prob", "Prob");
  l.push_back(cc9);
  
  TH1D *chi2 = new TH1D("chi2/n", "chi2/n", 2000, 0, 10);
  ReconstructionProfile *chi2_prof = (ReconstructionProfile *) gDirectory->Get("chi2_error_overn_N_vs_log_t_log_ksi_phi_1");
  chi2_prof->FillHistogram(*chi2);
  TCanvas *c10 = new TCanvas("chi2/n", "chi2/n");
  l.push_back(c10);
  
  TH2::AddDirectory(kFALSE);
  TCanvas *c11 = new TCanvas("Acceptance", "Acceptance");
  TH2F * accept = (TH2F*) gDirectory->Get("log_ksi_log_t_total_acceptance_1");
  accept->SetXTitle("Log_{10}(-t/GeV^{2})");
  accept->SetYTitle("Log_{10}(-#xi)");
  gStyle->SetPalette(1);
  accept->SetStats(0);
  accept->Draw("colz");
  l.push_back(c11);

  return l;

}

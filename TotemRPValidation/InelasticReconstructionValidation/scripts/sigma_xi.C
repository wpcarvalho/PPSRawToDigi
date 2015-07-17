void sigma_xi(ReconstructionProfile *p) {

   Double_t tmin = -1;                                             // t options
   Double_t tmax = 1;
   Int_t tbins = 5;
   Double_t tstep = fabs(tmin - tmax)/(tbins-1)/2;

   Double_t ximin = -3.5;                                          // xi options
   Double_t ximax = -1;
   Int_t xibins = 6;
   Double_t xistep = fabs(ximin - ximax)/(xibins-1)/2;

   Double_t phi1 = 0;                                              // phi options
   Double_t phi2 = 2*TMath::Pi();

// sigma_xi_log_xi ************************************

   TCanvas *c1 = new TCanvas("c1","#sigma(#xi)(log #xi)",200,0,700,500);

   c1->SetBorderMode(0);c1->SetFrameBorderMode(0);
   c1->SetLeftMargin(0.12);c1->SetRightMargin(0.05);   

   Double_t xiA[1000];        // xi array
   Double_t tA[1000];         // t array
   Double_t sigA[1000];       // sigma array
   Double_t evnA[1000];       // event entries - must not be integer because of TGraph:
                              // TGraph TGraph::TGraph(Int_t n,const Double_t* x,const Double_t* y);

   TH1F *hr1 = c1->DrawFrame(ximin-xistep,0,ximax+xistep,0.02);    
         hr1->SetXTitle("log #xi");hr1->SetYTitle("#sigma(#xi)");
         hr1->GetXaxis()->SetTitleOffset(1);hr1->GetYaxis()->SetTitleOffset(1.6);

   TLegend *leg1 = new TLegend(.89,.98,.98,.65);

   Int_t i;
   for (i=0;i<tbins;i++)
   {   tA[i] = tmin + 2*i*tstep;
       Int_t j;
       for (j=0;j<xibins;j++)
       {   xiA[j] = ximin + 2*j*xistep;
           cout <<p->Entries(tA[i]-tstep,tA[i]+tstep,xiA[j]-xistep,xiA[j]+xistep,phi1,phi2)<<", "
                <<i<<","<<j<<", "
                <<tA[i]-tstep << ", " << tA[i]+tstep << ", " 
                <<xiA[j]-xistep << ", " << xiA[j]+xistep << endl;
           sigA[j] = p->Sigma(tA[i]-tstep,tA[i]+tstep,xiA[j]-xistep,xiA[j]+xistep,phi1,phi2);
       }     
       TGraph *gr1 = new TGraph(xibins,xiA,sigA);
               gr1->SetMarkerColor(1);gr1->SetMarkerStyle(20+i);gr1->Draw("PL");

       Char_t tvalue[6];                   // t value for legend
       sprintf(tvalue, "%4.1f", tA[i]);

       leg1->AddEntry(gr1,tvalue,"p");     // legend
       leg1->SetHeader("t values");
       leg1->Draw();
   }

   parameters(phi1,phi2,tmin,tmax,ximin,ximax,tbins,tstep,xibins,xistep); // for title

// sigma_xi_log_xi_log_t ************************************

   TCanvas *c2 = new TCanvas("c2","#sigma(#xi)(log #xi,log t)",200,0,700,500);
            c2->SetBorderMode(0);c2->SetFrameBorderMode(0);
            c2->SetLeftMargin(0.15);c2->SetRightMargin(0.05);

   TProfile2D *h = new TProfile2D("#sigma(#xi)","",tbins,tmin-tstep,tmax+tstep,xibins,ximin-xistep,ximax+xistep,"s");   
               h->SetXTitle("log t");h->SetYTitle("log #xi");h->SetZTitle("#sigma(#xi)");
               h->GetXaxis()->SetTitleOffset(1.5);
               h->GetYaxis()->SetTitleOffset(1.7);
               h->GetZaxis()->SetTitleOffset(1.6);

   Int_t i=0;
   while (i < tbins)
   {     Int_t j=0;
         while (j < xibins)
      	 {     Double_t t = tmin + 2*i*tstep;
               Double_t xi = ximin + 2*j*xistep;
               Double_t sig = p->Sigma(t-tstep,t+tstep,xi-xistep,xi+xistep,phi1,phi2);
	       h->Fill(t,xi,sig);
	       j++;
         }     i++;
   }

   TPaveStats *ptstats = new TPaveStats(0.8,0.8,0.99,0.99,"bfNDC");
               ptstats->SetOptStat(1111);
               ptstats->Draw();

   h->GetListOfFunctions()->Add(ptstats);gStyle->SetPalette(1);h->Draw("LEGO2");

   cout<<"t bins:"<<tbins<<", xi bins:"<<xibins<< ", tstep:"<<tstep<<", xistep:"
       <<xistep<<endl<<"t min:"<<tmin<<", t max:"<<tmax
       <<", xi min:"<<ximin<<", xi max:"<<ximax<<endl;

   parameters(phi1,phi2,tmin,tmax,ximin,ximax,tbins,tstep,xibins,xistep); 

// entries_log_xi ************************************

   TCanvas *c3 = new TCanvas("c3","events (log #xi)",200,0,700,500);
   c3->SetBorderMode(0);c3->SetFrameBorderMode(0);
   c3->SetLeftMargin(0.1);c3->SetRightMargin(0.05);

   TH1F *hr3 = c3->DrawFrame(ximin-xistep,0,ximax+xistep,5000);    
         hr3->SetXTitle("log #xi");hr3->SetYTitle("events");
         hr3->GetXaxis()->SetTitleOffset(1);hr3->GetYaxis()->SetTitleOffset(1.3);

   TLegend *leg3 = new TLegend(.89,.98,.98,.65);

   Int_t i;
   for (i=0;i<tbins;i++)
   {   tA[i] = tmin + 2*i*tstep;
       Int_t j;
       for (j=0;j<xibins;j++)
       {   xiA[j] = ximin + 2*j*xistep;
           evnA[j] = p->Entries(tA[i]-tstep,tA[i]+tstep,xiA[j]-xistep,xiA[j]+xistep,phi1,phi2);
       }     
       TGraph *gr3 = new TGraph(xibins,xiA,evnA);
               gr3->SetMarkerColor(1);gr3->SetMarkerStyle(20+i);gr3->Draw("PL");

       Char_t tvalue[6];                  // t value for legend
       sprintf(tvalue, "%4.1f", tA[i]);

       leg1->Draw();
   }

   parameters(phi1,phi2,tmin,tmax,ximin,ximax,tbins,tstep,xibins,xistep);

}

// parameters in title ************************************

void parameters(phi1,phi2,tmin,tmax,ximin,ximax,tbins,tstep,xibins,xistep) {
     Char_t tminvalue[6];sprintf(tminvalue, "%1.2f", tmin);
     Char_t tmaxvalue[6];sprintf(tmaxvalue, "%1.2f", tmax);
     Char_t ximinvalue[6];sprintf(ximinvalue, "%1.2f", ximin);
     Char_t ximaxvalue[6];sprintf(ximaxvalue, "%1.2f", ximax);
     Char_t phi1value[6];sprintf(phi1value, "%1.4f", phi1);
     Char_t phi2value[6];sprintf(phi2value, "%1.4f", phi2);
     Char_t tbinsvalue[6];sprintf(tbinsvalue, "%2.0f", tbins);
     Char_t tstepvalue[6];sprintf(tstepvalue, "%2.3f", tstep);
     Char_t xibinsvalue[6];sprintf(xibinsvalue, "%2.0f", xibins);
     Char_t xistepvalue[6];sprintf(xistepvalue, "%2.3f", xistep);
     TPaveText *pt1 = new TPaveText(0.01,0.935,0.78,0.99,"NDC");pt1->SetBorderSize(0);
                pt1->AddText("#phi:              t:           xi:           tbins:      (#pm            )     xibins:       (#pm            )");
                pt1->Draw(); 
     TPaveText *pt2 = new TPaveText(0.07,0.93,0.12,0.99,"NDC");pt2->SetBorderSize(0);
                pt2->AddText(phi1value);pt2->AddText(phi2value);pt2->Draw();
     TPaveText *pt3 = new TPaveText(0.16,0.93,0.21,0.99,"NDC");pt3->SetBorderSize(0);
                pt3->AddText(tminvalue);pt3->AddText(tmaxvalue);pt3->Draw(); 
     TPaveText *pt4 = new TPaveText(0.25,0.93,0.30,0.99,"NDC");pt4->SetBorderSize(0);
                pt4->AddText(ximinvalue);pt4->AddText(ximaxvalue);pt4->Draw(); 
     TPaveText *pt5 = new TPaveText(0.375,0.94,0.405,0.99,"NDC");pt5->SetBorderSize(0);
                pt5->AddText(tbinsvalue);pt5->Draw(); 
     TPaveText *pt6 = new TPaveText(0.43,0.94,0.50,0.99,"NDC");pt6->SetBorderSize(0);
                pt6->AddText(tstepvalue);pt6->Draw(); 
     TPaveText *pt7 = new TPaveText(0.61,0.94,0.64,0.99,"NDC");pt7->SetBorderSize(0);
                pt7->AddText(xibinsvalue);pt7->Draw(); 
     TPaveText *pt8 = new TPaveText(0.67,0.94,0.74,0.99,"NDC");pt8->SetBorderSize(0);
                pt8->AddText(xistepvalue);pt8->Draw(); 
}

/****************************************************************************
*
* This is a part of the TOTEM offline software.
* Authors:
*   Frigyes Janos Nemes (frigyes.janos.nemes@cern.ch) 
*
* $$RCSfile: EdgeEfficiency.cc,v $: $
* $Revision: 1 $
* $Date: 2010-08-20 10:45:42       (Fri, 20 Aug 2010) $
*
****************************************************************************/

#include "TotemRPValidation/EdgeEfficiency/interface/EdgeEfficiency.h"

#include <iostream>
#include <fstream>

#include "DataFormats/Common/interface/DetSetVector.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPTrackCandidateCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "Geometry/TotemRecords/interface/RealGeometryRecord.h"
#include "Geometry/TotemRPDetTopology/interface/RPTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TFile.h"
#include "Math/SpecFunc.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TStyle.h"


using namespace std ;
using namespace edm ;

EdgeEfficiency::EdgeEfficiency(const edm::ParameterSet& ps) : 
	tolerance_in_v(ps.getParameter<double>("tolerance_in_v")),
	tolerance_in_angle(ps.getParameter<double>("tolerance_in_angle")),
	RP_angle_resolution(ps.getParameter<double>("RP_angle_resolution")),
	RP_size_along_the_z_axis(ps.getParameter<double>("RP_size_along_the_z_axis"))
{
	rPFittedTrackCollectionLabel = ps.getParameter<edm::InputTag>("RPFittedTrackCollectionLabel");
	detSetVectorRPRecoHitLabel = ps.getParameter<edm::InputTag>("DetSetVectorRPRecoHitLabel");
	Verbosity = ps.getParameter<unsigned int>("Verbosity");
	RootFileName = ps.getParameter<std::string>("RootFileName");
	selectedTest = ps.getParameter<unsigned int>("selectedTest");

	stripZeroPosition = RPTopology::last_strip_to_border_dist_ + (RPTopology::no_of_strips_-1)*RPTopology::pitch_ - RPTopology::y_width_/2.;

}


EdgeEfficiency::~EdgeEfficiency()
{
}


Double_t fitf(Double_t *x,Double_t *par) // This function (a Gaussian error function) will be fitted to the efficiency histogram. 
{
  Double_t arg = 0;
  if (par[2] != 0) arg = (x[0] - par[1])/par[2];
  Double_t fitval = (par[0]/2)*(1 + ROOT::Math::erf(arg));
  return fitval;
}


void EdgeEfficiency::beginRun(edm::Run const&, edm::EventSetup const& eSetup)
{
	eSetup.get<RealGeometryRecord>().get(geom); 			// The real geometrics is used.

        map<string,TH1D *> hists ; 					// This map will contain the "Accepted", "Expected", "Efficiency" strings as key values and the correspondig histrograms

	Double_t upper_edge_of_last_bin ; 			// [mm] 
	Double_t lower_edge_of_last_bin ; 			// [mm]

	if((selectedTest == 0) || (selectedTest == 1) || (selectedTest == 3) || (selectedTest == 4))
	{
		if(selectedTest == 4)
		{
			resolution_in_w = tolerance_in_v ; //RPTopology::pitch_ ;
		        upper_edge_of_last_bin = 700 ;                   // [mm]
                	lower_edge_of_last_bin = -700 ;                  // [mm]

		}
		else
		{
			resolution_in_w = tolerance_in_v ;
	                upper_edge_of_last_bin = 15 ;                   // [mm]
        	        lower_edge_of_last_bin = -30 ;                  // [mm]

		}
	}
	else if((selectedTest == 2) || (selectedTest == 5))
	{
		upper_edge_of_last_bin = 50 ;			// [mm] 
		lower_edge_of_last_bin = -50 ; 			// [mm]
		// resolution_in_w = RP_angle_resolution * RP_size_along_the_z_axis ;  // [mm] 
		resolution_in_w = 2 * RP_angle_resolution * 214 ;  // [mm] 
	}
	else {
		throw cms::Exception("RPReconstructedTracksValidation") << "selectedTest value should be between 0 and 5";
	}

	Int_t number_of_bins = (Int_t)((upper_edge_of_last_bin - lower_edge_of_last_bin ) / resolution_in_w) ;


	// Now we create three empty histograms for every detector : expected, accepted and efficiency
	TH1D::SetDefaultSumw2() ; // definition of the errors. THERE ARE HUGE ERRORS ?! TURNED OFF
	
	unsigned int start_arm ;

	if((selectedTest == 0) || (selectedTest == 1) || (selectedTest == 3) || (selectedTest == 4)) start_arm = 0 ;
	if((selectedTest == 2) || (selectedTest == 5)) start_arm = 0 ; // In a test beam geometry there is defined only one arm, it is arm=1.
	
	for(unsigned int arm = start_arm ; arm < 2 ; ++arm) 			// Loop over all of the arms.
	{
		const set<unsigned int> stations = geom->StationsInArm(arm);
		for (set<unsigned int>::const_iterator st = stations.begin(); st != stations.end(); ++st) 	// Loop over all the stations of the actual arm.
		{
        		const set<unsigned int> RPs = geom->RPsInStation(*st);
	       		for (set<unsigned int>::iterator rit = RPs.begin(); rit != RPs.end(); ++rit) 		// Loop over all the Roman Pots of the actual station.
        		{
        			const set<unsigned int> &dets = geom->DetsInRP(*rit);
	        		for (set<unsigned int>::iterator dit = dets.begin(); dit != dets.end(); ++dit) // Loop over all the detectors of the actual Roman Pot.
        			{
                        		unsigned int rawId = TotRPDetId::DecToRawId(*dit); // Get the rawId of the actual detector.

					char histname[40], histtitle[40] ;

					// Here we make three empty histograms, and we map them to the actual detector.
					sprintf(histname,"Expected"); sprintf(histtitle,"Detector ID: %i, Expected Hits",*dit);
					TH1D *expected = new TH1D(histname,histtitle,number_of_bins, lower_edge_of_last_bin, upper_edge_of_last_bin) ;

					if(selectedTest == 4) expected->SetXTitle("Strip number"); 
					else
						expected->SetXTitle("w [mm]"); 
					expected->SetYTitle("Number of hits");
					expected->SetLineColor(kGreen) ;
					expected->SetMarkerStyle(21);
					expected->SetMarkerColor(kBlack);
					expected->SetMarkerSize(0.15);


					sprintf(histname,"Accepted"); sprintf(histtitle,"Detector ID: %i, Accepted Hits",*dit);
					TH1D *accepted = new TH1D(histname,histtitle,number_of_bins, lower_edge_of_last_bin, upper_edge_of_last_bin) ;

					if(selectedTest == 4) accepted->SetXTitle("Strip number"); 
                                        else
                                                accepted->SetXTitle("w [mm]");

					accepted->SetYTitle("Number of hits");
					accepted->SetLineColor(kGreen) ;
                                        accepted->SetMarkerStyle(21);
                                        accepted->SetMarkerColor(kBlack);
					accepted->SetMarkerSize(0.15);


	                        	sprintf(histname,"Efficiency"); sprintf(histtitle,"Detector ID: %i, Efficiency",*dit);
	        	                TH1D *efficiency = new TH1D(histname,histtitle,number_of_bins, lower_edge_of_last_bin, upper_edge_of_last_bin) ;

					if(selectedTest == 4) efficiency->SetXTitle("Strip number");        
                                        else
                                                efficiency->SetXTitle("w [mm]");
					efficiency->SetYTitle("Efficiency");
					efficiency->SetLineColor(kGreen) ;
                                        efficiency->SetMarkerStyle(21);
                                        efficiency->SetMarkerColor(kBlack);
					efficiency->SetMarkerSize(0.15);


					hists.clear() ;

		        	        hists.insert(pair<string,TH1D*>("Expected",expected)) ;
		                	hists.insert(pair<string,TH1D*>("Accepted",accepted)) ;
		        	        hists.insert(pair<string,TH1D*>("Efficiency",efficiency)) ;

					// The key of the map contains the RPId too.
					// I only use this RPId information when I write out the histograms into the root file into directories, where the name of the directory is the RPId.
					pair <unsigned int,unsigned int> key (*rit,rawId);

					// Inserting the three type of histograms for the actual detector into the histogram map.
	        	                histograms.insert(pair< pair<unsigned int, unsigned int>, map<string,TH1D *> >(key,hists)) ;
        			}	
			}	
		}
	}

//	In the following lines we enumerate the reference detector pairs
	pair<unsigned int, unsigned int> reference_detectors ;
	
	reference_detectors.first = 22 ; reference_detectors.second = 23 ;
	map_from_tested_RP_to_its_reference_RPs.insert(pair<unsigned int, pair<unsigned int, unsigned int> >(20, reference_detectors)) ;
        map_from_tested_RP_to_its_reference_RPs.insert(pair<unsigned int, pair<unsigned int, unsigned int> >(21, reference_detectors)) ;
        map_from_tested_RP_to_its_reference_RPs.insert(pair<unsigned int, pair<unsigned int, unsigned int> >(24, reference_detectors)) ;
        map_from_tested_RP_to_its_reference_RPs.insert(pair<unsigned int, pair<unsigned int, unsigned int> >(25, reference_detectors)) ;

//      The RP with Id 22 is a horizontal one, so it has two reference detector pairs (20,24) and (21,25). We use the former, and comment the latter. 
        reference_detectors.first = 20 ; reference_detectors.second = 24 ;
        map_from_tested_RP_to_its_reference_RPs.insert(pair<unsigned int, pair<unsigned int, unsigned int> >(22, reference_detectors)) ;
//        reference_detectors.first = 21 ; reference_detectors.second = 25 ;
//        map_from_tested_RP_to_its_reference_RPs.insert(pair<unsigned int, pair<unsigned int, unsigned int> >(22, reference_detectors)) ;

//      The RP with Id 23 is a horizontal one, so it has two reference detector pairs (20,24) and (21,25). We use the former, and comment the latter. 
        reference_detectors.first = 20 ; reference_detectors.second = 24 ;
        map_from_tested_RP_to_its_reference_RPs.insert(pair<unsigned int, pair<unsigned int, unsigned int> >(23, reference_detectors)) ;
//        reference_detectors.first = 21 ; reference_detectors.second = 25 ;
//        map_from_tested_RP_to_its_reference_RPs.insert(pair<unsigned int, pair<unsigned int, unsigned int> >(23, reference_detectors)) ;

        reference_detectors.first = 122 ; reference_detectors.second = 123 ;
        map_from_tested_RP_to_its_reference_RPs.insert(pair<unsigned int, pair<unsigned int, unsigned int> >(120, reference_detectors)) ;
        map_from_tested_RP_to_its_reference_RPs.insert(pair<unsigned int, pair<unsigned int, unsigned int> >(121, reference_detectors)) ;
        map_from_tested_RP_to_its_reference_RPs.insert(pair<unsigned int, pair<unsigned int, unsigned int> >(124, reference_detectors)) ;
        map_from_tested_RP_to_its_reference_RPs.insert(pair<unsigned int, pair<unsigned int, unsigned int> >(125, reference_detectors)) ;

//      The RP with Id 122 is a horizontal one, so it has two reference detector pairs (120,124) and (121,125). We use the former, and comment the latter. 
        reference_detectors.first = 120 ; reference_detectors.second = 124 ;
        map_from_tested_RP_to_its_reference_RPs.insert(pair<unsigned int, pair<unsigned int, unsigned int> >(122, reference_detectors)) ;
//        reference_detectors.first = 121 ; reference_detectors.second = 125 ;
//        .insert(pair<unsigned int, pair<unsigned int, unsigned int> >(122, reference_detectors)) ;

//      The RP with Id 123 is a horizontal one, so it has two reference detector pairs (120,124) and (121,125). We use the former, and comment the latter. 
        reference_detectors.first = 120 ; reference_detectors.second = 124 ;
        map_from_tested_RP_to_its_reference_RPs.insert(pair<unsigned int, pair<unsigned int, unsigned int> >(123, reference_detectors)) ;
//        reference_detectors.first = 121 ; reference_detectors.second = 125 ;
//        .insert(pair<unsigned int, pair<unsigned int, unsigned int> >(123, reference_detectors)) ;

}

void EdgeEfficiency::endJob()
{
        TFile *f = new TFile(RootFileName.c_str(), "recreate") ; // The output file.
	f->cd("/") ;

	unsigned int prevRP = UINT_MAX ; 		// The prevRP number variable contains the Id of the previous RP in the loop below. This variable has to start from a NON_EXISTING !!! RP number.
	unsigned int prevDet = UINT_MAX ; 		// The same as above for the detectors.

//      string filenameTeX("/afs/cern.ch/user/f/fnemes/scratch/logs/detectors/") ;

	ostringstream resolution_str, tolerance_str ;
	resolution_str << resolution_in_w ;
	tolerance_str << tolerance_in_v ;


//      filenameTeX += "Edge_study_Bin_size_" ; filenameTeX += resolution_str.str() ; filenameTeX += "_mm_Tolerance_" ; filenameTeX += tolerance_str.str() ; filenameTeX += "_mm.tex" ;

//  	ofstream TeXfile ;
//	TeXfile.open(filenameTeX.c_str()) ;

//	TeXfile << "\\documentclass[10pt,landscape]{article}" << endl << "\\usepackage{graphicx}" << endl << "\\usepackage{float}" << endl << "\\usepackage{color}" << endl << 
//		"\\usepackage[top=0.2in, bottom=0.2in, left=0.3in, right=0.3in]{geometry}"  << endl << "\\usepackage{amsmath,amssymb}" << endl << endl << "\\begin{document}" << endl ;

	// Writing out the histograms from the histogram map to the root file.
        for(map<pair<unsigned int, unsigned int>, map<string,TH1D *> >::iterator mapit = histograms.begin() ; mapit != histograms.end() ; ++mapit) // Loop over the key(RPId,rawId) ordered map. 
        {
		// Creating a subdirectory in the root file for the given RP, if it is not exist yet.
		char directory_name[40] ;
		sprintf(directory_name,"RP%i",mapit->first.first) ;
		if(mapit->first.first != prevRP) 	// If the previous RP is not equal with the actual one, a new directory has to be made for the actual RP.
		{
			f->mkdir(directory_name) ;
			f->cd(directory_name) ;
			prevRP = mapit->first.first ;
		}
		else
			f->cd(directory_name) ;

		// Creating a subdirectory in the root file for the given silicon detector, if it is not exist yet. 
		char detector_name[40] ;
		unsigned int detId = TotRPDetId::RawToDecId(mapit->first.second);
		sprintf(detector_name,"%i",detId) ;

                if(mapit->first.second != prevDet) 	// If the previous Det not the same than the actual one, a new directory has to be made for the detector.
                {
                        gDirectory->mkdir(detector_name) ;
                        gDirectory->cd(detector_name) ;
                        prevDet = mapit->first.second ;
                }
                else
                        gDirectory->cd(detector_name) ;
			
		// Creating the canvas, drawing the expected and accepted histograms. 
		char canvas_name[60] ;
		sprintf(canvas_name,"The efficiency near the edge of the %i silicon detector",detId);

		TCanvas *c1 = new TCanvas(canvas_name,"Efficiency near the edge of the silicon detector",10,10,1400,1000);
                c1->SetFillColor(18);

   		
		TPad *pad1 = new TPad("pad1","The pad with the function",0.02,0.50,0.48,0.98,21);
		TPad *pad2 = new TPad("pad2","The pad with the histogram",0.02,0.02,0.48,0.48,21);
		TPad *pad3 = new TPad("pad3","The pad with the histogram",0.50,0.50,0.98,0.98,21);
		TPad *pad4 = new TPad("pad4","The pad with the histogram",0.50,0.02,0.98,0.48,21);

		pad1->GetFrame()->SetFillColor(42);
		pad2->GetFrame()->SetFillColor(42);
		pad3->GetFrame()->SetFillColor(42);
		pad4->GetFrame()->SetFillColor(42);

		pad1->GetFrame()->SetBorderMode(-1);
		pad2->GetFrame()->SetBorderMode(-1);
		pad3->GetFrame()->SetBorderMode(-1);
		pad4->GetFrame()->SetBorderMode(-1);

		pad1->GetFrame()->SetBorderSize(5);
		pad2->GetFrame()->SetBorderSize(5);
		pad3->GetFrame()->SetBorderSize(5);
		pad4->GetFrame()->SetBorderSize(5);

		pad1->Draw() ; pad2->Draw(); pad3->Draw() ; pad4->Draw() ;
		pad1->cd() ; 

              	if(Verbosity && ((mapit->second["Expected"]->GetEntries())>0)) cout << "Expected non-empty histogram RPiD and DetId: " << mapit->first.first << " " <<  TotRPDetId::RawToDecId(mapit->first.second) << endl ;
	      	mapit->second["Expected"]->Draw("E1") ;
                mapit->second["Expected"]->Write() ;
		pad2->cd() ;		
	      	if(Verbosity && ((mapit->second["Accepted"]->GetEntries())>0)) cout << "Accepted non-empty histogram RPiD and DetId: " << mapit->first.first << " " << TotRPDetId::RawToDecId(mapit->first.second) << endl ;
         	mapit->second["Accepted"]->Draw("E1") ;
         	mapit->second["Accepted"]->Write() ;

		pad3->cd() ;
		char histname[40] ;
		char histtitle[40] ;
			
		// Calculating the efficiency by a divison 
		mapit->second["Efficiency"]->Divide(mapit->second["Accepted"],mapit->second["Expected"],1,1,"B") ;  // The 3rd argument of the Divide is 100 in order to get the result in percent.
		mapit->second["Efficiency"]->Scale(100) ;

	        sprintf(histname,"Efficiency"); sprintf(histtitle,"Detector ID: %i, Efficiency",TotRPDetId::RawToDecId(mapit->first.second));

                mapit->second["Efficiency"]->SetName(histname) ;
                mapit->second["Efficiency"]->SetTitle(histtitle) ;
                mapit->second["Efficiency"]->SetYTitle("Efficiency [%]") ;

		/*
		// Fitting the efficiency histogram
		
		TF1 *func = new TF1("fit",fitf,-3,3,3) ;
		func->SetParameters(500,mapit->second["Efficiency"]->GetMean(),mapit->second["Efficiency"]->GetRMS());
		func->SetLineColor(kRed) ;
		func->SetLineWidth(3) ;
		func->SetLineStyle(2) ;
		
		*/
		// Give the parameters meaningful names.
		//func->SetParNames ("Constant","Mean_value","Sigma");
		// call TH1::Fit with the name of the TF1 object

		//mapit->second["Efficiency"]->Fit("fit","Q");
		//Double_t par[3];
		//func->GetParameters(par);

		mapit->second["Efficiency"]->SetAxisRange(0,110,"Y") ;
      		mapit->second["Efficiency"]->Draw("E1") ;
      		mapit->second["Efficiency"]->Write() ;

		CLHEP::Hep3Vector ep = geom->GetDetEdgePosition(mapit->first.second);
		CLHEP::Hep3Vector center = geom->LocalToGlobal(mapit->first.second, CLHEP::Hep3Vector(0, 0, 0));

		double distanceEdgeCenter = (ep - center).mag() ;
		TLine geometrical_edge_marker(-distanceEdgeCenter,0,-distanceEdgeCenter,100);

		if(selectedTest != 4)
		{	
        	        geometrical_edge_marker.SetLineWidth(1);
        	        geometrical_edge_marker.SetLineColor(kBlue);
	                geometrical_edge_marker.SetLineStyle(kDashed);
	                geometrical_edge_marker.Draw();
		}

		pad4->cd() ;

		TH1D *zoomedefficiency = (TH1D *)mapit->second["Efficiency"]->Clone() ;
	
		zoomedefficiency->SetAxisRange(-distanceEdgeCenter - 0.03,-distanceEdgeCenter + 0.27,"X") ;
		zoomedefficiency->SetAxisRange(0,110,"Y") ;
		zoomedefficiency->SetStats(kFALSE) ;

		zoomedefficiency->Draw("E1") ;
		zoomedefficiency->Write() ;
	        TLine geometrical_edge_marker2(-distanceEdgeCenter,0,-distanceEdgeCenter,100);

		if(selectedTest != 4)
		{
		        geometrical_edge_marker2.SetLineWidth(1);
        		geometrical_edge_marker2.SetLineColor(kBlue);
	                geometrical_edge_marker2.SetLineStyle(kDashed);
		        geometrical_edge_marker2.Draw();
		}

		//PaveText *text=new TPaveText(0.1,0.15,0.68,0.65);
		//char parameter_string[40] ;

		//text->AddText("The fit parameters of the Gaussian error function") ;
		//sprintf(parameter_string,"eta_0   = %f", par[0]); text->AddText(parameter_string) ;
		//sprintf(parameter_string,"w0      = %f", par[1]); text->AddText(parameter_string) ;
		//sprintf(parameter_string,"sigma_w = %f", par[2]); text->AddText(parameter_string) ;

		//text->Draw();
	

		// Everything is ready, writing out the canvas...
		c1->Write() ;

		/* Test. Remove from here...*/
/*
                string filenamePNG("/afs/cern.ch/user/f/fnemes/scratch/logs/detectors/") ;
		filenamePNG += detector_name ;

		string work ;
		work = resolution_str.str() ;
		
		size_t found;
  		found=work.find_first_of(".");
		while (found!=string::npos)
  		{
		    work[found]='0';
		    found=work.find_first_of(".",found+1);
		}
		filenamePNG += work ;

                work=tolerance_str.str() ;
                found=work.find_first_of(".");
                while (found!=string::npos)
                {
                    work[found]='0';
                    found=work.find_first_of(".",found+1);
                }
                filenamePNG += work ;

		filenamePNG += ".png" ; 

		TeXfile << "\\begin{figure}[H]" << endl << "\\centering" << endl << "\\includegraphics[scale=0.57]{" << filenamePNG << "}" << endl << "\\end{figure}" << endl ;

		c1->SaveAs(filenamePNG.c_str()) ;*/

		/* ... to here.*/
		f->cd() ;
	}

//	TeXfile << "\\end{document}" ;


//	TeXfile.close();
        f->Close() ;
	delete f ;

}

// In the test function below we test ALL of the detectors of 1 RomanPot. The Id of the Roman Pot is the second argumentum of the test(...) function.
void EdgeEfficiency::test(const edm::Event& e, unsigned int Id_of_the_tested_RP)
{
        edm::Handle< RPFittedTrackCollection > tracks_of_reference ; // RPFittedTrackCollection is part of the input. We will use only the point information X0,Y0,Z0.
        e.getByLabel(rPFittedTrackCollectionLabel, tracks_of_reference); 

        if (!tracks_of_reference.isValid())
                throw cms::Exception("RPReconstructedTracksValidation") << "edm::Handle< RPFittedTrackCollection > is invalid";


        RPFittedTrackCollection::const_iterator it1, it2 ;

	it1 = tracks_of_reference->find(map_from_tested_RP_to_its_reference_RPs[Id_of_the_tested_RP].first) ; // Determining the first and second reference detector  
	it2 = tracks_of_reference->find(map_from_tested_RP_to_its_reference_RPs[Id_of_the_tested_RP].second); // of the given tested detector. 

	if((it1 !=  tracks_of_reference->end()) && (it2 != tracks_of_reference->end()))
	{
		const set<unsigned int> &dets = geom->DetsInRP(Id_of_the_tested_RP);
                for (set<unsigned int>::iterator dit = dets.begin(); dit != dets.end(); ++dit) // Loop over the detectors of the tested RomanPot.
                {
                	// Calculating the intersecting point of the track and the detector plane in global coordinates.
                        unsigned int rawId = TotRPDetId::DecToRawId(*dit);
                        DetGeomDesc *ginfo = geom->GetDetector(rawId);			

			// Here we determine the line which connects the following two points: one point on the first, the another point is on the second reference detector.
                        double dx = (it2->second).X0() - (it1->second).X0() ;
                        double dy = (it2->second).Y0() - (it1->second).Y0() ;
			double dz = (it2->second).Z0() - (it1->second).Z0() ;

			double nx = dx / dz ; double ny = dy / dz ;
			if((nx > tolerance_in_angle) || (ny > tolerance_in_angle)) continue ; // If we connect the two point (on the different detectors) under too "large" angle, we drop this line
			double z = ginfo->translation().z() ;
                        double x = (it1->second).X0() + nx * (z - (it1->second).Z0()) ;
                        double y = (it1->second).Y0() + ny * (z - (it1->second).Z0()) ;

                        // Converting the common point into local coordinates
                        CLHEP::Hep3Vector hit(x, y, z);
                        hit = geom->GlobalToLocal(rawId, hit);
                        double u = hit.x(); double v = hit.y();


                        // Was this hit within the detector ?
                        if(RPTopology::IsHit(u, v)) 
			{
                        // Calculating the w coordinate of the hit and fill the histograms
			pair<unsigned int, unsigned int> key (Id_of_the_tested_RP,rawId) ;

                        double w = (u+v) / sqrt(2.) ; 			// Calculating the w coordinate...
			
			histograms[key]["Expected"]->Fill(w,1) ; 	// ...and filling out the histogram. 

		        Handle< edm::DetSetVector<RPRecoHit> > allHits ; // Now we get the hit information for the actual detector.
        		e.getByLabel(detSetVectorRPRecoHitLabel, allHits) ;

			// Here we loop over all detected (reconstructed) hits, and try to find at least one among them which is "close" to the actual expected hit.
			// This part of the algorithm would be faster.
			bool accepted = false ;	 // I use this variable to break the nested 2 loop if I accept the actual expected hit
                        for (DetSetVector<RPRecoHit>::const_iterator dit2 = allHits->begin(); (dit2 != allHits->end()) && !accepted ; ++dit2) 	// Loop over all DetSets.
                        	for(DetSet<RPRecoHit>::const_iterator hit2 = dit2->begin(); (hit2 != dit2->end()) && !accepted ; ++hit2)   	// Loop over all hits of one DetSet.
                                {
                                	if(hit2->DetId() != rawId) continue ; // If the given hit NOT belongs to the actual detector just continue.

                                       	if((v - tolerance_in_v) < (hit2->Position()) && ((hit2->Position()) < (v + tolerance_in_v))) // If the detected hit close enough to the exptected one, then...
					{
						histograms[key]["Accepted"]->Fill(w,1) ; 	// ...it is accepted.
                                        	accepted = true ; 					// I accept the actual expected hit just once time.
					}
                                        
                                }
			}
		}
	}
}

void EdgeEfficiency::testGivenPairOfRPs(int testedRP, int referenceRP, const edm::Event& e)
{
        edm::Handle< RPFittedTrackCollection > tracks_of_reference ; // RPFittedTrackCollection is part of the input. We will use only the point information X0,Y0,Z0.
        e.getByLabel(rPFittedTrackCollectionLabel, tracks_of_reference);

        if (!tracks_of_reference.isValid())
                throw cms::Exception("RPReconstructedTracksValidation") << "edm::Handle< RPFittedTrackCollection > is invalid";


        RPFittedTrackCollection::const_iterator it ;


        it = tracks_of_reference->find(referenceRP) ;   // we have only 1 detector for test

        if((it !=  tracks_of_reference->end()))
        {
                const set<unsigned int> &dets = geom->DetsInRP(testedRP);
                for (set<unsigned int>::iterator dit = dets.begin(); dit != dets.end(); ++dit) // Loop over the detectors of the tested RomanPot.
                {
                        // Calculating the intersecting point of the track and the detector plane in global coordinates.
                        unsigned int rawId = TotRPDetId::DecToRawId(*dit);
                        DetGeomDesc *ginfo = geom->GetDetector(rawId);

			TVector3 direction = (it->second).GetDirectionVector() ; // we use the angle information for the reference track

                        if(sqrt(pow(direction.X(),2) + pow(direction.Y(),2)) < tolerance_in_angle) // If the angle is too large, we drop this track
			{
	                        double z = ginfo->translation().z() ;
        	                double x = (it->second).X0() + direction.X() * (z - (it->second).Z0()) ;
                	        double y = (it->second).Y0() + direction.Y() * (z - (it->second).Z0()) ;

                        	// Converting the common point into local coordinates
	                        CLHEP::Hep3Vector hit(x, y, z);
        	                hit = geom->GlobalToLocal(rawId, hit);
                	        double u = hit.x(); double v = hit.y();

                        	// Was this hit within the detector ?
	                        if (RPTopology::IsHit(u, v))  // if it was not a hit, then get the next detector. Now we do not use this.
				{
        	                // Calculating the w coordinate of the hit and fill the histograms
	                	        pair<unsigned int, unsigned int> key (testedRP,rawId) ;

		                        double w = (u+v) / sqrt(2.) ;                   // Calculating the w coordinate...
        				double m = stripZeroPosition - v;
        				signed int strip = (int) floor(m / RPTopology::pitch_);


	        	                if(selectedTest == 1)
        	        	                histograms[key]["Expected"]->Fill(w,1) ;        // ...and filling out the histogram. 
	              	        	else if(selectedTest == 4)
        	              	        	histograms[key]["Expected"]->Fill(strip,1) ;        // ...and filling out the histogram. 

	                	        Handle< edm::DetSetVector<RPRecoHit> > allHits ; // Now we get the hit information for the actual detector.
        	                	e.getByLabel(detSetVectorRPRecoHitLabel, allHits) ;

	        	                // Here we loop over all detected (reconstructed) hits, and try to find at least one among them which is "close" to the actual expected hit.
        	        	        // This part of the algorithm would be faster.
                	        	bool accepted = false ;  // I use this variable to break the nested 2 loop if I accept the actual expected hit
	                        	for (DetSetVector<RPRecoHit>::const_iterator dit2 = allHits->begin(); (dit2 != allHits->end()) && !accepted ; ++dit2)   // Loop over all DetSets.
		                                for(DetSet<RPRecoHit>::const_iterator hit2 = dit2->begin(); (hit2 != dit2->end()) && !accepted ; ++hit2)        // Loop over all hits of one DetSet.
                	                	{
        	        	                        if(hit2->DetId() != rawId) continue ; // If the given hit NOT belongs to the actual detector just continue.
	
        	                	                if((v - 3 * sqrt(2)* resolution_in_w) < (hit2->Position()) && ((hit2->Position()) < (v + 3* sqrt(2)* resolution_in_w))) // If the detected hit close enough to the exptected one, then...
                	                	        {
	                                        	        if(selectedTest == 1)
        	                                        	        histograms[key]["Accepted"]->Fill(w,1) ;        // ...it is accepted.
		              	                                else if(selectedTest == 4)
        		              	                                histograms[key]["Accepted"]->Fill(strip,1) ;
                	                	                accepted = true ;                                       // I accept the actual expected hit just once time.
                        	                	}
                                		}
				}
			}
                }
        }
	
}

void EdgeEfficiency::testTwoRPTestBeam(const edm::Event& e)
{

// for test beam
//	testGivenPairOfRPs(121,120,e) ;
//	testGivenPairOfRPs(120,121,e) ;

//	if(selectedTest==4)
	{
                testGivenPairOfRPs(20,24,e) ;
                testGivenPairOfRPs(21,25,e) ;
                testGivenPairOfRPs(24,20,e) ;
                testGivenPairOfRPs(25,21,e) ;

                testGivenPairOfRPs(22,23,e) ;
                testGivenPairOfRPs(23,22,e) ;

                testGivenPairOfRPs(120,124,e) ;
                testGivenPairOfRPs(121,125,e) ;
                testGivenPairOfRPs(124,120,e) ;
                testGivenPairOfRPs(125,121,e) ;

                testGivenPairOfRPs(122,123,e) ;
                testGivenPairOfRPs(123,122,e) ;

	}
/*	else
	{
	        testGivenPairOfRPs(20,22,e) ;
        	testGivenPairOfRPs(21,22,e) ;
	        testGivenPairOfRPs(24,23,e) ;
        	testGivenPairOfRPs(25,23,e) ;

	        testGivenPairOfRPs(22,20,e) ;
        	testGivenPairOfRPs(23,24,e) ;

	        testGivenPairOfRPs(120,122,e) ;
        	testGivenPairOfRPs(121,122,e) ;
        	testGivenPairOfRPs(124,123,e) ;
	        testGivenPairOfRPs(125,123,e) ;

	        testGivenPairOfRPs(122,120,e) ;
        	testGivenPairOfRPs(123,124,e) ;
	}*/
}

void EdgeEfficiency::testPackages(const edm::Event& e)
{
        edm::Handle< RPFittedTrackCollection > tracks_of_reference ; 
        e.getByLabel(rPFittedTrackCollectionLabel, tracks_of_reference) ;

        if (!tracks_of_reference.isValid())
                throw cms::Exception("RPReconstructedTracksValidation") << "edm::Handle< RPFittedTrackCollection > is invalid";


        RPFittedTrackCollection::const_iterator it ;
        for(it = tracks_of_reference->begin(); it != tracks_of_reference->end()  ; ++it)
	{
                const set<unsigned int> &dets = geom->DetsInRP(it->first);
                for (set<unsigned int>::iterator dit = dets.begin(); dit != dets.end(); ++dit) 
                {
                        // Calculating the intersecting point of the track and the detector plane in global coordinates.
                        unsigned int rawId = TotRPDetId::DecToRawId(*dit);
                        DetGeomDesc *ginfo = geom->GetDetector(rawId);


			TVector3 direction(it->second.GetTx(), it->second.GetTy(), 1) ;

                        double z = ginfo->translation().z() ;
                        double x = (it->second).X0() + direction.X() * (z - (it->second).Z0()) ;
                        double y = (it->second).Y0() + direction.Y() * (z - (it->second).Z0()) ;

//			cout << "I am here" ;

                        // Converting the common point into local coordinates
                        CLHEP::Hep3Vector hit(x, y, z);
                        hit = geom->GlobalToLocal(rawId, hit);
                        double u = hit.x(); double v = hit.y();

			Handle <edm::DetSetVector<RPRecoHit> > allHits ; 
		        e.getByLabel(detSetVectorRPRecoHitLabel, allHits) ;

//                        if (RPTopology::IsHit(u, v))  // if it was not a hit, then get the next detector. Now we do not use this.
                        {
	                        pair<unsigned int, unsigned int> key (it->first,rawId) ;

        	                double w = (u+v) / sqrt(2.) ;                   // Calculating the w coordinate...

//				cout << "I am here2" ;

                        	histograms[key]["Expected"]->Fill(w,1) ;        // ...and filling out the histogram. 

//				cout << "I am here3" ;

        	                bool accepted = false ;  // I use this variable to break the nested 2 loop if I accept the actual expected hit
                	        for (DetSetVector<RPRecoHit>::const_iterator dit2 = allHits->begin(); (dit2 != allHits->end()) && !accepted ; ++dit2)   
                        	        for(DetSet<RPRecoHit>::const_iterator hit2 = dit2->begin(); (hit2 != dit2->end()) && !accepted ; ++hit2)        
                                	{
                                        	if(hit2->DetId() != rawId) continue ; 

	                                        if((v - 2*tolerance_in_v) < (hit2->Position()) && ((hit2->Position()) < (v + 2*tolerance_in_v))) 
        	                                {
                        	                        histograms[key]["Accepted"]->Fill(w,1) ;        // ...it is accepted.
                	                                accepted = true ;                                       // I accept the actual expected hit just once time.
                                	        }

	                                }
				if(accepted==false)
				{	
					// cout << "Event number: " << e.id().event() << " is not accepted by detector: " << *dit <<  endl ; 
				}
			}
		}
        }
}

void EdgeEfficiency::analyze(const edm::Event& e, const edm::EventSetup& eSetup)
{
/*Test !!! REMOVE IT FROM here...

        edm::Handle< RPFittedTrackCollection > tracks_of_reference ; // RPFittedTrackCollection is part of the input. We will use only the point information X0,Y0,Z0.
        e.getByLabel(rPFittedTrackCollectionLabel, tracks_of_reference);

        if (!tracks_of_reference.isValid())
                throw cms::Exception("RPReconstructedTracksValidation") << "edm::Handle< RPFittedTrackCollection > is invalid";


        RPFittedTrackCollection::const_iterator it ;
//	ofstream outfile("/afs/cern.ch/user/f/fnemes/scratch/EdgeEfficieny.log",ios::out | ios::app) ;	

	for(it = tracks_of_reference->begin(); it != tracks_of_reference->end()  ; ++it)
	{
		cout << "RP id: " <<  it->first << endl ;
	}

        Handle< edm::DetSetVector<RPRecoHit> > allHits ; // Now we get the hit information for the actual detector.
        e.getByLabel(detSetVectorRPRecoHitLabel, allHits) ;

        for (DetSetVector<RPRecoHit>::const_iterator dit2 = allHits->begin(); dit2 != allHits->end() ; ++dit2)   // Loop over all DetSets.
        	for(DetSet<RPRecoHit>::const_iterator hit2 = dit2->begin(); hit2 != dit2->end() ; ++hit2)        // Loop over all hits of one DetSet.
                {
			cout << "RecoHit DetId: " << hit2->DetId() << endl ;
			cout << "RecoHit RawId: " << TotRPDetId::RawToDecId(hit2->DetId()) << endl ;
                }

...TO HERE*/

	if((selectedTest == 0) || (selectedTest == 3))
	{
		map<unsigned int, pair<unsigned int, unsigned int> >::iterator tested_RP_it ; // This iterator points to the actual RP which is tested.

		// Loop over all of the detectors which we want to test.
		for(tested_RP_it = map_from_tested_RP_to_its_reference_RPs.begin(); tested_RP_it != map_from_tested_RP_to_its_reference_RPs.end()  ; ++tested_RP_it)
			test(e, tested_RP_it->first) ;
	}
	else if((selectedTest == 1) || (selectedTest == 4))
	{ 
		testTwoRPTestBeam(e) ;
	}
	else if((selectedTest == 2) || (selectedTest == 5))
	{
		testPackages(e) ;
	}
}

DEFINE_FWK_MODULE(EdgeEfficiency);

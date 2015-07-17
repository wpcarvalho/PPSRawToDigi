/****************************************************************************
 *
 * This is a part of TOTEM testbeam/monitoring/reconstruction software.
 * Authors: 
 *	 Giuseppe Latino
 *    
 * $Id: T2DetClustReconst.cc,v 1.2 2008/09/16 15:08:58 lgrzanka Exp $
 * $Revision: 1.2 $
 * $Date: 2008/09/16 15:08:58 $
 *
 ****************************************************************************/
//Created by Giuseppe Latino: reconstruction of clusters (September 06).
//Hubert Niewiadomski: cluster storage added (November 06).
//Giuseppe Latino: geometrical reconstruction of clusters (May 07).


#include "TMath.h"
#include "TotemT1T2Validation/T2GeometryValidation/interface/T2DetClustReconsT.h"


using namespace std;

T2DetClustReconsT::T2DetClustReconsT()
{

} // T2DetClustReconsT Constructor


T2DetClustReconsT::~T2DetClustReconsT()
{
} // T2DetClustReconsT Destructor


void T2DetClustReconsT::SetStripHits(vector<int> strip_ch)
{
  //by HN
  StrHits.clear();
  
  // StripChaNumb.Fill(strip_ch.size());
  for(unsigned int i=0; i<strip_ch.size(); ++i) {
    //   StripChaID.Fill(strip_ch[i]);
    if(strip_ch[i]<=256) {            // originale
      //if(strip_ch[i]<=255) {
      StrHits[strip_ch[i]-1][0] = 1; // originale
      // StrHits[strip_ch[i]][0]=1;
    }else{
      StrHits[strip_ch[i]-257][1] = 1; // originale
      //StrHits[strip_ch[i]-256][1]=1;   
    }
  }
}

void T2DetClustReconsT::SetPadHits(vector<int> pad_ch)
{  
  //by HN
  PadHits.clear();
 
  int ipad, jpad;
  float fpad_ch;
  //  PadChaNumb.Fill(pad_ch.size());
  for(unsigned int i=0; i<pad_ch.size(); ++i) {
    // PadChaID.Fill(pad_ch[i]);
    fpad_ch = pad_ch[i];
    if(pad_ch[i] <= 24) {
      ipad = pad_ch[i];
    }else{
      ipad = int( fmod(fpad_ch,24) );
      if(ipad == 0) ipad = 24;
    }
    jpad = (pad_ch[i]-ipad)/24 + 1;
    PadHits[ipad-1][jpad-1] = 1;
  }
}

void T2DetClustReconsT::BuildClusters() 
{ 
  ClearClusterContainers();
  int iClu = 0, CluSize = 0; 
  int ir = 0;
  bool clus = 0;
  bool goclus[3][3];
  vector<int> StrCluSize, PadsCluSize;
  //
  // 
  // **** Strip Clustering **** 
  //
  // cout<<"Etrato in BuildCluster Strip before = "<<StrHits.size()<<endl;
  for(int j=0; j<2; ++j) {
    for(int i=0; i<256; ++i) {
      if(StrHits[i][j]) {
	T2Cluster strip_cluster;
	strip_cluster.SetDetID(det_id_);
	strip_cluster.SetClusterType(T2Cluster::strip);     
	StrCluStrRawID[iClu][CluSize] = i;
	StrCluStrColID[iClu][CluSize] = j;
	strip_cluster.AddEntry(cluster_entry(i,j));
	CluSize  += 1;
	StrHits[i][j] = 0;
	ir = 1;
	do{
	  clus = 0;
	  if(i+ir < 256 && StrHits[i+ir][j]) {
	    StrCluStrRawID[iClu][CluSize] = i+ir;
	    StrCluStrColID[iClu][CluSize] = j;
	    CluSize += 1;
	    strip_cluster.AddEntry(cluster_entry(i+ir,j));
	    StrHits[i+ir][j] = 0;
	    clus = 1;
	  }
	  ir++;
	} while(clus);
	iClu++;
	StrCluSize.push_back(CluSize);
	CluSize = 0;
	strip_cluster.ComputeClusterParams();
	strip_clusters.push_back(strip_cluster);
      }
    } // end loop on strip raws
  }  // end loop on strip columns 
  //
  //  StripCluNumb.Fill(StrCluSize.size());
  // for(unsigned int i=0; i<StrCluSize.size(); ++i){
  //  StripCluSize.Fill(StrCluSize[i]);
  //}
  //
  //


  // cout<<"FROM T2DETCLUSTRECONST  # of Strip Clusters after  = "<<StrCluSize.size()<<endl;
  /*  cout<<""<<endl;
      for(int i=0; i<StrCluSize.size(); ++i) {
      cout<<"Strip Clus. # "<<i+1<<" => Strip # = "<<StrCluSize[i]<<endl;
      for(int j=0; j<StrCluSize[i]; j++ ) {
      cout<<"Str. # "<<j+1<<" => StripRawID = "<<StrCluStrRawID[i][j]<<
      "\t"<<"StripColID = "<<StrCluStrColID[i][j]<<endl;
      }
      cout<<""<<endl;
      }
  */

  //
  //
  // **** Pad Clustering **** (6X6 window, at least...)
  //
  //cout<<"Pad before = "<<PadHits.size()<<endl;
  iClu = 0;
  CluSize = 0;
  for(Int_t j = 0; j < 65; j++) {
    for(Int_t i = 0; i < 24; i++) {
      if(PadHits[i][j]){
	T2Cluster pad_cluster;
	pad_cluster.SetDetID(det_id_);
	pad_cluster.SetClusterType(T2Cluster::pad);
	PadCluPadRawID[iClu][CluSize] = i;
	PadCluPadColID[iClu][CluSize] = j;
	pad_cluster.AddEntry(cluster_entry(i,j));
	CluSize += 1;
	PadHits[i][j] = 0;
	for(Int_t jr = 0; jr < 3; jr++){
	  for(Int_t jc = -1; jc < 2; jc++){
	    goclus[jr][jc+1] = 0;
	    if(i+jr < 24 && j+jc > -1 && j+jc < 65){
	      if(PadHits[i+jr][j+jc]){
		if(jr < 2) goclus[jr][jc+1] = 1;
		if(jr == 2){
		  if(PadHits[i+1][j+jc]) goclus[jr][jc+1] = 1;
		  if(jc == -1 && PadHits[i+1][i+jc+1] && j+jc+1 < 65) goclus[jr][jc+1] = 1;
		  if(jc == 0 && PadHits[i+1][i+jc+1] && j+jc+1 < 65) goclus[jr][jc+1] = 1;
		  if(jc == 0 && PadHits[i+1][i+jc-1] && j+jc-1 > -1) goclus[jr][jc+1] = 1;
		  if(jc == 1 && PadHits[i+1][i+jc-1] && j+jc-1 > -1) goclus[jr][jc+1] = 1;
		}
	      }
	    }
	  }
	}
	for(Int_t jr = 0; jr < 3; jr++){
	  for(Int_t jc = -1; jc < 2; jc++){
	    if(goclus[jr][jc+1]){
	      for(Int_t l = -1; l < 2; l++){
		if(i+jr+l > -1 && i+jr+l < 24 && j+jc > -1 && j+jc < 65){
		  if(PadHits[i+jr+l][j+jc]){
		    PadCluPadRawID[iClu][CluSize] = i+jr+l;
		    PadCluPadColID[iClu][CluSize] = j+jc;
		    pad_cluster.AddEntry(cluster_entry(i+jr+l,j+jc));
		    CluSize += 1;
		    PadHits[i+jr+l][j+jc] = 0;
		  }
		}
		if(i+jr > -1 && i+jr < 24 && j+jc+l > -1 && j+jc+l < 65){
		  if(PadHits[i+jr][j+jc+l]){
		    PadCluPadRawID[iClu][CluSize] = i+jr;
		    PadCluPadColID[iClu][CluSize] = j+jc+l;
		    pad_cluster.AddEntry(cluster_entry(i+jr,j+jc+l));
		    CluSize += 1;
		    PadHits[i+jr][j+jc+l] = 0;
		  }
		}
		if(i+jr+l > -1 && i+jr+l < 24 && j+jc+l > -1 && j+jc+l < 65){
		  if(PadHits[i+jr+l][j+jc+l]){
		    PadCluPadRawID[iClu][CluSize] = i+jr+l;
		    PadCluPadColID[iClu][CluSize] = j+jc+l;
		    pad_cluster.AddEntry(cluster_entry(i+jr+l,j+jc+l));
		    CluSize += 1;
		    PadHits[i+jr+l][j+jc+l] = 0;
		  }
		}
		if(i+jr-l > -1 && i+jr-l < 24 && j+jc+l > -1 && j+jc+l < 65){
		  if(PadHits[i+jr-l][j+jc+l]){
		    PadCluPadRawID[iClu][CluSize] = i+jr-l;
		    PadCluPadColID[iClu][CluSize] = j+jc+l;
		    pad_cluster.AddEntry(cluster_entry(i+jr-l,j+jc+l));
		    CluSize += 1;
		    PadHits[i+jr-l][j+jc+l] = 0;
		  }
		}
	      }
	    }
	  }
	}
	iClu++;
	PadsCluSize.push_back(CluSize);
	CluSize = 0;
	pad_cluster.ComputeClusterParams();
	pad_clusters.push_back(pad_cluster);
      }  // endif on "seed" channel for pad clustering
    }   // end loop on pad raws
  }    // end loop on pad columns
  //
  // PadCluNumb.Fill(PadsCluSize.size());
  //for(unsigned int i=0; i<PadsCluSize.size(); ++i){
  // PadCluSize.Fill(PadsCluSize[i]);
  // }
  //
  //
  //
  // Evaluate T2 Plane "Hits" (from which tracks will be reconstructed)
  //
  /*
    double StripCluRmin, StripCluRmax, StripCluPhimin, StripCluPhimax;
    double PadCluRmin, PadCluRmax, PadCluPhimin, PadCluPhimax;
    for(unsigned int l=0; l<strip_clusters.size(); ++l){
    T2RecHit Hit;
    Hit.AddEntry(strip_clusters[l]);
    StripCluRmin = strip_clusters[l].GetClusterR() - strip_clusters[l].GetClusterDR();
    StripCluRmax = strip_clusters[l].GetClusterR() + strip_clusters[l].GetClusterDR();
    StripCluPhimin = strip_clusters[l].GetClusterPhi() - strip_clusters[l].GetClusterDPhi();
    StripCluPhimax = strip_clusters[l].GetClusterPhi() + strip_clusters[l].GetClusterDPhi();
    for(unsigned int k=0; k<pad_clusters.size(); ++k){
    PadCluRmin = pad_clusters[k].GetClusterR() - pad_clusters[k].GetClusterDR();
    PadCluRmax = pad_clusters[k].GetClusterR() + pad_clusters[k].GetClusterDR();
    PadCluPhimin = pad_clusters[k].GetClusterPhi() - pad_clusters[k].GetClusterDPhi();
    PadCluPhimax = pad_clusters[k].GetClusterPhi() + pad_clusters[k].GetClusterDPhi();
    //
    if(StripCluRmax <= PadCluRmax || StripCluRmin >= PadCluRmin){
    Hit.AddEntry(pad_clusters[k]);
    }else{
    
    }
    //
    }
    }
  */


  //  cout<<"da T2DetClusReconst: # of Pad Clusters  after = "<<PadsCluSize.size()<<endl;
  /*  cout<<""<<endl;
      for(int i=0; i<PadsCluSize.size(); ++i) {
      cout<<"Pad Clus. # "<<i+1<<" => Pad # = "<<PadsCluSize[i]<<endl;
      for(int j=0; j<PadsCluSize[i]; j++ ) {
      cout<<"Pad # "<<j+1<<" => PadRawID = "<<PadCluPadRawID[i][j]<<
      "\t"<<"PadColID = "<<PadCluPadColID[i][j]<<endl;
      }
      cout<<""<<endl;
      }
  */
  //
  //
}


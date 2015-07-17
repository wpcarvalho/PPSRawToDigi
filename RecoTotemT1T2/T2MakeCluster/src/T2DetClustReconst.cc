/****************************************************************************
 *
 * This is a part of TOTEM testbeam/monitoring/reconstruction software.
 * Authors: 
 *	 Giuseppe Latino
 *       Mirko Berretti (new pad clusterization)
 *    
 * $Id: T2DetClustReconst.cc,v 1.2 2008/09/16 15:08:58 lgrzanka Exp $
 * $Revision: 1.2 $
 * $Date: 2008/09/16 15:08:58 $
 *
 ****************************************************************************/
//Created by Giuseppe Latino: reconstruction of clusters (September 06).
//Hubert Niewiadomski: cluster storage added (November 06).
//Giuseppe Latino: geometrical reconstruction of clusters (May 07).
//Mirko Berretti: Pad clusterization based on nearest neighbour approach (Dec 10) for track reco purpose

#include "TMath.h"
#include "RecoTotemT1T2/T2MakeCluster/interface/T2DetClustReconst.h"


using namespace std;

T2DetClustReconst::T2DetClustReconst()
{
  
} // T2DetClustReconst Constructor




T2DetClustReconst::~T2DetClustReconst()
{
} // T2DetClustReconst Destructor


void T2DetClustReconst::SetStripHits(vector<int> strip_ch)
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

void T2DetClustReconst::SetPadHits(vector<int> pad_ch)
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

void T2DetClustReconst::BuildClusters() 
{ 
  ClearClusterContainers();
  int iClu = 0, CluSize = 0; 
  int ir = 0;
  bool clus = 0;
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
 

  //clusterId, vector::associateduniquepad_id;
  std::map<unsigned int, std::vector<unsigned int> > PadClusterMap; 
  unsigned int clustercounter=0;
  std::map<unsigned int, std::vector<unsigned int> >::iterator PadClusterMapIt;
  std::map<unsigned int, std::vector<unsigned int> >::iterator PadClusterMapIter1;
  std::map<unsigned int, std::vector<unsigned int> >::iterator PadClusterMapIter2;
  
  
  unsigned int paduniqueid=0;
  bool padassociated=false;
  int vecti=0;
  int vectj=0;
  std::vector<unsigned int> newcluster;
  // std::cout<<"CLUSTERING PAD:"<<std::endl;
  for(int i=0;i<24;i++){ //Start iterating rows (to always find the beginning of the next cluster to be reco'd)
    for(int j=0;j<65;j++){     
      if(PadHits[i][j]){
	//std::cout<<"looking in the plane, pad: i-j"<<i<<"-"<<j<<std::endl;
	 paduniqueid=i+65*j;	       
	 padassociated=false;	
	 // std::cout<<"Map size:"<<PadClusterMap.size()<<" NumPadClu="<<clustercounter<<std::endl;
	 for(PadClusterMapIt=PadClusterMap.begin();PadClusterMapIt!=PadClusterMap.end();PadClusterMapIt++){
	   //std::vector<unsigned int> uniqueidvect=PadClusterMapIt->second;
	  for(unsigned int u=0;u<PadClusterMapIt->second.size();u++){
	      vecti=PadClusterMapIt->second.at(u)%65;
	      vectj=PadClusterMapIt->second.at(u)/65;
	      //std::cout<<"Map have already: "<<vecti<<"-"<<vectj<<std::endl;
	      //change this condition to change cluster building criteria
	      if(((vecti==i)&&(abs(vectj-j)==1)) ||((vectj==j)&&(abs(vecti-i)==1)))        
		{
		  PadClusterMapIt->second.push_back(paduniqueid);
		  // std::vector<unsigned int>::iterator it= find(PadClusterMapIt->second.begin(),PadClusterMapIt->second.end(),paduniqueid);
		  //std::cout<<"clusterized"<<std::endl;
		  padassociated=true;
		  break;
		}
	    }
	  if(padassociated==true)
	    break; 
	}
	
	if(padassociated==false)
	  {
	    newcluster.clear();
	    newcluster.push_back(paduniqueid);
	    clustercounter++; 
	    //pair <unsigned int, std::vector<unsigned int> > Pair_to_insert;
	    //PadClusterMap[clustercounter]=Pair_to_insert;
	    PadClusterMap[clustercounter]=newcluster;
	  }	      
      }
    }
  }

  //The clustering is not finished since you can have broken some clusters.example if XXX is read in the loop like:
  //XXX
  //OOX
  //XXX
  //Here going from up to down he will find 2 separated cluster.
  
  unsigned int size1=0;unsigned int size2=0; std::vector<unsigned int> mergedV;
  std::map<unsigned int, std::vector<unsigned int> > PadClusterMapRewrite;
 
  bool foundclustertomerge=true;unsigned int newclusterid=0;
  while(foundclustertomerge==true){
    //each loop will rewrite the full map, at each loop a single merging will be done.
    
    newclusterid++;

    PadClusterMapRewrite=PadClusterMap;

    foundclustertomerge=false;
    for(PadClusterMapIter1=PadClusterMapRewrite.begin();PadClusterMapIter1!=PadClusterMapRewrite.end();PadClusterMapIter1++){
      
      for(PadClusterMapIter2=PadClusterMapIter1;PadClusterMapIter2!=PadClusterMapRewrite.end();PadClusterMapIter2++){
	if(PadClusterMapIter1!=PadClusterMapIter2)
	  {
	    //Lets see if the 2 clusters touch each other
	    size1=PadClusterMapIter1->second.size();size2=PadClusterMapIter2->second.size();
	    for(unsigned int c1=0;c1<size1;c1++){
	      for(unsigned int c2=0;c2<size2;c2++)
		{
		  int i1=PadClusterMapIter1->second.at(c1)%65;
		  int j1=PadClusterMapIter1->second.at(c1)/65;
		  
		  int i2=PadClusterMapIter2->second.at(c2)%65;
		  int j2=PadClusterMapIter2->second.at(c2)/65;
		 
		  if(((i1==i2)&&(abs(j1-j2)==1))||((j1==j2)&&(abs(i1-i2)==1)))
		    {
		      foundclustertomerge=true;
		      break;
		    }		 
		}
	      if(foundclustertomerge)
		break;
	    }

	    if(foundclustertomerge)
	       {
		 mergedV.clear();
		 mergedV=PadClusterMapIter1->second;
		 for(unsigned int c2=0;c2<size2;c2++)
		   mergedV.push_back(PadClusterMapIter2->second.at(c2));
		 break;
	       }
	    
	  }

	if(foundclustertomerge)
	  break;
      }

      if(foundclustertomerge)
	{
	  //Save the merged cluster and remove the wrongs in the original map. Return to the while loop
	  int indextoremove1=PadClusterMapIter1->first;
	  int indextoremove2=PadClusterMapIter2->first;
	  PadClusterMap.erase(indextoremove1);PadClusterMap.erase(indextoremove2);
	  PadClusterMap[indextoremove1]=mergedV;
	  //PadClusterMapRewrite[indextoremove1]=mergedV;
	  break;
	}
      //else
      //{
	  //Save the cluster pointed by PadClusterMapIter1
      //  PadClusterMapRewrite[newclusterid]=mergedV;
      //}
    }
    
  }  
  
  
  //Rewrite the pad cluster collection in a T2 Cluster object.
  //std::map<unsigned int, std::vector<unsigned int> >::iterator PadClusterMapIt;
  unsigned int numpadinclu=0;int rowpad=0;int colpad=0;
  for(PadClusterMapIt=PadClusterMapRewrite.begin();PadClusterMapIt!=PadClusterMapRewrite.end();PadClusterMapIt++){
    T2Cluster pad_cluster;
   
    pad_cluster.SetDetID(det_id_);
    pad_cluster.SetClusterType(T2Cluster::pad);
    numpadinclu=PadClusterMapIt->second.size();
    for(unsigned int pn=0;pn<numpadinclu;pn++)
      {
	int uniquepadid=PadClusterMapIt->second.at(pn);
	rowpad=uniquepadid%65;
	colpad=uniquepadid/65;
	pad_cluster.AddEntry(cluster_entry(rowpad,colpad));
      }
    
    //real padrow col will be 0..23 - 0..64 but I used 65 to make the unique pad-id.
   
   
    pad_cluster.ComputeClusterParamsForTracking(Projectionthreshold,BlobMinSize);
    //pad_cluster.ComputeClusterParams();
    pad_clusters.push_back(pad_cluster); //Add cluster to collection
  }


}
















 /*
//OLD clustering
  std::vector<int> NextRow; //Vectors for storing which pads need to be checked for the cluster under progress
  std::vector<int> NextCol;
  bool notdone=true;
  int ThisRow=0; //Row being checked
  int ThisCol=0; //Column being checked
  for(int i=0;i<24;i++){ //Start iterating rows (to always find the beginning of the next cluster to be reco'd)
    for(int j=0;j<65;j++){ //Start iterating columns
      if(PadHits[i][j]){ //Start building cluster if pad is on
	notdone=true;
        ThisRow=i;
        ThisCol=j;
	PadHits[ThisRow][ThisCol]=0; //Mark pad as off once in a cluster to avoid it being added to another cluster (should not happen, of course...) or multiple times to the same one
	T2Cluster pad_cluster;
        pad_cluster.SetDetID(det_id_);
        pad_cluster.SetClusterType(T2Cluster::pad);
        pad_cluster.AddEntry(cluster_entry(ThisRow,ThisCol)); //Add pad to cluster
	//std::cout << std::endl << "New cluster starts at: " << ThisRow << ", " << ThisCol << std::endl;
	int CluSize=1;
	do{ //Start checking all of the neighbouring pads
	  if(ThisRow!=0){ //These are to stop the algorithm from "jumping over the edge" causing a seg fault
	    if(PadHits[ThisRow-1][ThisCol]){//Up
	      pad_cluster.AddEntry(cluster_entry(ThisRow-1,ThisCol));
	      PadHits[ThisRow-1][ThisCol]=0;
	      NextRow.push_back(ThisRow-1);
	      NextCol.push_back(ThisCol);
	      CluSize++;
	    }
	    if(ThisCol!=0){
	      if(PadHits[ThisRow-1][ThisCol-1]){//Up&Left
		pad_cluster.AddEntry(cluster_entry(ThisRow-1,ThisCol-1));
		PadHits[ThisRow-1][ThisCol-1]=0;
		NextRow.push_back(ThisRow-1);
		NextCol.push_back(ThisCol-1);
		CluSize++;
	      }
	    }
	  }
	  if(ThisCol!=0){
	    if(PadHits[ThisRow][ThisCol-1]){//Left
	      pad_cluster.AddEntry(cluster_entry(ThisRow,ThisCol-1));
	      PadHits[ThisRow][ThisCol-1]=0;
	      NextRow.push_back(ThisRow);
	      NextCol.push_back(ThisCol-1);
	      CluSize++;
	    }
	    if(ThisRow!=23){
	      if(PadHits[ThisRow+1][ThisCol-1]){//Down&Left
		pad_cluster.AddEntry(cluster_entry(ThisRow+1,ThisCol-1));
		PadHits[ThisRow+1][ThisCol-1]=0;
		NextRow.push_back(ThisRow+1);
		NextCol.push_back(ThisCol-1);
		CluSize++;
	      }
	    }
	  }
	  if(ThisRow!=23){
	    if(PadHits[ThisRow+1][ThisCol]){//Down
	      pad_cluster.AddEntry(cluster_entry(ThisRow+1,ThisCol));
	      PadHits[ThisRow+1][ThisCol]=0;
	      NextRow.push_back(ThisRow+1);
	      NextCol.push_back(ThisCol);
	      CluSize++;
	    }
	    if(ThisCol!=64){
	      if(PadHits[ThisRow+1][ThisCol+1]){//Down&Right
		pad_cluster.AddEntry(cluster_entry(ThisRow+1,ThisCol+1));
		PadHits[ThisRow+1][ThisCol+1]=0;
		NextRow.push_back(ThisRow+1);
		NextCol.push_back(ThisCol+1);
		CluSize++;
	      }
	    }
	  }
	  if(ThisCol!=64){
	    if(PadHits[ThisRow][ThisCol+1]){//Right
	      pad_cluster.AddEntry(cluster_entry(ThisRow,ThisCol+1));
	      PadHits[ThisRow][ThisCol+1]=0;
	      NextRow.push_back(ThisRow);
	      NextCol.push_back(ThisCol+1);
	      CluSize++;
	    }
	    if(ThisRow!=0){
	      if(PadHits[ThisRow-1][ThisCol+1]){//Up&Right
		pad_cluster.AddEntry(cluster_entry(ThisRow-1,ThisCol+1));
		PadHits[ThisRow-1][ThisCol+1]=0;
		NextRow.push_back(ThisRow-1);
		NextCol.push_back(ThisCol+1);
		CluSize++;
	      }
	    }
	  }
	  if(NextRow.size()!=0){
	    ThisRow=NextRow[0]; //One pad & neighbours checked now, go on to the next pad in the cluster that has unchecked neighbouts
	    ThisCol=NextCol[0];
	    NextRow.erase(NextRow.begin()); //Remove this from the list of pads to be checked
	    NextCol.erase(NextCol.begin());
	    //std::cout << "Moving on to: " << ThisRow << ", " << ThisCol << std::endl;
	  }
	  else notdone=false; //nothing to check anymore, cluster done!
	}while(notdone);
	//std::cout << CluSize << std::endl;
        pad_cluster.ComputeClusterParams();
        pad_clusters.push_back(pad_cluster); //Add cluster to collection
      }
    }
  }
  */

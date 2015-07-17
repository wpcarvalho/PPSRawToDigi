/****************************************************************************
*
* 
* Authors: Mirko Berretti
* mirko.berretti@gmail.com	
* University of Siena   
* $Id: T2DetHitReconst.cc,v 1.4.2.1 2009/11/07 20:09:23 berretti Exp $ 
* $Date: 2009/11/07 20:09:23 $
*
****************************************************************************/


#include "TMath.h"

#include "RecoTotemT1T2/T2RecHit/interface/T2DetHitReconst.h"




T2DetHitReconst::T2DetHitReconst(uint32_t det_id, unsigned int  Cl1MaxPad_, unsigned int  Cl1MaxStrip_)
{
  det_id_ = det_id;
  //currentDet = new T2DetId(det_id);
  currentDet=T2DetId(det_id);
  //std::cout<<"Constructor rawid "<<det_id<<std::endl;
  Cl1MaxPad= Cl1MaxPad_;
  Cl1MaxStrip= Cl1MaxStrip_;
 
  // this->ClearClusterContainers();
  //this->CalculateZDet();



   //CalculateZDet
  
  double z1=13817;//As written in Z13100_Forward.pdf // was13828.3 //14035.605; //first Gem first drift gas zone (mm)
  double planedist= 91;//91;//86.0;//  è 91 !!
  double btbdist=24.6;   //25.0;
  double ovdist=46;//43.0;
  double zinsidedet= 0;//1.5;  //Put 0!
  int myhalftele= currentDet.halfTelescope();
  int myplaneside= currentDet.planeSide();
  int myplane= currentDet.plane();
  double zdetshift=myplane*planedist+myplaneside*btbdist;
  
  // ht==0 innner; ht==1 outer respect to the ring center.

  if((myhalftele==0)&&(currentDet.arm()==0))
    zdetshift=zdetshift+ovdist;
 
  if((myhalftele==1)&&(currentDet.arm()==1))
    zdetshift=zdetshift+ovdist;


  Zdet=zdetshift+z1+zinsidedet;

  if(currentDet.arm()==1)
    Zdet=Zdet*(-1.0);

  DZdet=1.5;
  


  //ClearClusterContainers
  
  pad_clusters.clear(); 
  strip_clusters.clear(); 



}

 


T2DetHitReconst::~T2DetHitReconst()
{
} 




void T2DetHitReconst::FindHits()                       
{
  
  std::vector<std::vector<int> > IndCombinatedCluster; IndCombinatedCluster.clear();
  std::vector<unsigned int> combinatedStripClusterindex; combinatedStripClusterindex.clear();
  
  std::vector<int> rowClusterindex;
  for (unsigned int i = 0; i < pad_clusters.size(); i++) {
    rowClusterindex.clear();
    rowClusterindex.push_back(i);

    
    for (unsigned int k = 0; k < strip_clusters.size(); k++) {
      
      //std::cout<<"T2RecHit2: Pad RPhi-Strip RPhi:"<<pad_clusters[i].GetClusterR()<<" "<<pad_clusters[i].GetClusterPhi()<<" - "<<strip_clusters[k].GetClusterR()<<" "<<strip_clusters[k].GetClusterPhi()<<std::endl;

      if(((fabs(pad_clusters[i].GetClusterPhi()>179.5))&&(fabs(pad_clusters[i].GetClusterPhi()<180.5)))||((fabs(pad_clusters[i].GetClusterPhi()<0.5))||(fabs(pad_clusters[i].GetClusterPhi()>359.5))))
	{
	
	   if(fabs(pad_clusters[i].GetClusterR() - strip_clusters[k].GetClusterR())< pad_clusters[i].GetClusterDR())
	      {
		rowClusterindex.push_back(k);
		combinatedStripClusterindex.push_back(k);
	      }	 
	}
	else
	  if(fabs(pad_clusters[i].GetClusterPhi() - strip_clusters[k].GetClusterPhi())< strip_clusters[k].GetClusterDPhi())
	    if(fabs(pad_clusters[i].GetClusterR() - strip_clusters[k].GetClusterR())< pad_clusters[i].GetClusterDR())
	      {
		/*
		if(pad_clusters[i].ClustId_unique == 3378){
		  std::cout<<"Cluster 3378 Processed in FindHits i="<<i<<" rowClusterindex.size:"<<rowClusterindex.size()<<" IndCombinatedCluster.size():"<< IndCombinatedCluster.size() <<"Cross check: "<<pad_clusters[rowClusterindex[0]].ClustId_unique<<"Has a trip associated with Id:"<<strip_clusters[k].ClustId_unique<<std::endl;      
		  }
		*/

		rowClusterindex.push_back(k);
		combinatedStripClusterindex.push_back(k);
	      }	     
    }
    
   
    
    IndCombinatedCluster.push_back(rowClusterindex);

    

  }

  // Add the unmatched strip cluster index
  
  bool stripactive;
  std::vector<unsigned int> unmatchedstripindex;
  //unmatchedstripindex.clear();
  for (unsigned int k = 0; k < strip_clusters.size(); k++) {
    stripactive=false;	  
    for (unsigned int l = 0; l < combinatedStripClusterindex.size(); l++) {
      if (k==combinatedStripClusterindex[l])
	stripactive=true;
    }
    if(stripactive==false)
      {
	unmatchedstripindex.push_back(k);
      }  
  }
	
  //	std::cout<<"new hit test"<<std::endl;
	//OnePlaneHits
  
	
	uint32_t det_id;
	T2DetId giveNum;



	// class 0=best; 1=pad 2=strip;
	for (unsigned int k = 0; k < IndCombinatedCluster.size(); k++) {
	  
	  if(IndCombinatedCluster[k].size()==1)// hit made by  1 pad clu
	    {
	      
	      T2Hit onehit;
	      std::vector<int> onerowmatrix = IndCombinatedCluster[k];
	      /*
	      if(pad_clusters[onerowmatrix[0]].ClustId_unique == 3378)
		std::cout<<"\\\\\\ The problematic cluster Pad Is processed inside T2DetHitRec.cc A"<<std::endl;
	      */

	      onehit.HitUniqueId=std::make_pair(pad_clusters[onerowmatrix[0]].ClustId_unique,-1);
	      onehit.SetHitR(pad_clusters[onerowmatrix[0]].GetClusterR());
	      onehit.SetHitPhi(pad_clusters[onerowmatrix[0]].GetClusterPhi());
	      onehit.SetHitZ(Zdet);   
	      onehit.SetHitDR(pad_clusters[onerowmatrix[0]].GetClusterDR());
	      onehit.SetHitDPhi(pad_clusters[onerowmatrix[0]].GetClusterDPhi());
	      onehit.SetHitDZ(DZdet); // to check
	      onehit.SetHitClass(2);
	      onehit.SetHitNumStrip(0);           
	      onehit.SetHitNumPad(pad_clusters[onerowmatrix[0]].GetNoOfEntries());  
	      // std::cout<<"a-"<<pad_clusters[onerowmatrix[0]].GetNoOfEntries()<<"||"<<std::endl;
	      //insert for noise studies
	      onehit.AddCluPadEntries(pad_clusters[onerowmatrix[0]].GetEntries());
	      //end

	      onehit.SetHitX(onehit.GetHitR()*cos(onehit.GetHitPhi()*3.14159265/180.0));
	      onehit.SetHitY(onehit.GetHitR()*sin(onehit.GetHitPhi()*3.14159265/180.0));
	      //The error on X,Y are from propagation formula	  	 
	      onehit.SetHitDX(sqrt(onehit.GetHitDR()*onehit.GetHitDR()*cos(onehit.GetHitPhi()*3.14159265/180.0)*cos(onehit.GetHitPhi()*3.14159265/180.0)+onehit.GetHitDPhi()*3.14159265/180.0*onehit.GetHitDPhi()*3.14159265/180.0*onehit.GetHitR()*onehit.GetHitR()*sin(onehit.GetHitPhi()*3.14159265/180.0)*sin(onehit.GetHitPhi()*3.14159265/180.0)));         	  
	      onehit.SetHitDY(sqrt(onehit.GetHitDR()*onehit.GetHitDR()*sin(onehit.GetHitPhi()*3.14159265/180.0)*sin(onehit.GetHitPhi()*3.14159265/180.0)+onehit.GetHitDPhi()*3.14159265/180.0*onehit.GetHitDPhi()*3.14159265/180.0*onehit.GetHitR()*onehit.GetHitR()*cos(onehit.GetHitPhi()*3.14159265/180.0)*cos(onehit.GetHitPhi()*3.14159265/180.0)));

	      det_id= pad_clusters[onerowmatrix[0]].GetDetID();
	      onehit.SetHitArm((giveNum.arm(det_id))); 
	      onehit.SetHitHalftele((giveNum.halfTelescope(det_id))); 
	      onehit.SetHitPlane((giveNum.plane(det_id))); 
	      onehit.SetHitPlaneSide((giveNum.planeSide(det_id)));
	      onehit.SetHitDetRawId(det_id);
	  

	      theHitV.push_back(onehit);
	      
	    }
	      
	  
	  if(IndCombinatedCluster[k].size()==2)
	    {
	      // hit made by  1 pad e 1 strip
	      T2Hit onehit;
	      //  if(pad_clusters[onerowmatrix[0]].ClustId_unique == 3378)
	      //std::cout<<"\\\\\\ The problematic cluster Pad Is processed inside T2DetHitRec.cc B"<<std::endl;

	      std::vector<int> onerowmatrix = IndCombinatedCluster[k];
	      onehit.HitUniqueId=std::make_pair(pad_clusters[onerowmatrix[0]].ClustId_unique,strip_clusters[onerowmatrix[1]].ClustId_unique);
	      onehit.SetHitR(strip_clusters[onerowmatrix[1]].GetClusterR());
	      onehit.SetHitPhi(pad_clusters[onerowmatrix[0]].GetClusterPhi());
	      onehit.SetHitZ(Zdet); 
	      onehit.SetHitDR(strip_clusters[onerowmatrix[1]].GetClusterDR());
	      onehit.SetHitDPhi(pad_clusters[onerowmatrix[0]].GetClusterDPhi());
	      onehit.SetHitDZ(DZdet); // to check
	      onehit.SetHitClass(1);
    	      onehit.SetHitNumStrip(strip_clusters[onerowmatrix[1]].GetNoOfEntries());           
	      onehit.SetHitNumPad(pad_clusters[onerowmatrix[0]].GetNoOfEntries());  
	      //   std::cout<<"b-"<<pad_clusters[onerowmatrix[0]].GetNoOfEntries()<<"||"<<std::endl;
	      //insert for noise studies
	      onehit.AddCluPadEntries(pad_clusters[onerowmatrix[0]].GetEntries());
	      //std::cout<<"T2Hit: Added a vector of pad of "<<onehit.ClusterPad_entries.size()<<std::endl;
	      // for(unsigned int m=0; m<onehit.ClusterPad_entries.size();m++)
	      //	{
		  //std::cout<<"RadCord: "<<onehit.ClusterPad_entries.at(m).rad_coord<<std::endl;
		  //std::cout<<"ColCord: "<<onehit.ClusterPad_entries.at(m).ang_coord<<std::endl;
	      //	}
	      onehit.AddCluStripEntries(strip_clusters[onerowmatrix[1]].GetEntries());
	      //end

	      //Removed on 28/12/2010.
	      //if((onehit.GetHitNumStrip()>Cl1MaxStrip)||(onehit.GetHitNumPad()>Cl1MaxPad))
	      //{
	      //  onehit.SetHitClass(0);
	      //}

	 
	      onehit.SetHitX(onehit.GetHitR()*cos(onehit.GetHitPhi()*3.14159265/180.0));
	      onehit.SetHitY(onehit.GetHitR()*sin(onehit.GetHitPhi()*3.14159265/180.0));
	      //The error on X,Y are from propagation formula	  	 
	      onehit.SetHitDX(sqrt(onehit.GetHitDR()*onehit.GetHitDR()*cos(onehit.GetHitPhi()*3.14159265/180.0)*cos(onehit.GetHitPhi()*3.14159265/180.0)+onehit.GetHitDPhi()*3.14159265/180.0*onehit.GetHitDPhi()*3.14159265/180.0*onehit.GetHitR()*onehit.GetHitR()*sin(onehit.GetHitPhi()*3.14159265/180.0)*sin(onehit.GetHitPhi()*3.14159265/180.0)));         	  
	      onehit.SetHitDY(sqrt(onehit.GetHitDR()*onehit.GetHitDR()*sin(onehit.GetHitPhi()*3.14159265/180.0)*sin(onehit.GetHitPhi()*3.14159265/180.0)+onehit.GetHitDPhi()*3.14159265/180.0*onehit.GetHitDPhi()*3.14159265/180.0*onehit.GetHitR()*onehit.GetHitR()*cos(onehit.GetHitPhi()*3.14159265/180.0)*cos(onehit.GetHitPhi()*3.14159265/180.0)));


	      det_id= pad_clusters[onerowmatrix[0]].GetDetID();
	      onehit.SetHitArm((giveNum.arm(det_id))); 
	      onehit.SetHitHalftele((giveNum.halfTelescope(det_id))); 
	      onehit.SetHitPlane((giveNum.plane(det_id))); 
	      onehit.SetHitPlaneSide((giveNum.planeSide(det_id)));
	      onehit.SetHitDetRawId(det_id);
	  


	      theHitV.push_back(onehit);	      
	      
	    }
	  
	  if(IndCombinatedCluster[k].size()>2)
	    {
	      // Reconstruction of more than one hit built with 1 pad e 1 strip cluster starting from 1 Pad & 1 or more strip
	      std::vector<int> onerowmatrix = IndCombinatedCluster[k];
	      
	      
	      //This part increase a lot the Hit pad clu size: try ll<2 and you will get pad cluster results
	      
	      for(unsigned int ll=1;ll<IndCombinatedCluster[k].size();ll++)
		{
		  T2Hit onehit;
		  onehit.HitUniqueId=std::make_pair(pad_clusters[onerowmatrix[0]].ClustId_unique,strip_clusters[onerowmatrix[ll]].ClustId_unique);
		  onehit.SetHitR(strip_clusters[onerowmatrix[ll]].GetClusterR());
		  onehit.SetHitPhi(pad_clusters[onerowmatrix[0]].GetClusterPhi());
		  		  

		  onehit.SetHitZ(Zdet);
		  onehit.SetHitDR(strip_clusters[onerowmatrix[ll]].GetClusterDR());
		  onehit.SetHitDPhi(pad_clusters[onerowmatrix[0]].GetClusterDPhi());
		  onehit.SetHitDZ(DZdet);// da mettere
		  onehit.SetHitClass(1);
		  onehit.SetHitNumStrip(strip_clusters[onerowmatrix[ll]].GetNoOfEntries());   //bug fix: before l<->1 09-06-11       
		  onehit.SetHitNumPad(pad_clusters[onerowmatrix[0]].GetNoOfEntries());

		  //insert for noise studies
		  //onehit.ClusterStrip_entries
		  // std::cout<<"c-"<<pad_clusters[onerowmatrix[0]].GetNoOfEntries()<<"||"<<std::endl;
		  onehit.AddCluStripEntries(strip_clusters[onerowmatrix[ll]].GetEntries());
		  onehit.AddCluPadEntries(pad_clusters[onerowmatrix[0]].GetEntries());
		  //end

		  //Removed on 28/12/2010.
		  //if((onehit.GetHitNumStrip()>Cl1MaxStrip)||(onehit.GetHitNumPad()>Cl1MaxPad))
		  // {
		  //   onehit.SetHitClass(0);
		  // }				  

		  onehit.SetHitX(onehit.GetHitR()*cos(onehit.GetHitPhi()*3.14159265/180.0));
		  onehit.SetHitY(onehit.GetHitR()*sin(onehit.GetHitPhi()*3.14159265/180.0));
		  //The error on X,Y are from propagation formula	  	 
		  onehit.SetHitDX(sqrt(onehit.GetHitDR()*onehit.GetHitDR()*cos(onehit.GetHitPhi()*3.14159265/180.0)*cos(onehit.GetHitPhi()*3.14159265/180.0)+onehit.GetHitDPhi()*3.14159265/180.0*onehit.GetHitDPhi()*3.14159265/180.0*onehit.GetHitR()*onehit.GetHitR()*sin(onehit.GetHitPhi()*3.14159265/180.0)*sin(onehit.GetHitPhi()*3.14159265/180.0)));         	  
		  onehit.SetHitDY(sqrt(onehit.GetHitDR()*onehit.GetHitDR()*sin(onehit.GetHitPhi()*3.14159265/180.0)*sin(onehit.GetHitPhi()*3.14159265/180.0)+onehit.GetHitDPhi()*3.14159265/180.0*onehit.GetHitDPhi()*3.14159265/180.0*onehit.GetHitR()*onehit.GetHitR()*cos(onehit.GetHitPhi()*3.14159265/180.0)*cos(onehit.GetHitPhi()*3.14159265/180.0)));
		  
		  det_id= pad_clusters[onerowmatrix[0]].GetDetID();
		  onehit.SetHitArm((giveNum.arm(det_id))); 
		  onehit.SetHitHalftele((giveNum.halfTelescope(det_id))); 
		  onehit.SetHitPlane((giveNum.plane(det_id))); 
		  onehit.SetHitPlaneSide((giveNum.planeSide(det_id)));
		  onehit.SetHitDetRawId(det_id);
	
		  theHitV.push_back(onehit);
		 
		}
	    }
	  	 	      	    	 	    
	}

	
	// Hit Reconstruction from 1 strip cluster
	  for (unsigned int k = 0; k < unmatchedstripindex.size(); k++) {  
	    T2Hit onehit;	 
	    
	    onehit.HitUniqueId=std::make_pair(-1,strip_clusters[unmatchedstripindex[k]].ClustId_unique);
	    

       	      onehit.SetHitR(strip_clusters[unmatchedstripindex[k]].GetClusterR());
	      onehit.SetHitPhi(strip_clusters[unmatchedstripindex[k]].GetClusterPhi());
	      onehit.SetHitZ(Zdet);   // da mettere
	      onehit.SetHitDR(strip_clusters[unmatchedstripindex[k]].GetClusterDR());
	      onehit.SetHitDPhi(strip_clusters[unmatchedstripindex[k]].GetClusterDPhi());
	      onehit.SetHitDZ(DZdet); //da mettere 
	      onehit.SetHitClass(2);
	      
	      onehit.SetHitNumStrip(strip_clusters[unmatchedstripindex[k]].GetNoOfEntries());           
	      onehit.SetHitNumPad(0);
	      //std::cout<<"c-"<<strip_clusters[unmatchedstripindex[k]]..GetNoOfEntries()<<"||"<<std::endl;
	      //insert for noise studies
	      onehit.AddCluStripEntries(strip_clusters[unmatchedstripindex[k]].GetEntries());
	      //end

	      onehit.SetHitX(onehit.GetHitR()*cos(onehit.GetHitPhi()*3.14159265/180.0));
	      onehit.SetHitY(onehit.GetHitR()*sin(onehit.GetHitPhi()*3.14159265/180.0));
	      //The error on X,Y are from propagation formula	  	 
	      onehit.SetHitDX(sqrt(onehit.GetHitDR()*onehit.GetHitDR()*cos(onehit.GetHitPhi()*3.14159265/180.0)*cos(onehit.GetHitPhi()*3.14159265/180.0)+onehit.GetHitDPhi()*3.14159265/180.0*onehit.GetHitDPhi()*3.14159265/180.0*onehit.GetHitR()*onehit.GetHitR()*sin(onehit.GetHitPhi()*3.14159265/180.0)*sin(onehit.GetHitPhi()*3.14159265/180.0)));         	  
	      onehit.SetHitDY(sqrt(onehit.GetHitDR()*onehit.GetHitDR()*sin(onehit.GetHitPhi()*3.14159265/180.0)*sin(onehit.GetHitPhi()*3.14159265/180.0)+onehit.GetHitDPhi()*3.14159265/180.0*onehit.GetHitDPhi()*3.14159265/180.0*onehit.GetHitR()*onehit.GetHitR()*cos(onehit.GetHitPhi()*3.14159265/180.0)*cos(onehit.GetHitPhi()*3.14159265/180.0)));
	      

	      det_id= strip_clusters[unmatchedstripindex[k]].GetDetID();
	      onehit.SetHitArm((giveNum.arm(det_id))); 
	      onehit.SetHitHalftele((giveNum.halfTelescope(det_id))); 
	      onehit.SetHitPlane((giveNum.plane(det_id))); 
	      onehit.SetHitPlaneSide((giveNum.planeSide(det_id)));
	      onehit.SetHitDetRawId(det_id);

	      
	      theHitV.push_back(onehit);	 
	      
	      //AlleventPlaneHits.push_back(onehit);
	    }


      //ahit.ComputeHit();

     
}



	    


void T2DetHitReconst::CalculateZDet()
{

  //CalculateZDet
  /*
double z1=13828.3;          //14035.605; //first Gem first drift gas zone (mm)
double planedist= 86.0;
double btbdist=24.6;   //25.0;
double ovdist=43.0;
 double zinsidedet= 1.5;  
int myhalftele= currentDet.halfTelescope();
int myplaneside= currentDet.planeSide();
int myplane= currentDet.plane();
double zdetshift=myplane*planedist+myplaneside*btbdist;

// ht==0 innner; ht==1 outer respect to the ring center.

 if((myhalftele==0)&&(currentDet.arm()==0))
   zdetshift=zdetshift+ovdist;
 
 if((myhalftele==1)&&(currentDet.arm()==1))
   zdetshift=zdetshift+ovdist;


 Zdet=zdetshift+z1+zinsidedet;

 if(currentDet.arm()==1)
   Zdet=Zdet*(-1.0);

DZdet=1.5;
  */
}

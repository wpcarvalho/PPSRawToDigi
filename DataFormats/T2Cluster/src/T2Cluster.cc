/****************************************************************************
 *
 * This is a part of TOTEM testbeam/monitoring/reconstruction software.
 * Authors: 
 *	Hubert Niewiadomski
 *    
 * $Id: T2Cluster.cc,v 1.2.2.1 2009/11/07 20:04:35 berretti Exp $
 * $Revision: 1.2.2.1 $
 * $Date: 2009/11/07 20:04:35 $
 *
 ****************************************************************************/
//Hubert Niewiadomski: the GEM cluster element (Nov. 06)
//Giuseppe Latino: geometrical reconstruction of GEM clusters added (May 07).
//Mirko Berretti ComputeClusterParamsForTracking (Dec 010).

#include "DataFormats/T2Cluster/interface/T2Cluster.h"
#include "DataFormats/T2Cluster/interface/T2ROGeometry.h"
#include <TMath.h>


//Projection Threshold Default 0.1
void T2Cluster::ComputeClusterParamsForTracking(double Projectionthreshold,int BlobMinSize)
{
  //Projectionthreshold means that row or column are not counted in the average if they contribute less than
  // Projectionthreshold to the integral (cluster size). Utilized only for pad clusters.

  double StrPitch = 0.400;
  entry_numb_ = cluster_entries.size();


  if(entry_numb_==0)
    {
      radial_spread_ = angular_spread_ = 0;
      radial_centre_ = angular_centre_ = 0.0;
      std::cout<<"Warning, problem in ComputeClusterParamsForTracking"<<std::endl;
      return;
    }



  bool isblob=false;
  if(type_ == pad)
    {
      if(entry_numb_>=BlobMinSize)
	{
	  isblob=true;
	}
    }
  

  if(isblob)//Make particular row-col weight;
    {

      unsigned int row_beg=0;  unsigned int row_end=0; unsigned int col_beg=0;  unsigned int col_end=0;
      row_end = cluster_entries[0].rad_coord;
      col_end = cluster_entries[0].ang_coord;
      row_beg = row_end;
      col_beg = col_end;

      radial_centre_ = cluster_entries[0].rad_coord;
      angular_centre_ = cluster_entries[0].ang_coord;

     double theintegral=(double)entry_numb_;
     std::map<int,vector<int> > clurowsMap;  
     std::map<int,vector<int> > clucolsMap;
     clurowsMap.clear(); 
     clucolsMap.clear();
      //std::map<int,vector<int> >::iterator itmapr; std::map<int,vector<int> >::iterator itmapc;
      for(unsigned int i=0; i<cluster_entries.size(); i++)
	{	  	  
	    clurowsMap[cluster_entries[i].rad_coord].push_back(cluster_entries[i].ang_coord);
	    clucolsMap[cluster_entries[i].ang_coord].push_back(cluster_entries[i].rad_coord);
	}

    //Now the Blob Map is loaded. I remove from the map the row-col which is contributing less 
    //than Projectionthreshold to the integral (clusize). If we will end with nothing the threshold will be reduced by factor 2.
      double truncatingThrPadRow=Projectionthreshold;
      double truncatingThrPadCol=Projectionthreshold;
      
      bool blobcleaned=false;
      std::map<int,vector<int> >  FinalCleanclurowsMap;  
      std::map<int,vector<int> >  FinalCleanclucolsMap; 
      
      
      while(blobcleaned==false)
	{
	 
	  std::map<int,vector<int> >::iterator itmapr; 	  
	  FinalCleanclurowsMap.clear();  
	  std::vector<int> rowtoberemoved;  
	 
	  for(itmapr=clurowsMap.begin();itmapr!=clurowsMap.end();itmapr++)
	    {
	      int entries=(*itmapr).second.size();
	      if(entries<truncatingThrPadRow*theintegral)
		{
		  rowtoberemoved.push_back((*itmapr).first);		 
		}
	      else
		FinalCleanclurowsMap[(*itmapr).first]=(*itmapr).second;
	    }


	  std::map<int,vector<int> >::iterator itmapc;
	  std::vector<int> coltoberemoved;  
	  FinalCleanclucolsMap.clear();
	  
	  for(itmapc=clucolsMap.begin();itmapc!=clucolsMap.end();itmapc++)
	    {
	      int entries=(*itmapc).second.size();
	      if(entries<truncatingThrPadCol*theintegral)
		{
		  coltoberemoved.push_back((*itmapc).first);
		}
	      else
		FinalCleanclucolsMap[(*itmapc).first]=(*itmapc).second;
	    }

	  
	  if((FinalCleanclurowsMap.size()>0)&&(FinalCleanclucolsMap.size()>0))
	    {
	      blobcleaned=true;	      
	    }
	  else{
	    if(FinalCleanclurowsMap.size()==0)
	      truncatingThrPadRow=(truncatingThrPadRow/2.0);
	    if(FinalCleanclucolsMap.size()==0)
	      truncatingThrPadCol=(truncatingThrPadCol/2.0);
	  }
	  
	} //while blobcleaned==false

      //Now the blob is cleaned and there is a map <rowNumb,padcols> <colNumb,padrows> from which start the geo convertion
      //(std::map<int,int> FinalCleanclurowsMap;  std::map<int,int> FinalCleanclucolsMap)

      

      T2ROGeometry convert(det_id_);
       std::map<int,vector<int> >::iterator itmapr; 
       std::vector<double> AllRadialDR; std::vector<double> AllRadialcenters; std::vector<double> AllRowscenters;
       for(itmapr=FinalCleanclurowsMap.begin();itmapr!=FinalCleanclurowsMap.end();itmapr++){
	 int thisrow=itmapr->first;
	 std::vector<int> associatedcol =itmapr->second;
	 
	 
	 for(unsigned int w=0;w<associatedcol.size();w++)
	   {
	     int thiscol=associatedcol.at(w);
	     double myR_centre_  = (convert.GetPadRMin(thisrow,thiscol)+ convert.GetPadRMax(thisrow,thiscol))/2.;
	     double myDR_centre_  = fabs(convert.GetPadRMax(thisrow,thiscol) - convert.GetPadRMin(thisrow,thiscol))/2.;
	     AllRadialcenters.push_back(myR_centre_);
	     AllRadialDR.push_back(myDR_centre_);
	     AllRowscenters.push_back(thisrow);
	   }
	 
       }

       std::map<int,vector<int> >::iterator itmapc;
       std::vector<double> AllAzimDphi;
       std::vector<double> AllAzimuthcenters;
       std::vector<double> AllColscenters;
       for(itmapc=FinalCleanclucolsMap.begin();itmapc!=FinalCleanclucolsMap.end();itmapc++){
	 int thiscol=itmapc->first;
	 std::vector<int> associatedrows=itmapc->second;

	 for(unsigned int w=0;w<associatedrows.size();w++)
	   {
	     int thisrow=associatedrows.at(w);
	     double myPhi_centre_=0.;
	     double myDPhi_centre_=0.;

	     if(convert.GetPadPhiMin(thisrow,thiscol) < convert.GetPadPhiMax(thisrow,thiscol))
	       {
		 myPhi_centre_  = (convert.GetPadPhiMin(thisrow,thiscol) + convert.GetPadPhiMax(thisrow,thiscol))/2.;
		 myDPhi_centre_  = (convert.GetPadPhiMax(thisrow,thiscol) - convert.GetPadPhiMin(thisrow,thiscol))/2.;              
	       }
	     else{	   
	       myPhi_centre_  = (convert.GetPadPhiMin(thisrow,thiscol) + convert.GetPadPhiMax(thisrow,thiscol)+360.)/2.;
	       if(myPhi_centre_ > 360.)
		 myPhi_centre_ = myPhi_centre_ -360.;
	       
	       myDPhi_centre_  = (convert.GetPadPhiMax(thisrow,thiscol)+360. - convert.GetPadPhiMin(thisrow,thiscol))/2.;
	     }

	     AllAzimDphi.push_back(myDPhi_centre_); 
	     AllAzimuthcenters.push_back(myPhi_centre_); 
	     AllColscenters.push_back(thiscol);
	   }
	 
       }

       
       BlobR_DR(AllRadialcenters,AllRadialDR,R_centre_,DR_centre_);
       BlobPhi_Dphi(AllAzimuthcenters,AllAzimDphi,Phi_centre_,DPhi_centre_);
       Blobradial_centre_spread(AllRowscenters,radial_centre_,radial_spread_);
       Blobangular_centre_spread(AllColscenters,angular_centre_,angular_spread_);
       
       //Begin of Cluster-Id-Assignment calculation
      	unsigned int symbplane=0; T2DetId converter; 
	long int uniquePadIDinEvt=0;
       converter=T2DetId(det_id_);
       symbplane= converter.plane()*2 + converter.planeSide() + converter.arm()*20 + converter.halfTelescope()*10; 
    	   
       //uniquePadIDinEvt=symbplane*1560+col_beg*24+row_beg;
       uniquePadIDinEvt=symbplane*1560+col_beg*24+row_beg;

       if((uniquePadIDinEvt<0)||(uniquePadIDinEvt>150000))
	 std::cout<<"Bug at the RECO-clustering level, fix the code"<<std::endl;
       ClustId_unique=uniquePadIDinEvt;	 
       //End of Cluster-Id-Assignment calculation

       if((DR_centre_<0.1)||(DPhi_centre_<1.)){
	 std::cout<<"T2Cluster.cc Blob Err: Dr:"<<DR_centre_<<"DPhi:"<<DPhi_centre_<<"Sizes:"<<AllRadialcenters.size()<<" "<<AllAzimuthcenters.size()<<std::endl;
	 for(unsigned int t=0;t<AllRadialcenters.size();t++)
	   {
	     std::cout<<AllRadialcenters.at(t)<<std::endl;
	   }
       }
    }
  else      //Do as usual, as ComputeClusterParams
    {

       unsigned int row_beg, row_end, col_beg, col_end;
       row_beg = row_end = cluster_entries[0].rad_coord;
       col_beg = col_end = cluster_entries[0].ang_coord;
       radial_centre_ = cluster_entries[0].rad_coord;
       angular_centre_ = cluster_entries[0].ang_coord;
       //  
       for(unsigned int i=1; i<cluster_entries.size(); ++i)
	 {
	   if(row_beg>cluster_entries[i].rad_coord)
	     row_beg = cluster_entries[i].rad_coord;
	   if(row_end<cluster_entries[i].rad_coord)
	     row_end = cluster_entries[i].rad_coord;
	   if(col_beg>cluster_entries[i].ang_coord)
	     col_beg = cluster_entries[i].ang_coord;
	   if(col_end<cluster_entries[i].ang_coord)
	     col_end = cluster_entries[i].ang_coord;
      
	   radial_centre_ += cluster_entries[i].rad_coord;
	   angular_centre_ += cluster_entries[i].ang_coord;
	 }
       //  
       radial_centre_ = radial_centre_/(float) cluster_entries.size();
       angular_centre_ = angular_centre_/(float) cluster_entries.size();
       radial_spread_ = row_end - row_beg + 1;
       angular_spread_ = col_end - col_beg + 1;


        //Begin of Cluster-Id-Assignment calculation
       unsigned int symbplane=0; T2DetId converter; 
       long int uniquePadIDinEvt=0;
       converter=T2DetId(det_id_);
       symbplane= converter.plane()*2 + converter.planeSide() + converter.arm()*20 + converter.halfTelescope()*10;
       uniquePadIDinEvt=symbplane*1560+col_beg*24+row_beg;
       if((uniquePadIDinEvt<0)||(uniquePadIDinEvt>150000))
	 std::cout<<"Bug at the RECO-clustering level, fix the code"<<std::endl;
       ClustId_unique=uniquePadIDinEvt;
	 
       //End of Cluster-Id-Assignment calculation



       // 
       // GL: evaluate geometrical position (R-Phi) of strip/pad clusters
       //
       T2ROGeometry convert(det_id_);
       if(type_ == strip){
	 if(col_beg == col_end){
	   R_centre_  = (convert.GetStripRMin(row_beg,col_beg) + convert.GetStripRMax(row_end,col_beg))/2.;
	   DR_centre_  = StrPitch/sqrt(12.); // standard definition of error for digital readout  
	   //cout<<"FROM T2Cluster: R Centre strip ="<< R_centre_ <<endl;

	   if(convert.GetStripPhiMin(row_beg,col_beg) < convert.GetStripPhiMax(row_end,col_beg)){
	     Phi_centre_  = (convert.GetStripPhiMin(row_beg,col_beg) + convert.GetStripPhiMax(row_end,col_beg))/2.;
	     DPhi_centre_  = (convert.GetStripPhiMax(row_end,col_beg) - convert.GetStripPhiMin(row_beg,col_beg))/2.;
	   }else{
	     Phi_centre_  = (convert.GetStripPhiMin(row_beg,col_beg) + convert.GetStripPhiMax(row_end,col_beg)+360.)/2.;
	     if(Phi_centre_ > 360.) Phi_centre_ = Phi_centre_ -360.;
	     DPhi_centre_  = (convert.GetStripPhiMax(row_end,col_beg)+360. - convert.GetStripPhiMin(row_beg,col_beg))/2.;
	   }
	 }else{
	   // cout<<"T2Cluster: Strip Cluster on 2 Columns !"<<endl;
	 }
       }
  //
       if(type_ == pad){
	 R_centre_  = (convert.GetPadRMin(row_beg,col_beg) + convert.GetPadRMax(row_end,col_beg))/2.;
	 DR_centre_  = (convert.GetPadRMax(row_end,col_beg) - convert.GetPadRMin(row_beg,col_beg))/2.;
	 
	 if(convert.GetPadPhiMin(row_beg,col_beg) < convert.GetPadPhiMax(row_end,col_beg))
	   {
	     //std::cout<<"---Min<Max: Pad Phimin-Phimax:"<<convert.GetPadPhiMin(row_beg,col_beg)<<"-"<<convert.GetPadPhiMax(row_end,col_beg)<<std::endl;
	     
	     double padphimin=convert.GetPadPhiMin(row_beg,col_beg); 
	     double padphimax=convert.GetPadPhiMax(row_end,col_beg);
	     
	     double padphiminbeg=convert.GetPadPhiMin(row_beg,col_beg); 	    	     
	     double padphimaxbeg=convert.GetPadPhiMax(row_beg,col_beg);
	     double padphiminend=convert.GetPadPhiMin(row_beg,col_end); 
	     double padphimaxend=convert.GetPadPhiMax(row_beg,col_end);
	     
	     if((fabs(padphiminbeg-padphiminend)>20.)||(fabs(padphimaxbeg-padphimaxend)>20.)||(fabs(padphiminbeg-padphimaxbeg)>20.)||(fabs(padphimaxbeg-padphiminend)>20.)){
	       if(padphiminbeg<10.)
		 padphiminbeg+=360.;
	       if(padphimaxbeg<10.)
		 padphimaxbeg+=360.;
	       if(padphiminend<10.)
		 padphiminend+=360.;
	       if(padphimaxend<10.)
		 padphimaxend+=360.;
	     }

	     padphimin=padphiminbeg;
	     if(padphimaxbeg<padphimin) padphimin=padphimaxbeg;
	     if(padphiminend<padphimin) padphimin=padphiminend;
	     if(padphimaxend<padphimin) padphimin=padphimaxend;

	     padphimax=padphiminbeg;
	     if(padphimaxbeg>padphimax) padphimax=padphimaxbeg;
	     if(padphiminend>padphimax) padphimax=padphiminend;
	     if(padphimaxend>padphimax) padphimax=padphimaxend;


	     Phi_centre_  = (padphimin+ padphimax)/2.;
	     if(Phi_centre_ > 360.) 
	       Phi_centre_ = Phi_centre_ -360.;
	     DPhi_centre_  = (padphimax- padphimin)/2.;



	     /*
	     Phi_centre_  = (convert.GetPadPhiMin(row_beg,col_beg) + convert.GetPadPhiMax(row_beg,col_end))/2.;
	     DPhi_centre_  = (convert.GetPadPhiMax(row_beg,col_beg) - convert.GetPadPhiMin(row_beg,col_end))/2.;         
	     */
	     //if(DPhi_centre_<1.0)
	     //	 cout<<"T2Cluster: minch rb,re,cb,ce =" <<row_beg  <<"-"<< row_end <<"-" <<  col_beg <<"-"<< col_end << endl;
	     //cout<<"GetPadPhiMax(row_beg,col_beg) - GetPadPhiMin(row_beg,col_end)" <<convert.GetPadPhiMax(row_beg,col_beg)<<" - "<<convert.GetPadPhiMin(row_beg,col_end) << endl        
	     //	cout<<"1- Phi_centre_  = "<<Phi_centre_<<" = (convert.GetPadPhiMin(row_beg,col_beg) + convert.GetPadPhiMax(row_beg,col_end))/2. = ("<<convert.GetPadPhiMin(row_beg,col_beg)<<" + "<<convert.GetPadPhiMax(row_beg,col_end)<<")/2."<<endl;
	   }
	 else{
	   //  std::cout<<"---Min>=Max: Pad Phimin-Phimax:"<<convert.GetPadPhiMin(row_beg,col_beg)<<"-"<<convert.GetPadPhiMax(row_end,col_beg)<<std::endl;
	     
	   
	   // cout<<"2- Phi_centre_  = "<<Phi_centre_<<" = (convert.GetPadPhiMin(row_beg,col_beg) + convert.GetPadPhiMax(row_end,col_beg) +360.)/2. = ("<<convert.GetPadPhiMin(row_beg,col_beg)<<" + "<<convert.GetPadPhiMax(row_end,col_beg)<<"+360)/2."<<endl;


	   //Try to adjust the problem around X axis


	   double padphimin=convert.GetPadPhiMin(row_beg,col_beg); 
	     double padphimax=convert.GetPadPhiMax(row_end,col_beg);
	     
	     double padphiminbeg=convert.GetPadPhiMin(row_beg,col_beg); 	    	     
	     double padphimaxbeg=convert.GetPadPhiMax(row_beg,col_beg);
	     double padphiminend=convert.GetPadPhiMin(row_beg,col_end); 
	     double padphimaxend=convert.GetPadPhiMax(row_beg,col_end);
	     
	     if((fabs(padphiminbeg-padphiminend)>20.)||(fabs(padphimaxbeg-padphimaxend)>20.)||(fabs(padphiminbeg-padphimaxbeg)>20.)||(fabs(padphimaxbeg-padphiminend)>20.)){
	       if(padphiminbeg<10.)
		 padphiminbeg+=360.;
	       if(padphimaxbeg<10.)
		 padphimaxbeg+=360.;
	       if(padphiminend<10.)
		 padphiminend+=360.;
	       if(padphimaxend<10.)
		 padphimaxend+=360.;
	     }

	     padphimin=padphiminbeg;
	     if(padphimaxbeg<padphimin) padphimin=padphimaxbeg;
	     if(padphiminend<padphimin) padphimin=padphiminend;
	     if(padphimaxend<padphimin) padphimin=padphimaxend;

	     padphimax=padphiminbeg;
	     if(padphimaxbeg>padphimax) padphimax=padphimaxbeg;
	     if(padphiminend>padphimax) padphimax=padphiminend;
	     if(padphimaxend>padphimax) padphimax=padphimaxend;


	     Phi_centre_  = (padphimin+ padphimax)/2.;
	     if(Phi_centre_ > 360.) 
	       Phi_centre_ = Phi_centre_ -360.;
	     DPhi_centre_  = (padphimax- padphimin)/2.;

	     /*     
	   double padphimin=convert.GetPadPhiMin(row_beg,col_beg); 
	   double padphimax=convert.GetPadPhiMax(row_end,col_beg);

	   if(fabs(padphimin-padphimax)>20.)
	     {
	       if((padphimin<10)&&(padphimax>350))
		 padphimin=360.0+padphimin;
	       
	       if((padphimax<10)&&(padphimin>350))
		 padphimax=360.0+padphimax;
	       
	       if(padphimin>padphimax){
		 double swp=padphimax;
		 padphimax=padphimin;
		 padphimin=swp;
	       }
	     
	     }
	   Phi_centre_  = (padphimin+ padphimax)/2.;
	   if(Phi_centre_ > 360.) 
	      Phi_centre_ = Phi_centre_ -360.;
	   DPhi_centre_  = (padphimax- padphimin)/2.;
	     */

	   /*
	   Phi_centre_  = (convert.GetPadPhiMin(row_beg,col_beg) + convert.GetPadPhiMax(row_end,col_beg)+360.)/2.;
	   if(Phi_centre_ > 360.) Phi_centre_ = Phi_centre_ -360.;
	   DPhi_centre_  = (convert.GetPadPhiMax(row_end,col_beg)+360. - convert.GetPadPhiMin(row_beg,col_beg))/2.;
	   */ 

	 }
       }
    
    }
  
}



//BlobR_DR(&AllRadialcenters,&AllRadialDR,&R_centre_,&DR_centre_);
//BlobPhi_Dphi(&AllAzimuthcenters,&AllAzimDphi,&Phi_centre_,&DPhi_centre_);
//Blobradial_centre_spread(&AllRowscenters,&radial_centre_,&radial_spread_);
//Blobangular_centre_spread(&AllColscenters,&angular_centre_,&angular_spread_);

void T2Cluster::BlobR_DR(std::vector<double> &AllRadialcenters,std::vector<double> &AllRadialDR,float &R_centre_,float &DR_centre_)
{
  float toret=0.;
  for(unsigned u=0;u<AllRadialcenters.size();u++)
    {
      toret+=AllRadialcenters.at(u);
    }
  toret=toret/(double(AllRadialcenters.size()));
  R_centre_=toret;
  float toret2=0.;
  for(unsigned u=0;u<AllRadialcenters.size();u++)
    {
      toret2=toret2+(R_centre_-AllRadialcenters.at(u))*(R_centre_-AllRadialcenters.at(u));
    }
  toret2=sqrt((toret2)/(double(AllRadialcenters.size())));

  if(toret2<0.115) //could be for a "single strips" pad triplets.
    toret2=0.115;

  DR_centre_=toret2;
}

void T2Cluster::BlobPhi_Dphi(std::vector<double> &AllAzimuthcenters,std::vector<double> &AllAzimDphi,float &Phi_centre_,float &DPhi_centre_)
{
  float toret=0.; int numphiup=0;int numphidown=0;
  float toretA=0.;float toretB=0.;

  for(unsigned u=0;u<AllAzimuthcenters.size();u++)
    {
      toretB+=(AllAzimuthcenters.at(u));
      if(AllAzimuthcenters.at(u)<110.){
	numphiup++;	
	toretA+=(AllAzimuthcenters.at(u)+360.);
      }
      if(AllAzimuthcenters.at(u)>260.){
	toretA+=(AllAzimuthcenters.at(u));
	numphidown++;
      }
      
    }


  if((numphiup>0)&&(numphidown>0))
    {
      toret=toretA;     
    }
  else
    toret=toretB;

 
  toret=toret/(double(AllAzimuthcenters.size()));
  
  float toret2=0.;
  for(unsigned u=0;u<AllAzimuthcenters.size();u++)
    {
      double touse=AllAzimuthcenters.at(u);
      
      if((numphiup>0)&&(numphidown>0))
	if(touse<110.)
	  touse=touse+360.;

      toret2=toret2+(toret-touse)*(toret-touse);
    }
  
  if(toret>360.)
    toret=toret-360.;
  
  Phi_centre_=toret;
  
  if(fabs(Phi_centre_)<0.00001)
    std::cout<<"Warning in BlobPhi_Dphi T2Cluster.cc. Phi="<<Phi_centre_<<" Azim center loaded: "<<AllAzimuthcenters.size()<<std::endl;

  toret2=sqrt((toret2)/(double(AllAzimuthcenters.size())));
  
  if(toret2<1.75) //could be for a "single strips" pad triplets.
    toret2=1.75;
  
  DPhi_centre_=toret2;
}

void T2Cluster::Blobradial_centre_spread(std::vector<double> &AllRowscenters,float &radial_centre_,unsigned short &radial_spread_)
{
  
  float toret=0.;
  for(unsigned u=0;u<AllRowscenters.size();u++)
    {
      toret+=AllRowscenters.at(u);
    }

  radial_centre_=float(toret/(float(AllRowscenters.size())));
  radial_spread_=((unsigned short)(AllRowscenters.size()));
}

void T2Cluster::Blobangular_centre_spread(std::vector<double> &AllColscenters,float &angular_centre_,unsigned short &angular_spread_)
{
 float toret=0.;
 for(unsigned u=0;u<AllColscenters.size();u++)
   {
     toret+=AllColscenters.at(u);
   }

 angular_centre_=float(toret/(float(AllColscenters.size())));
 angular_spread_=((unsigned short)(AllColscenters.size()));
}








void T2Cluster::ComputeClusterParams()
{
  double StrPitch = 0.400;
  entry_numb_ = cluster_entries.size();
  if(entry_numb_==0)
    {
      radial_spread_ = angular_spread_ = 0;
      radial_centre_ = angular_centre_ = 0.0;
      return;
    }

 


  unsigned int row_beg, row_end, col_beg, col_end;
  row_beg = row_end = cluster_entries[0].rad_coord;
  col_beg = col_end = cluster_entries[0].ang_coord;
  radial_centre_ = cluster_entries[0].rad_coord;
  angular_centre_ = cluster_entries[0].ang_coord;
  //  
  for(unsigned int i=1; i<cluster_entries.size(); ++i)
    {
      if(row_beg>cluster_entries[i].rad_coord)
	row_beg = cluster_entries[i].rad_coord;
      if(row_end<cluster_entries[i].rad_coord)
	row_end = cluster_entries[i].rad_coord;
      if(col_beg>cluster_entries[i].ang_coord)
	col_beg = cluster_entries[i].ang_coord;
      if(col_end<cluster_entries[i].ang_coord)
	col_end = cluster_entries[i].ang_coord;
      
      radial_centre_ += cluster_entries[i].rad_coord;
      angular_centre_ += cluster_entries[i].ang_coord;
    }
  //  
  radial_centre_ = radial_centre_/(float) cluster_entries.size();
  angular_centre_ = angular_centre_/(float) cluster_entries.size();
  radial_spread_ = row_end - row_beg + 1;
  angular_spread_ = col_end - col_beg + 1;


  //Begin of Cluster-Id-Assignment calculation This method should also contain the strip case.
  unsigned int symbplane=0; T2DetId converter; 
  long int uniquePadIDinEvt=0;
  converter=T2DetId(det_id_);
  symbplane= converter.plane()*2 + converter.planeSide() + converter.arm()*20 + converter.halfTelescope()*10;
  if(type_ == strip){
    
     uniquePadIDinEvt=100000+symbplane*512+col_beg*128+row_beg;
    if((uniquePadIDinEvt<0)||(uniquePadIDinEvt>150000))
      std::cout<<"Bug at the RECO-clustering level, fix the code"<<std::endl;
    ClustId_unique=uniquePadIDinEvt;	
    
  }
  else
    {
      
       uniquePadIDinEvt=symbplane*1560+col_beg*24+row_beg;
      if((uniquePadIDinEvt<0)||(uniquePadIDinEvt>150000))
	std::cout<<"Bug at the RECO-clustering level, fix the code"<<std::endl;
      ClustId_unique=uniquePadIDinEvt;
    }
  //End of Cluster-Id-Assignment calculation







  // 
  // GL: evaluate geometrical position (R-Phi) of strip/pad clusters
  //
  T2ROGeometry convert(det_id_);
  if(type_ == strip){
    if(col_beg == col_end){
      R_centre_  = (convert.GetStripRMin(row_beg,col_beg) + convert.GetStripRMax(row_end,col_beg))/2.;
      DR_centre_  = StrPitch/sqrt(12.); // standard definition of error for digital readout  
      //cout<<"FROM T2Cluster: R Centre strip ="<< R_centre_ <<endl;

      if(convert.GetStripPhiMin(row_beg,col_beg) < convert.GetStripPhiMax(row_end,col_beg)){
	Phi_centre_  = (convert.GetStripPhiMin(row_beg,col_beg) + convert.GetStripPhiMax(row_end,col_beg))/2.;
	DPhi_centre_  = (convert.GetStripPhiMax(row_end,col_beg) - convert.GetStripPhiMin(row_beg,col_beg))/2.;
      }else{
	Phi_centre_  = (convert.GetStripPhiMin(row_beg,col_beg) + convert.GetStripPhiMax(row_end,col_beg)+360.)/2.;
	if(Phi_centre_ > 360.) Phi_centre_ = Phi_centre_ -360.;
	DPhi_centre_  = (convert.GetStripPhiMax(row_end,col_beg)+360. - convert.GetStripPhiMin(row_beg,col_beg))/2.;
      }
    }else{
      // cout<<"T2Cluster: Strip Cluster on 2 Columns !"<<endl;
    }
  }
  //
  if(type_ == pad){
    R_centre_  = (convert.GetPadRMin(row_beg,col_beg) + convert.GetPadRMax(row_end,col_beg))/2.;
    DR_centre_  = (convert.GetPadRMax(row_end,col_beg) - convert.GetPadRMin(row_beg,col_beg))/2.;

    if(convert.GetPadPhiMin(row_beg,col_beg) < convert.GetPadPhiMax(row_end,col_beg))
      {
	Phi_centre_  = (convert.GetPadPhiMin(row_beg,col_beg) + convert.GetPadPhiMax(row_beg,col_end))/2.;
	DPhi_centre_  = (convert.GetPadPhiMax(row_beg,col_beg) - convert.GetPadPhiMin(row_beg,col_end))/2.;         
	//if(DPhi_centre_<1.0)
	//	 cout<<"T2Cluster: minch rb,re,cb,ce =" <<row_beg  <<"-"<< row_end <<"-" <<  col_beg <<"-"<< col_end << endl;
	//cout<<"GetPadPhiMax(row_beg,col_beg) - GetPadPhiMin(row_beg,col_end)" <<convert.GetPadPhiMax(row_beg,col_beg)<<" - "<<convert.GetPadPhiMin(row_beg,col_end) << endl        
	//	cout<<"1- Phi_centre_  = "<<Phi_centre_<<" = (convert.GetPadPhiMin(row_beg,col_beg) + convert.GetPadPhiMax(row_beg,col_end))/2. = ("<<convert.GetPadPhiMin(row_beg,col_beg)<<" + "<<convert.GetPadPhiMax(row_beg,col_end)<<")/2."<<endl;
      }
    else{
      // cout<<"2- Phi_centre_  = "<<Phi_centre_<<" = (convert.GetPadPhiMin(row_beg,col_beg) + convert.GetPadPhiMax(row_end,col_beg) +360.)/2. = ("<<convert.GetPadPhiMin(row_beg,col_beg)<<" + "<<convert.GetPadPhiMax(row_end,col_beg)<<"+360)/2."<<endl;

      Phi_centre_  = (convert.GetPadPhiMin(row_beg,col_beg) + convert.GetPadPhiMax(row_end,col_beg)+360.)/2.;
      if(Phi_centre_ > 360.) Phi_centre_ = Phi_centre_ -360.;
      DPhi_centre_  = (convert.GetPadPhiMax(row_end,col_beg)+360. - convert.GetPadPhiMin(row_beg,col_beg))/2.;
    }
  }
  //
}





ostream & operator<<(ostream & out, const T2Cluster& clust)
{
  out<<"==Cluster info=="<<endl  
     <<"Entries="<<clust.GetNoOfEntries()<<endl
     <<"Rad spread="<<clust.GetRadialSpread()<<endl
     <<"Ang. spread="<<clust.GetAngularSpread()<<endl
     <<"Rad. centre="<<clust.GetRadialCentrePos()<<endl
     <<"Ang. centre="<<clust.GetAngularCentrePos()<<endl;
  out<<"Cluster points:"<<endl;
  
  for(unsigned int i=0; i<clust.GetEntries().size(); ++i)
    {
      out<<"{"<<clust.GetEntries()[i].rad_coord<<","<<clust.GetEntries()[i].ang_coord<<"}, ";
    }
  out<<endl<<endl;
  return out;
}


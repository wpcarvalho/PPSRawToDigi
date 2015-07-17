/** 
 * Author: Mirko Berretti
 * Email: mirko.berretti@gmail.com
 *
 */

#include "TotemAnalysis/T2Cuts/interface/T2SelectionCutUtils.h"
#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"
#include "TMath.h"

T2SelectionCutUtils::T2SelectionCutUtils()
{}

T2SelectionCutUtils::~T2SelectionCutUtils()
{}


void T2SelectionCutUtils::SetCuts(double T2_TrkEtamin_,double T2_TrkEtaMAX_,int T2_trkMultimin_,int T2_trkMultiMAX_,double T2_DZMultiplicator_,double T2_PhiChiProbCut_,double T2_RChiProbCut_,std::vector<int> T2_QuarterUsed_,bool XYFitUsed_){
  T2_TrkEtamin = T2_TrkEtamin_;
  T2_TrkEtaMAX = T2_TrkEtaMAX_;
  T2_trkMultimin = T2_trkMultimin_;
  T2_trkMultiMAX = T2_trkMultiMAX_ ; 
  T2_DZMultiplicator = T2_DZMultiplicator_ ; 
  T2_PhiChiProbCut =  T2_PhiChiProbCut_;
  T2_RChiProbCut = T2_RChiProbCut_;
  T2_QuarterUsed = T2_QuarterUsed_;
  XYFitUsed= XYFitUsed_;
  IgnoredSmallAnglePrimarySlopeCut=0.0;

  std::cout<<">> T2SelectionCutUtils::SetCuts"<<std::endl;
  std::cout<<"    |Etamin:"<< T2_TrkEtamin <<"|EtaMAX:"<<T2_TrkEtaMAX <<"|trkMultimin:"<<T2_trkMultimin<<"|trkMultiMAX:"<<T2_trkMultiMAX<<"|T2_DZMult:"<<T2_DZMultiplicator<<"|T2_PhiChiProbCut:"<<T2_PhiChiProbCut<<"|T2_RChiProbCut:"<<T2_RChiProbCut<<"|IgnoredSmallAnglePrimarySlopeCut:"<<IgnoredSmallAnglePrimarySlopeCut<<"|XYFitUsed:"<<XYFitUsed;

 dxMisalError=0.;
 dyMisalError=0.;


}

void T2SelectionCutUtils::SetCuts(double T2_TrkEtamin_,double T2_TrkEtaMAX_,int T2_trkMultimin_,int T2_trkMultiMAX_,double T2_DZMultiplicator_,double T2_PhiChiProbCut_,double T2_RChiProbCut_,std::vector<int> T2_QuarterUsed_,double IgnoredSmallAnglePrimarySlopeCut_,bool XYFitUsed_) {
  T2_TrkEtamin = T2_TrkEtamin_;
  T2_TrkEtaMAX = T2_TrkEtaMAX_;
  T2_trkMultimin = T2_trkMultimin_;
  T2_trkMultiMAX = T2_trkMultiMAX_ ; 
  T2_DZMultiplicator = T2_DZMultiplicator_ ; 
  T2_PhiChiProbCut =  T2_PhiChiProbCut_;
  T2_RChiProbCut = T2_RChiProbCut_;
  IgnoredSmallAnglePrimarySlopeCut  = IgnoredSmallAnglePrimarySlopeCut_;
  T2_QuarterUsed = T2_QuarterUsed_;
  XYFitUsed= XYFitUsed_;

  std::cout<<"T2SelectionCutUtils Init:.."<<std::endl;
  std::cout<<"|Etamin:"<< T2_TrkEtamin <<"|EtaMAX:"<<T2_TrkEtaMAX <<"|trkMultimin:"<<T2_trkMultimin<<"|trkMultiMAX:"<<T2_trkMultiMAX<<"|T2_DZMult:"<<T2_DZMultiplicator<<"|T2_PhiChiProbCut:"<<T2_PhiChiProbCut<<"|T2_RChiProbCut:"<<T2_RChiProbCut<<"|IgnoredSmallAnglePrimarySlopeCut:"<<IgnoredSmallAnglePrimarySlopeCut<<"|XYFitUsed:"<<XYFitUsed<<" Quarter used: ";

 for(unsigned int iii=0;iii<T2_QuarterUsed_.size();iii++)
   std::cout<<T2_QuarterUsed_.at(iii)<<" ";

 dxMisalError=0.;
 dyMisalError=0.;
  std::cout<<".. done:"<<std::endl;
}


double T2SelectionCutUtils::FirstPlaneForThishit(T2Hit hit)
{
  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo;
  planeinfo=conv.GetT2Info(hit.GetHitDetRawId());
  unsigned int quarterid=(planeinfo.symb/10); 
  planeinfo=conv.GetT2Info(quarterid*10);  
  double firstplanez=planeinfo.Zdet;

  return firstplanez;
}

/*

  double GetInternalMisErrorDX()
    {
      return dxMisalError;
    }

  double GetInternalMisErrorDY()
    {
      return dyMisalError;
    }

  double SetInternalMisErrorDX(double dxmisalError_)
    {
      dxMisalError=dxmisalError_;
    }

  double SetInternalMisErrorDY(double dxmisalError_)
    {
      dyMisalError=dymisalError_;
    }
*/

std::vector<T2Hit> T2SelectionCutUtils::GetHitsInQuarterFrame(std::vector<T2Hit> hitvec)
{
  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo;
  planeinfo=conv.GetT2Info(hitvec.at(0).GetHitDetRawId());
  unsigned int quarterid=(planeinfo.symb/10); 
  planeinfo=conv.GetT2Info(quarterid*10);  
  double firstplanez=planeinfo.Zdet;
  unsigned int  sizeHitv=hitvec.size();
  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      hitvec.at(jj).SetHitZ(hitvec.at(jj).GetHitZ()-firstplanez);
    }

  return hitvec;
}

std::vector<double> T2SelectionCutUtils::MyLinearfitXYSimple(std::vector<T2Hit> hitvec2)
{  

double Sr=0.;
double Srz=0.;
double Szz_r=0.;
double Sz_r=0.; 
double S0_r=0.; 

double Sphi=0.;
double Sphiz=0.;
double Szz_phi=0.;
double Sz_phi=0.; 
double S0_phi=0.; 

 std::vector<T2Hit> hitvec;

 for(unsigned int jj =0; jj<hitvec2.size(); jj++)
   {
     if(hitvec2.at(jj).GetHitClass()==1)
       {
	 hitvec.push_back(hitvec2.at(jj));
       }
   }

unsigned int  sizeHitv=hitvec.size();

 if(sizeHitv<2)
   {
     std::cout<<" T2SelectionCutUtils::MyLinearfitXYSimple problem: Track with less than 2 Cl1 hits!! Continue the fitting.."<<sizeHitv <<" begining:"<<hitvec2.size()<<std::endl;
     hitvec.clear();
     hitvec=hitvec2;
     sizeHitv=hitvec.size();
   }

 std::vector<double> r;   
std::vector<double> z;
std::vector<double> phi;
std::vector<double> er;
 std::vector<double> ez;
std::vector<double> ephi;
 double a_rz, b_rz, a_phiz, b_phiz;

  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      //   std::cout<<hitvec[jj].GetHitX()<<" "<<hitvec[jj].GetHitPhi()<<" "<<hitvec[jj].GetHitZ()<<std::endl;
      //r->x   phi->y
      r.push_back(hitvec[jj].GetHitX());
      
      phi.push_back(hitvec[jj].GetHitY());
      z.push_back(hitvec[jj].GetHitZ());    
 
      double phirad=hitvec[jj].GetHitPhi()*3.14159/180.0;
      
      //er.push_back(0.1);
      //dalla Formula di JointProb viene sigma2x=r*sigma2R
      //   if(UseJointProb==1)
      //sigmax=((0.12+0.05)*hitvec[jj].GetHitR())*((0.12+0.05)/**hitvec[jj].GetHitR()*/);
      //if(UseJointProb==1)
      //sigmay=((0.12+0.05)*hitvec[jj].GetHitR())*((0.12+0.05)/**hitvec[jj].GetHitR()*/);
      double sigmax;
	double sigmay;
      if(hitvec[jj].GetHitClass()==1)
	{
	  sigmax=cos(phirad)*cos(phirad)*(0.12+0.05)*(0.12+0.05)+sin(phirad)*sin(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
         
	  sigmay=sin(phirad)*sin(phirad)*(0.12+0.05)*(0.12+0.05)+cos(phirad)*cos(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
	}
      else
	{

	  double dr=hitvec[jj].GetHitDR();
	   if(hitvec[jj].GetHitClass()!=9)
	     {
	       sigmax=cos(phirad)*cos(phirad)*dr*dr+sin(phirad)*sin(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();	  
	       sigmay=sin(phirad)*sin(phirad)*dr*dr+cos(phirad)*cos(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();	
	     }
	   else
	     {
	       sigmax=hitvec[jj].GetHitDX()*hitvec[jj].GetHitDX();
	       sigmay=hitvec[jj].GetHitDY()*hitvec[jj].GetHitDY();
	     }
	}
      er.push_back(sqrt(sigmax));  
      ephi.push_back(sqrt(sigmay));
      
     
      //ez.push_back(hitvec[jj].GetHitDZ());
            
      Srz += r[jj]*z[jj]/er[jj]/er[jj];
      Szz_r += z[jj]*z[jj]/er[jj]/er[jj];
      Sz_r += z[jj]/er[jj]/er[jj];
      Sr += r[jj]/er[jj]/er[jj];
      S0_r += 1.0/er[jj]/er[jj];



      Sphiz += phi[jj]*z[jj]/ephi[jj]/ephi[jj];
      Szz_phi += z[jj]*z[jj]/ephi[jj]/ephi[jj];
      Sz_phi += z[jj]/ephi[jj]/ephi[jj];
      Sphi += phi[jj]/ephi[jj]/ephi[jj];
      S0_phi += 1.0/ephi[jj]/ephi[jj];




  }


a_rz = (Srz*S0_r - Sz_r*Sr) / (Szz_r*S0_r - Sz_r*Sz_r);   // angular coefficient
b_rz = (Sr*Szz_r - Sz_r*Srz) / (Szz_r*S0_r - Sz_r*Sz_r);  // intercept   X=(a_rz)Z + b_rz



a_phiz = (Sphiz*S0_phi - Sz_phi*Sphi) / (Szz_phi*S0_phi - Sz_phi*Sz_phi);   // angular coefficient
b_phiz = (Sphi*Szz_phi - Sz_phi*Sphiz) / (Szz_phi*S0_phi - Sz_phi*Sz_phi);  // intercept   


 double e_a_rz = sqrt( S0_r / (S0_r*Szz_r - Sz_r*Sz_r) );
 double e_b_rz = sqrt( Szz_r / (S0_r*Szz_r - Sz_r*Sz_r) );

 double e_a_phiz = sqrt( S0_phi / (S0_phi*Szz_phi - Sz_phi*Sz_phi) );
 double e_b_phiz = sqrt( Szz_phi / (S0_phi*Szz_phi - Sz_phi*Sz_phi) );

 std::vector<double> vect;
 for (unsigned int i=0;i<10;i++)
   vect.push_back(0.0);

 vect[0]=a_rz; //ax
 vect[1]=b_rz;  //bx
 vect[2]=a_phiz; //ay
 vect[3]=b_phiz; //by

 vect[4]=e_a_rz; //Eax
 vect[5]=e_b_rz;  //Ebx
 vect[6]=e_a_phiz; //Eay
 vect[7]=e_b_phiz;  //Eby
 
 double correl=(-1.0)*Sz_r*(1.0/(S0_r*Szz_r-(Sz_r*Sz_r)));
 vect[8]=correl;
 correl=(-1.0)*Sz_phi*(1.0/(S0_phi*Szz_phi-(Sz_phi*Sz_phi)));
 vect[9]=correl;

 //std::cout<<vect[0]<<" "<<vect[1]<<" "<<vect[2]<<" "<<vect[3]<<" "<<std::endl;
 return vect;


}
  



//should return ax,bx,ay,by,...
std::vector<double> T2SelectionCutUtils::MyLinearfitCorr(std::vector<T2Hit> hitvec2,TMatrixD &par_covariance_matrix,double &chi2_)
{


 std::vector<T2Hit> hitvec;
   for(unsigned int jj =0; jj<hitvec2.size(); jj++)
   {
     if(hitvec2.at(jj).GetHitClass()==1)
       {
	 hitvec.push_back(hitvec2.at(jj));
       }
   }
  //std::cout<<"Inside MyLinearfitCorr .. init"<<std::endl;

  double sigmaR=(0.12+0.05);
  double sigmaPhi=0.015;
  unsigned int  sizeHitv=hitvec.size();

  if(sizeHitv<2)
   {
     std::cout<<" T2SelectionCutUtils::MyLinearfitCorr problem: Track with less than 2 Cl1 hits!! Continue the fitting.."<<std::endl;
     hitvec.clear();
     hitvec=hitvec2;
     sizeHitv=hitvec.size();
   }



  TMatrixD ParCov_matr(4,4); //matrice di covarianza dei parametri trovati;
  unsigned int sizeArighe=sizeHitv*2;
  TMatrixD A(sizeArighe,4);
  TMatrixD At(4,sizeArighe);
  //A ?? la matrice per cui Mis=A(Param)
  TMatrixD Vy(sizeArighe,sizeArighe); //matrice di covarianza delle misure (una per ogni xy, quindi ?? diag a blocchi);
  TMatrixD Vym1(sizeArighe,sizeArighe); 

  TMatrixD Ls(4,4);//((A^T)(Vy^-1)A)
  // TMatrixD Cs(4,4);//(A^T)(Vy^-1)
  TMatrixD Cs(4,sizeArighe);//(A^T)(Vy^-1)

  TVectorD FittedParam(4);
 
  TVectorD Yvect(sizeArighe);

  std::vector<TVectorD> allXY; //vettore degli (xi,yi);
  
  for(unsigned int k =0; k<sizeHitv; k++)
    {
      TVectorD insert(2);
      insert[0]=hitvec[k].GetHitX();
      insert[1]=hitvec[k].GetHitY();
      allXY.push_back(insert);
    }

  Ls.Zero();
  TMatrixD Ls00(2,2);//un quarto della matrice LS.
  TMatrixD Ls01(2,2);
  TMatrixD Ls10(2,2);
  TMatrixD Ls11(2,2);

  TMatrixD Mia(2,2);
  TMatrixD Mib(2,2);

  TMatrixD MiaT(2,2);
  TMatrixD MibT(2,2);

  Ls00.Zero();
  Ls01.Zero();
  Ls10.Zero();
  Ls11.Zero();

  Vy.Zero();
  Vym1.Zero();

  //std::cout<<"MyLinearfitCorr Start computation  .. "<<std::endl;

for(unsigned int k =0; k<sizeHitv; k++)
 {
   TMatrixD OneVy(2,2); 
   OneVy.Zero();
   double phirad=hitvec[k].GetHitPhi()*3.14159/180.0;
   double r=hitvec[k].GetHitR();

   if(hitvec[k].GetHitClass()!=1)
     {
       sigmaR=hitvec[k].GetHitDR();
     }
   
   OneVy(0,0)=sigmaR*sigmaR*cos(phirad)*cos(phirad)+r*r*sin(phirad)*sin(phirad)*sigmaPhi*sigmaPhi;
   OneVy(0,1)=cos(phirad)*sin(phirad)*(sigmaR*sigmaR-r*r*sigmaPhi*sigmaPhi);  
   OneVy(1,0)=OneVy(0,1);
   OneVy(1,1)=sigmaR*sigmaR*sin(phirad)*sin(phirad)+ r*r*cos(phirad)*cos(phirad)*sigmaPhi*sigmaPhi;
   
   //Invert Vy matrix
   TMatrixD OneVym1(2,2); 
   OneVym1=OneVy;
   Double_t deti;	
   OneVym1.Invert(&deti);
   
   //OneVy.Print();
   if(fabs(deti)<0.001)     
     std::cout<<"Possible Vy Zero Determinant!!!"<<std::endl;
    
   Vym1(2*k,2*k)=OneVym1(0,0);
   Vym1(2*k,2*k+1)= OneVym1(0,1); 
   Vym1(2*k+1,2*k)=OneVym1(1,0);
   Vym1(2*k+1,2*k+1)=OneVym1(1,1);
   

   Vy(2*k,2*k)=OneVy(0,0);
   Vy(2*k,2*k+1)= OneVy(0,1); 
   Vy(2*k+1,2*k)=OneVy(1,0);
   Vy(2*k+1,2*k+1)=OneVy(1,1);

   



   Yvect(2*k)=hitvec[k].GetHitX();
   Yvect(2*k+1)=hitvec[k].GetHitY();  

   A(2*k,0)=hitvec[k].GetHitZ();
   A(2*k,1)=1.0;
   A(2*k,2)=0.0;
   A(2*k,3)=0.0;
   A(2*k+1,0)=0.0;
   A(2*k+1,1)=0.0;
   A(2*k+1,2)=hitvec[k].GetHitZ();
   A(2*k+1,3)=1.0;

     

   Mia.Zero();
   Mib.Zero();
   
   
   Mia(0,0)=hitvec[k].GetHitZ();
   Mia(0,1)=0.;
   Mia(1,0)=1.;
   Mia(1,1)=0.;
   Mib(0,0)=0.;
   Mib(0,1)=hitvec[k].GetHitZ();
   Mib(1,0)=0.;
   Mib(1,1)=1.;
   
   MiaT=Mia;
   MibT=Mib;
   MiaT.Transpose(MiaT);
   MibT.Transpose(MibT);  
      

   Ls00=Ls00+Mia*OneVym1*MiaT;
   
   Ls01=Ls01+Mia*OneVym1*MibT;
   
   Ls10=Ls10+Mib*OneVym1*MiaT;
   
   Ls11=Ls11+Mib*OneVym1*MibT;  


 }

//std::cout<<"MyLinearfitCorr Parameter calculation  .. "<<std::endl;

//Nota che qui Vy ?? in realt?? Vy^-1 

 Ls(0,0)=Ls00(0,0);  Ls(0,1)=Ls00(0,1);   Ls(0,2)=Ls01(0,0);   Ls(0,3)=Ls01(0,1);
 Ls(1,0)=Ls00(1,0);  Ls(1,1)=Ls00(1,1);   Ls(1,2)=Ls01(1,0);   Ls(1,3)=Ls01(1,1); 

 Ls(2,0)=Ls10(0,0);  Ls(2,1)=Ls10(0,1);   Ls(2,2)=Ls11(0,0);   Ls(2,3)=Ls11(0,1);
 Ls(3,0)=Ls10(1,0);  Ls(3,1)=Ls10(1,1);   Ls(3,2)=Ls11(1,0);   Ls(3,3)=Ls11(1,1);

 
 Double_t determ;	
 Ls.Invert(&determ);
 if(fabs(determ)<0.001)
     std::cout<<"Possible Zero LS  Determinant!!!"<<std::endl;

 At.Zero();
 At.Transpose(A);
 Cs.Zero();
 Cs=At*Vym1;


 
 FittedParam= Ls*Cs*Yvect;
 
 //std::cout<<"MyLinearfitCorr FittedParam: "<<std::endl;
 //FittedParam.Print();

 ParCov_matr=At*Vym1*A;


 ParCov_matr.Invert(&determ);
 //std::cout<<"MyLinearfitCorr Error Matrix: "<<std::endl;
 //ParCov_matr.Print();

 if(fabs(determ)<0.001)
   std::cout<<"Possible Zero ParCov_matr  Determinant!!!"<<std::endl;  

 //Calcolo chi2
 double chi2=0.;
 TVectorD Residui(sizeArighe);
 Residui=(Yvect-A*FittedParam);
 
 //TVectorD ResiduiT;
 //ResiduiT.Transpose(Residui);
 TVectorD VdotRes(sizeArighe);
  VdotRes=(Vym1*Residui);
 chi2=Residui*VdotRes;

 //----------------------
 // Save the output
 //----------------------
 chi2_=chi2;
 
 if(chi2<0)
   {
     std::cout<<"WARNING: MyLinearfitCorr Chi2: "<<chi2<<std::endl;
     std::cout<<"Residui"<<std::endl;
     Residui.Print();
     std::cout<<"Cov Matrix Vy"<<std::endl;
     Vy.Print();
     std::cout<<"Vym1*Residui"<<std::endl;
     VdotRes.Print();

   }
 par_covariance_matrix=ParCov_matr;
 
 //std::cout<<"MyLinearfitCorr Error Matrix: "<<std::endl;
 //par_covariance_matrix.Print();


 std::vector<double> vect;

 for (unsigned int i=0;i<10;i++)
   vect.push_back(0.);

 vect[0]=FittedParam(0);
 vect[1]=FittedParam(1);
 vect[2]=FittedParam(2);
 vect[3]=FittedParam(3);
 vect[4]=ParCov_matr(0,0);
 vect[5]=ParCov_matr(1,1);
 vect[6]=ParCov_matr(2,2);
 vect[7]=ParCov_matr(3,3);

 double correlabX=ParCov_matr(0,1);
 vect[8]=correlabX;
 double correlabY=ParCov_matr(2,3);
 vect[9]=correlabY;



 return vect;

}

double T2SelectionCutUtils::PhiFitAverage(std::vector<T2Hit> hitvec2)
{

  std::vector<T2Hit> hitvec;
  for(unsigned int jj =0; jj<hitvec2.size(); jj++)
    {   
      if(hitvec2[jj].GetHitClass()!=9)
	hitvec.push_back(hitvec2[jj]);
    }

  unsigned int sizeHitv=hitvec.size();
  unsigned int num0=0;
  unsigned int num360=0;
 
  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      if((hitvec[jj].GetHitPhi()>0.)&&(hitvec[jj].GetHitPhi()<10.))
	num0++;
      if((hitvec[jj].GetHitPhi()>350.)&&(hitvec[jj].GetHitPhi()<360.))
	num360++;
    }


  if((num0>0) && (num360>0)) 
    {
      for(unsigned int jj =0; jj<sizeHitv; jj++)
	{
	  if((hitvec[jj].GetHitPhi()>0.)&&(hitvec[jj].GetHitPhi()<10.))
	    hitvec[jj].SetHitPhi(hitvec[jj].GetHitPhi()+360.);
	}
    }

  float Swphi=0;
  float Sphiwphi=0; 
  std::vector<float> ephi;

  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      ephi.push_back(0.8);
      Swphi += 1.0/ephi[jj]/ephi[jj];
      Sphiwphi+= (hitvec[jj].GetHitPhi())/ephi[jj]/ephi[jj];
    }

  double phim=Sphiwphi/Swphi;
//  double e_phim= 1.0/sqrt(Swphi);

  return phim;
}




double T2SelectionCutUtils::EtaFromAveragesVtxCorr(T1T2Track thetrk,double vtx_x,double vtx_y,double vtx_z)
{

double thetaRHitAverage=0.;double hemis=0.;
  double counter=0.;

  if(thetrk.GetHitT2(1).GetHitClass()!=9)
    hemis=thetrk.GetHitT2(1).GetHitZ()/fabs(thetrk.GetHitT2(1).GetHitZ());
  else
    std::cout<<"Error in T2SelectionCutUtils::EtaFromAverages"<<std::endl;

  for(unsigned int mm=0;mm<thetrk.GetHitEntries();mm++)
    {
      // hitvector.push_back(thetrk.GetHitT2(mm));
       if(thetrk.GetHitT2(mm).GetHitClass()!=9)
	 {

	   double realR= (thetrk.GetHitT2(mm).GetHitX()-vtx_x)*(thetrk.GetHitT2(mm).GetHitX()-vtx_x)+(thetrk.GetHitT2(mm).GetHitY()-vtx_y)*(thetrk.GetHitT2(mm).GetHitY()-vtx_y);
	   realR=sqrt(realR);
	   double realZ=(thetrk.GetHitT2(mm).GetHitZ()-vtx_z);
	   thetaRHitAverage += fabs(realR / realZ);
	   counter=counter+1.0;
	 }
    }
  
  thetaRHitAverage=thetaRHitAverage/counter;
  double TrkEtaFromHitAverages_4Hits = (hemis*(-1.0)*log(thetaRHitAverage/2.0)); //low angle approximation
  
 return TrkEtaFromHitAverages_4Hits;

}

double T2SelectionCutUtils::EtaFromAverages(T1T2Track thetrk)
{

  double thetaRHitAverage=0.;double hemis=0.;
  double counter=0.;

  if(thetrk.GetHitT2(1).GetHitClass()!=9)
    hemis=thetrk.GetHitT2(1).GetHitZ()/fabs(thetrk.GetHitT2(1).GetHitZ());
  else
    std::cout<<"Error in T2SelectionCutUtils::EtaFromAverages"<<std::endl;

  for(unsigned int mm=0;mm<thetrk.GetHitEntries();mm++)
    {
      // hitvector.push_back(thetrk.GetHitT2(mm));
       if(thetrk.GetHitT2(mm).GetHitClass()!=9)
	 {
	   thetaRHitAverage += fabs(thetrk.GetHitT2(mm).GetHitR() / thetrk.GetHitT2(mm).GetHitZ());
	   counter=counter+1.0;
	 }
    }
  
  thetaRHitAverage=thetaRHitAverage/counter;
  double TrkEtaFromHitAverages_4Hits = (hemis*(-1.0)*log(thetaRHitAverage/2.0)); //low angle approximation
  
 return TrkEtaFromHitAverages_4Hits;
}




T1T2Track T2SelectionCutUtils::RPhiFit(std::vector<T2Hit> hitvec2)
{

  T1T2Track trk(2);
  trk.SetDet(2);

  //std::cout<<"T1T2Track BEGIN T2SelectionCutUtils::RPhiFit"<<std::endl;

  std::vector<T2Hit> hitvec; //Only CL1 hits
  std::vector<T2Hit> hitvecR;  //CL1 hits and Vtx
  

  for(unsigned int jj =0; jj<hitvec2.size(); jj++)
   {
     if(hitvec2.at(jj).GetHitClass()==1)
       {
	 hitvec.push_back(hitvec2.at(jj));
       }

     if((hitvec2.at(jj).GetHitClass()==1)||(hitvec2.at(jj).GetHitClass()==9))
       {
	 hitvecR.push_back(hitvec2.at(jj));
       }
   }
  //std::cout<<"Inside MyLinearfitCorr .. init"<<std::endl;

  
  unsigned int  sizeHitv=hitvec.size();
  unsigned int  sizeHitvR=hitvecR.size();
  if(sizeHitv<2)
   {
     std::cout<<" T2SelectionCutUtils::RPhiFit problem: Track with less than 2 Cl1 hits!! Continue the fitting.."<<std::endl;
     hitvec.clear();
     hitvec=hitvec2;
     sizeHitv=hitvec.size();
   }


  std::vector<T2Hit> hitvecraw;
  
  std::vector<T2Hit> doublehitz;
  std::vector<T2Hit> hits0360;


 unsigned int num0=0;
 unsigned int num360=0;
 

for(unsigned int jj =0; jj<sizeHitv; jj++)
  { //std::cout<<"Hitphi "<<hitvec[jj].GetHitPhi()<<std::endl;
    if((hitvec[jj].GetHitPhi()>0.)&&(hitvec[jj].GetHitPhi()<10.))
      num0++;
    if((hitvec[jj].GetHitPhi()>350.)&&(hitvec[jj].GetHitPhi()<360.))
      num360++;
  }


 if((num0>0) && (num360>0)) 
   {
     for(unsigned int jj =0; jj<sizeHitv; jj++)
       {	
	 if((hitvec[jj].GetHitPhi()>0.)&&(hitvec[jj].GetHitPhi()<10.))
	   hitvec[jj].SetHitPhi(hitvec[jj].GetHitPhi()+360.);
       }

     for(unsigned int jj =0; jj<sizeHitvR; jj++)
       {	
	 if((hitvecR[jj].GetHitPhi()>0.)&&(hitvecR[jj].GetHitPhi()<10.))
	   hitvecR[jj].SetHitPhi(hitvecR[jj].GetHitPhi()+360.);
       }
   }



 int hemisphere = (int)(hitvec[0].GetHitZ()/fabs(hitvec[0].GetHitZ()));

 std::vector<float> r;
 std::vector<float> phi;
 std::vector<float> z;
 std::vector<float> er;
 std::vector<float> ephi;
 std::vector<float> ez;


float Sr=0;
float Srz=0;
float Szz_r=0;
float Sz_r=0; 
float S0_r=0; 
float Swphi=0;
float Sphiwphi=0;
 
 float cl1err=(0.12+0.05);
 
 //loop on cl1hits+vtx
  for(unsigned int jj =0; jj<sizeHitvR; jj++)
    {
      
      r.push_back(hitvecR[jj].GetHitR());
      
      z.push_back(hitvecR[jj].GetHitZ());    
    //  er.push_back(hitvecR[jj].GetHitDR());
      if(hitvecR[jj].GetHitClass()==1)
	 er.push_back(cl1err);
      else
	er.push_back(hitvecR[jj].GetHitDR());

      //	er.push_back(0.12+0.05);

      //ephi.push_back(hitvecR[jj].GetHitDPhi());
      
      ez.push_back(hitvecR[jj].GetHitDZ());
      
    
      
      Srz += r[jj]*z[jj]/er[jj]/er[jj];
      Szz_r += z[jj]*z[jj]/er[jj]/er[jj];
      Sz_r += z[jj]/er[jj]/er[jj];
      Sr += r[jj]/er[jj]/er[jj];
      S0_r += 1.0/er[jj]/er[jj];


      if(hitvecR.at(jj).GetHitClass()!=9){
	phi.push_back(hitvecR[jj].GetHitPhi());
	ephi.push_back(0.8);
	Swphi += 1.0/ephi.at(ephi.size()-1)/ephi.at(ephi.size()-1);
	Sphiwphi+=phi.at(phi.size()-1)/ephi.at(ephi.size()-1)/ephi.at(ephi.size()-1);
	//	std::cout<<phi.at(phi.size()-1)<<"  ->"<<Sphiwphi<<std::endl;
	  //Swphi += 1.0/ephi[jj]/ephi[jj];
	  //Sphiwphi+= phi[jj]/ephi[jj]/ephi[jj];
      }
  }


double a_rz = (Srz*S0_r - Sz_r*Sr) / (Szz_r*S0_r - Sz_r*Sz_r);   // angular coefficient
double b_rz = (Sr*Szz_r - Sz_r*Srz) / (Szz_r*S0_r - Sz_r*Sz_r);  // intercept   R=(a_rz)Z + b_rz
double e_a_rz = sqrt( S0_r / (S0_r*Szz_r - Sz_r*Sz_r) );         
double e_b_rz = sqrt( Szz_r / (S0_r*Szz_r - Sz_r*Sz_r) );
double phim=Sphiwphi/Swphi;
double e_phim= 1.0/sqrt(Swphi);

//std::cout<<"->->---->"<<phim<<std::endl;
double covab= - (Sz_r) / (Szz_r*S0_r - Sz_r*Sz_r);

double chi2r = 0.0;
double chi2p = 0.0;
double chi2= 0.0;
double normchi2red=0.0;


 unsigned int cl1count=0;

 for(unsigned int jj =0; jj<sizeHitvR; jj++)
   {
     if(hitvecR[jj].GetHitClass()==1)
       {
	 chi2r += (a_rz*hitvecR[jj].GetHitZ()+b_rz - hitvecR[jj].GetHitR())*(a_rz*hitvecR[jj].GetHitZ()+b_rz - hitvecR[jj].GetHitR())/cl1err/cl1err;
	 cl1count++;
       }
     else
       chi2r += (a_rz*hitvecR[jj].GetHitZ()+b_rz - hitvecR[jj].GetHitR())*(a_rz*hitvecR[jj].GetHitZ()+b_rz - hitvecR[jj].GetHitR())/hitvecR[jj].GetHitDR()/hitvecR[jj].GetHitDR();
       
   }



  for(unsigned int jjj =0; jjj<sizeHitv; jjj++){          
      
      chi2p += (phim - phi[jjj])*(phim - phi[jjj])/ephi[jjj]/ephi[jjj];
    }
  
  if(cl1count>=2) 
    {      
      chi2=  (chi2p+chi2r)/2.0;
      normchi2red=(chi2p/((double)(cl1count)-1)+ chi2r/((double)(cl1count)-2))/2.0 ;
    }
  else
    {
      chi2=1.0;
      normchi2red=1.0;
      std::cout<<"Error: T2SelectionCutUtils RPhiFit without at least 2 hits"<<std::endl;
    }


  if(phim>360.0)                // This could happen when num0>0 AND num360>0
    phim=phim-360.0;
    

    phim=phim*3.14159265/180.0;
    e_phim=e_phim*3.14159265/180.0;

    //std::cout<<"rphifit: "<<phim*180/3.14159265<<std::endl;
    

    TVectorD vect(4);
    for(int oo=0; oo<4; oo++)
      vect[oo]=0.;
   
  
    vect[0]=phim;    // phi intercept
    vect[1]=0.;      // phi ang. coeff.
    vect[2]=b_rz;    // r intercept
    vect[3]=a_rz;    // r ang. coeff.
   

  
    TMatrixD mat(4,4);
    for(unsigned int oo=0; oo<4; oo++)
      for(unsigned int ooo=0;ooo<4; ooo++)
	mat[oo][ooo]=0.;


      mat[0][0]=e_phim*e_phim;
      mat[0][1]=0.;
      mat[0][2]=0.;
      mat[0][3]=0.;
      mat[1][0]=0.;
      mat[1][1]=0.;
      mat[1][2]=0.;
      mat[1][3]=0.;
      mat[2][0]=0.;
      mat[2][1]=0.;
      mat[2][2]=e_b_rz*e_b_rz;
      mat[2][3]=covab;
      mat[3][0]=0.;
      mat[3][1]=0.;
      mat[3][2]=covab;
      mat[3][3]=e_a_rz*e_a_rz;
    


//   double thetheta =atan(a_rz);
//   double theeta = fabs(-log(tan(thetheta/2.))) * (double)hemisphere;

   //std::cout<<"calc phim "<<phim*180.0/3.14159265<<" "<<hitvec2[0].GetHitPhi()<<" "<<hitvec[0].GetHitPhi()<<std::endl;
   
    T1T2Track fittedtrack(vect,mat,a_rz, b_rz, phim, e_a_rz, e_b_rz, e_phim, chi2,chi2r,chi2p, normchi2red, hemisphere,2);   
    for(unsigned int jj =0; jj<hitvec2.size(); jj++)
      {
	//fittedtrack.AddHit(hitvec[jj]);
	fittedtrack.AddHit(hitvec2[jj]);
      }
    


    return fittedtrack; 


}








 T1T2Track T2SelectionCutUtils::TrackFromHits(bool DoXYTrkFit, std::vector<T2Hit> hitvec)
{

  T1T2Track toreturn;
  if(hitvec.size()<=2)
    std::cout<<"TrackFromHits (hitvec.size()<=2) "<<hitvec.size()<<std::endl;

  if(DoXYTrkFit)
    {
      //The commented part belove has been used until 21 August 2010. Then I decided to put fitTrackspecialXY
      //because chi2X-chi2Y was not saved correctly.

      /*
      // std::cout<<"a0 "<<std::endl;      
       TMatrixD covmat(4,4);
       // std::cout<<"a1 "<<std::endl;
       covmat.Zero();
       double chi2corr;	      
       
       double chi2X= 1.; //to be properly implemented 
       double chi2Y= 1.; //to be properly implemented 
       double chi2 = 1.; //to be properly implemented 
       //std::cout<<"a2 "<<std::endl;
       int hemisphere = (int)(hitvec[0].GetHitZ()/fabs(hitvec[0].GetHitZ()));	             
       std::vector<double> corrFit=MyLinearfitXYSimple(hitvec);
       TVectorD vect(4);
       vect[0] = corrFit[1];//b_xz;
       vect[1] = corrFit[3];//b_yz;
       vect[2] = corrFit[0];//a_xz;
       vect[3] = corrFit[2];//a_yz;
       covmat(0,0)=corrFit[5]*corrFit[5];
       covmat(1,1)=corrFit[7]*corrFit[7];
       covmat(2,2)=corrFit[4]*corrFit[4];
       covmat(3,3)=corrFit[6]*corrFit[6];

       // std::cout<<"a5 "<<std::endl;
       //( const TVectorD & track_params_vector, const TMatrixD &par_covariance_matrix, double chiSquared, int hemi, int det)
       // T1T2Track fittedtrack(vect,covmat,chi2,chi2X,chi2Y,hemisphere,2);
       T1T2Track fittedtrack(vect,covmat,chi2,hemisphere,2);
       for(unsigned int i=0;i<hitvec.size();i++)
         fittedtrack.AddHit(hitvec.at(i));
       //std::cout<<"a6 "<<std::endl;
       */
      //  std::cout<<"TrackFromHits (hitvec.size()<=2) "<<hitvec.size()<<std::endl;

      T1T2Track fittedtrack =  fitTrackspecialXY(hitvec, false);

       toreturn=fittedtrack;
    }
  else
    {
       toreturn=RPhiFit(hitvec);
       //for(unsigned int i=0;i<hitvec.size();i++)
       // toreturn.AddHit(hitvec.at(i));
      //T1T2Track fittedtrack(vect,mat,a_rz, b_rz, phim, e_a_rz, e_b_rz, e_phim, chi2,chi2r,chi2p, normchi2red, hemisphere,2);  
    }

  


  return toreturn;
}



  
T1T2Track T2SelectionCutUtils::TrackXYFromHitsExcludingOneplane(unsigned int RawIdToExclude, T1T2Track origtrk)
{

  std::vector<T2Hit>  hitv;
  for(unsigned int m=0;m<origtrk.GetHitEntries();m++)
    {
      hitv.push_back(origtrk.GetHitT2(m));
    }

  std::vector<T2Hit> hitv2=GiveMeTheHitForResidual(RawIdToExclude, hitv);
  T1T2Track thetrk=TrackFromHits(true, hitv2);
  return thetrk;
  
}



std::vector<T2Hit> T2SelectionCutUtils::GiveMeTheHitForResidual(unsigned int RawIdToExclude, std::vector<T2Hit>  hitv)
{

  std::vector<T2Hit> returnedVect;
  for(unsigned int i=0;i<hitv.size();i++)
    {
      if(hitv.at(i).GetHitDetRawId()!=RawIdToExclude)
	returnedVect.push_back(hitv.at(i));
    }

  //std::cout<<"Sizebef: "<<hitv.size()<<"Sizeaft: "<<returnedVect.size()<<"  |   ";
  return returnedVect;
}




double T2SelectionCutUtils::sigmaZfunctE10(double eta)
{
  double theDZ=0.0; 
  double nepero=10.0;
 

 if(fabs(eta)<5.37)
     theDZ=3200.0;

  if(fabs(eta)>6.25)
     theDZ=4400.0;


   if((fabs(eta)>5.37)&&(fabs(eta)<6.25))
     {     
       theDZ=1.09980*pow(nepero,5.0) -3.84187*pow(nepero,4.0)*eta+  3.44411*pow(nepero,3.0)*eta*eta;
     }

   //if((fabs(eta)>5.42)&&(fabs(eta)<5.78))
   if((fabs(eta)>5.5)&&(fabs(eta)<5.75))  
     theDZ=3300.0;
 

  return theDZ;
}

bool T2SelectionCutUtils::IsinT2(double eta){
  bool flag=false;

  if ((fabs(eta)>T2_TrkEtamin)&&(fabs(eta)<T2_TrkEtaMAX))
    flag=true;

  //if(!flag)
  //std::cout<<"IsinT2 failed"<<std::endl;

  return flag;
}



bool T2SelectionCutUtils::T2TrkMultCond(unsigned int trkmultiplicity){
  bool flag=false;

  if(((int)trkmultiplicity >= T2_trkMultimin)&&((int)trkmultiplicity <= T2_trkMultiMAX))
    flag=true;

  //if(!flag)
  //std::cout<<"T2TrkMultCond failed: trk mult="<<trkmultiplicity<<"  min="<<T2_trkMultimin<<" max="<<T2_trkMultimin<<std::endl;

  return flag;
}



double T2SelectionCutUtils::Chi2R_Prob_Calculator(T1T2Track trk){
  double chi2RProb=0.;
  T1T2Track rztrk=TrackFromHits(false, trk.GetT2TrackHits());
  chi2RProb=TMath::Prob(rztrk.ChiSquaredR(),(rztrk.GetHitEntries()-2));  
  return chi2RProb;
}



bool T2SelectionCutUtils::ChiCutCond(T1T2Track trk, bool xyfitused,double prob1, double prob2){
  bool flag=false;
  
  double chiPhiProb=0.,chiRProb=0.;

  chiRProb=TMath::Prob(trk.ChiSquaredR(),(trk.GetHitEntries()-2));
  chiPhiProb=TMath::Prob(trk.ChiSquaredPhi(),(trk.GetHitEntries()-1)); 
 
  if((chiRProb>prob1)&&(chiPhiProb>prob2))  
    flag=true;

  if(xyfitused)
    {
      flag=false;
      //x=phi;y=R
      double chiXProb=0.,chiYProb=0.;
      chiXProb=TMath::Prob(trk.ChiSquaredX(),(trk.GetHitEntries()-2));
      chiYProb=TMath::Prob(trk.ChiSquaredY(),(trk.GetHitEntries()-2));
      if((chiXProb>=prob1)&&(chiYProb>=prob1)) 
	flag=true;
    }
  

  // if(!flag)
  //std::cout<<"ChiCutCond failed"<<chiRProb<<" "<<chiPhiProb<<" "<<trk.ChiSquaredR()<<trk.GetHitEntries()<<std::endl;

  return flag;
}



//Compare trk eta with the eta expected from the position of the first hit
bool EtaTrkCompatibleWithImpact(T1T2Track trk, double cut=0.1)
{
  bool flag=false;
  double etatrk=trk.Eta();
  if(trk.GetHitT2(1).GetHitClass()!=9)
    {
      double hittantheta=fabs(trk.GetHitT2(1).GetHitR()/trk.GetHitT2(1).GetHitZ());
      double hemi=trk.GetHitT2(1).GetHitZ()/fabs(trk.GetHitT2(1).GetHitZ());
      double hiteta=fabs(-log(hittantheta/2.0))*hemi; //Ok, it is approximated...
      if(fabs(hiteta-etatrk)<cut)
	flag=true;
    }
  else
    std::cout<<"EtaTrkCompatibleWithImpact Warning"<<std::endl;
    
  return flag;
}







bool T2SelectionCutUtils::ChiCutCond(T1T2Track trk){
  bool flag=false;
  
  double chiPhiProb=0.,chiRProb=0.;

  chiRProb=TMath::Prob(trk.ChiSquaredR(),(trk.GetHitEntries()-2));
  chiPhiProb=TMath::Prob(trk.ChiSquaredPhi(),(trk.GetHitEntries()-1)); 
 
  if((chiPhiProb>T2_PhiChiProbCut)&&(chiRProb>T2_RChiProbCut))  
    flag=true;

  if(XYFitUsed)
    {
      flag=false;
      //x=phi;y=R
      double chiXProb=0.,chiYProb=0.;
      chiXProb=TMath::Prob(trk.ChiSquaredX(),(trk.GetHitEntries()-2));
      chiYProb=TMath::Prob(trk.ChiSquaredY(),(trk.GetHitEntries()-2));
      if((chiXProb>=T2_PhiChiProbCut)&&(chiYProb>=T2_RChiProbCut)) 
	flag=true;
      //std::cout<<chiXProb<<" "<<chiYProb<<" chi2X: "<<trk.ChiSquaredX()<<std::endl;
    }
  

  //if(!flag)
  //std::cout<<"ChiCutCond failed"<<std::endl;

  return flag;
}


bool T2SelectionCutUtils::TrkInQuarter(T1T2Track trk,unsigned int selectedquarter)
{

  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo;
  bool t2trkcutPassed=false;

  unsigned int startindex;
  if(trk.GetHitT2(0).GetHitClass()!=9)
    startindex=0;
  else
    startindex=1;
  
  if((trk.GetHitEntries()>=(3+startindex))&&(trk.GetHitEntries()<=(10+startindex))) //Otherwise is an overlap
    {	 
      
	
      planeinfo=conv.GetT2Info(trk.GetHitT2(0+startindex).GetHitDetRawId());
      unsigned int symbh0=planeinfo.arm*2+planeinfo.ht; 
      planeinfo=conv.GetT2Info(trk.GetHitT2(1+startindex).GetHitDetRawId());  
      unsigned int symbh1=planeinfo.arm*2+planeinfo.ht; 
      planeinfo=conv.GetT2Info(trk.GetHitT2(2+startindex).GetHitDetRawId());	  
      unsigned int symbh2=planeinfo.arm*2+planeinfo.ht;  

      if(symbh0==selectedquarter)
	if(symbh1==selectedquarter)
	  if(symbh2==selectedquarter)
	    t2trkcutPassed=true;  
	
      //else
      //	std::cout<<"Error in T2SelectionCutUtils::TrkInQuarter"<<std::endl;
    }
 
  //if(!t2trkcutPassed)
  //std::cout<<"TrkInQuarter failed"<<std::endl;


  return t2trkcutPassed;

}


bool T2SelectionCutUtils::TrkInQuarter(T1T2Track trk)
{

  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo;
  bool t2trkcutPassed=false;
  unsigned int startindex;
  if(trk.GetHitT2(0).GetHitClass()!=9)
    startindex=0;
  else
    startindex=1;
  if(trk.GetHitEntries()>=(3+startindex))
    {	 
      planeinfo=conv.GetT2Info(trk.GetHitT2(0+startindex).GetHitDetRawId());
      unsigned int symbh0=planeinfo.arm*2+planeinfo.ht; 
      planeinfo=conv.GetT2Info(trk.GetHitT2(1+startindex).GetHitDetRawId());  
      unsigned int symbh1=planeinfo.arm*2+planeinfo.ht; 
      planeinfo=conv.GetT2Info(trk.GetHitT2(2+startindex).GetHitDetRawId());	  
      unsigned int symbh2=planeinfo.arm*2+planeinfo.ht;    

      //     std::cout<<"symbh012  "<<symbh0<<"-"<<symbh1<<"-"<<symbh2<<std::endl;
      if(std::find(T2_QuarterUsed.begin(), T2_QuarterUsed.end(), symbh0)!=T2_QuarterUsed.end())
	if(std::find(T2_QuarterUsed.begin(), T2_QuarterUsed.end(), symbh1)!=T2_QuarterUsed.end())
	  if(std::find(T2_QuarterUsed.begin(), T2_QuarterUsed.end(), symbh2)!=T2_QuarterUsed.end())
	    t2trkcutPassed=true;
    }
  
  //if(!t2trkcutPassed)
  //std::cout<<"TrkInQuarter failed"<<std::endl;
  return t2trkcutPassed;
}

bool T2SelectionCutUtils::TrkInQuarterTight(T1T2Track trk,unsigned int selectedquarter)
{
  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo;
  bool t2trkcutPassed=true;
  unsigned int symb;
  for(unsigned int i=0;i<trk.GetHitEntries();i++)
    {
      if(trk.GetHitT2(i).GetHitClass()!=9){
      planeinfo=conv.GetT2Info(trk.GetHitT2(i).GetHitDetRawId());
      symb=planeinfo.arm*2+planeinfo.ht; 
      if(symb!=selectedquarter)
	 t2trkcutPassed=false;
      }
    }
  
  
  //if(!t2trkcutPassed)
  //std::cout<<"TrkInQuarter failed"<<std::endl;
  return t2trkcutPassed;
}

//start from here
bool T2SelectionCutUtils::TrkInQuarterTight(T1T2Track trk)
{

  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo;
  bool t2trkcutPassed=true;
  unsigned int symb; 

  for(unsigned int i=0;i<trk.GetHitEntries();i++)
    {
      if(trk.GetHitT2(i).GetHitClass()!=9){
      planeinfo=conv.GetT2Info(trk.GetHitT2(i).GetHitDetRawId());
      symb=planeinfo.arm*2+planeinfo.ht; 
       if(std::find(T2_QuarterUsed.begin(), T2_QuarterUsed.end(), symb)==T2_QuarterUsed.end())
	 t2trkcutPassed=false;
      }
    }
  
  
  //if(!t2trkcutPassed)
  //std::cout<<"TrkInQuarter failed"<<std::endl;
  return t2trkcutPassed;
}






bool T2SelectionCutUtils::TrkAlsoInQuarter(T1T2Track trk)
{

  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo;
  bool t2trkcutPassed=false;
  unsigned int startindex;
  if(trk.GetHitT2(0).GetHitClass()!=9)
    startindex=0;
  else
    startindex=1;

  if(trk.GetHitEntries()>=(3+startindex))
    {	 
      planeinfo=conv.GetT2Info(trk.GetHitT2(0+startindex).GetHitDetRawId());
      unsigned int symbh0=planeinfo.arm*2+planeinfo.ht; 
      planeinfo=conv.GetT2Info(trk.GetHitT2(1+startindex).GetHitDetRawId());  
      unsigned int symbh1=planeinfo.arm*2+planeinfo.ht; 
      planeinfo=conv.GetT2Info(trk.GetHitT2(2+startindex).GetHitDetRawId());	  
      unsigned int symbh2=planeinfo.arm*2+planeinfo.ht;    

      //     std::cout<<"symbh012  "<<symbh0<<"-"<<symbh1<<"-"<<symbh2<<std::endl;
      if((std::find(T2_QuarterUsed.begin(), T2_QuarterUsed.end(), symbh0)!=T2_QuarterUsed.end())||(std::find(T2_QuarterUsed.begin(), T2_QuarterUsed.end(), symbh1)!=T2_QuarterUsed.end())||(std::find(T2_QuarterUsed.begin(), T2_QuarterUsed.end(), symbh2)!=T2_QuarterUsed.end())) 	    t2trkcutPassed=true;

      
      if(trk.GetHitEntries()>=(4+startindex))
	{
	  planeinfo=conv.GetT2Info(trk.GetHitT2(3+startindex).GetHitDetRawId());	  
	  unsigned int symbh4=planeinfo.arm*2+planeinfo.ht;    
	  if(std::find(T2_QuarterUsed.begin(), T2_QuarterUsed.end(), symbh4)!=T2_QuarterUsed.end())
	     t2trkcutPassed=true;
	}

    }
  
  //if(!t2trkcutPassed)
  //std::cout<<"TrkInQuarter failed"<<std::endl;
  return t2trkcutPassed;
}


bool T2SelectionCutUtils::TrkAlsoInQuarter(T1T2Track trk, std::vector<int>  QuarterWanted)
{
 T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo;
  bool t2trkcutPassed=false;
  unsigned int startindex;
  if(trk.GetHitT2(0).GetHitClass()!=9)
    startindex=0;
  else
    startindex=1;

  if(trk.GetHitEntries()>=(3+startindex))
    {	 
      planeinfo=conv.GetT2Info(trk.GetHitT2(0+startindex).GetHitDetRawId());
      unsigned int symbh0=planeinfo.arm*2+planeinfo.ht; 
      planeinfo=conv.GetT2Info(trk.GetHitT2(1+startindex).GetHitDetRawId());  
      unsigned int symbh1=planeinfo.arm*2+planeinfo.ht; 
      planeinfo=conv.GetT2Info(trk.GetHitT2(2+startindex).GetHitDetRawId());	  
      unsigned int symbh2=planeinfo.arm*2+planeinfo.ht;    

      //     std::cout<<"symbh012  "<<symbh0<<"-"<<symbh1<<"-"<<symbh2<<std::endl;
      if((std::find(QuarterWanted.begin(), QuarterWanted.end(), symbh0)!=QuarterWanted.end())||(std::find(QuarterWanted.begin(), QuarterWanted.end(), symbh1)!=QuarterWanted.end())||(std::find(QuarterWanted.begin(), QuarterWanted.end(), symbh2)!=QuarterWanted.end())) 	   
	t2trkcutPassed=true;

      
      if(trk.GetHitEntries()>=(4+startindex))
	{
	  planeinfo=conv.GetT2Info(trk.GetHitT2(3+startindex).GetHitDetRawId());	  
	  unsigned int symbh4=planeinfo.arm*2+planeinfo.ht;    
	  if(std::find(QuarterWanted.begin(), QuarterWanted.end(), symbh4)!=QuarterWanted.end())
	     t2trkcutPassed=true;
	}

    }
  
  //if(!t2trkcutPassed)
  //std::cout<<"TrkInQuarter failed"<<std::endl;
  return t2trkcutPassed;


}





bool T2SelectionCutUtils::HitVectorInAQuarter(std::vector<T2Hit> hitvect, unsigned int selectedquarter){
  
  bool t2trkcutPassed=true; 
  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo;
 

  for(unsigned int i=0;i<hitvect.size();i++)
    {      
      if(hitvect.at(i).GetHitClass()!=9)
	{
	  planeinfo=conv.GetT2Info(hitvect.at(i).GetHitDetRawId());
	  unsigned int symbh0=planeinfo.arm*2+planeinfo.ht; 
	  if(symbh0!=selectedquarter)
	    t2trkcutPassed=false;
	}
    }

   if(hitvect.size()<3)
     t2trkcutPassed=false;
   else
     if(((hitvect.at(0).GetHitClass()==9)&&(hitvect.size()==1)))
    t2trkcutPassed=false;

  return t2trkcutPassed;
}




bool T2SelectionCutUtils::ThetaAsPrimary(T1T2Track trk, bool usexytrack)
{
  bool flag=false;

  // double thephi=trk.Phi();
  //thephi=thephi*180.0/3.14159;
  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo;
  
  double signoftrkX=0.0;
  double signoftrkY=0.0;

  double arm=0.0;
  //Look at the sign of the trk hit X
//  int hemis=0;
  for(unsigned int n=0;n<trk.GetHitEntries();n++)
    {
      if(trk.GetHitT2(n).GetHitClass()!=9)
	{
	  planeinfo=conv.GetT2Info(trk.GetHitT2(n).GetHitDetRawId());
	  arm+=planeinfo.arm;
	  signoftrkX+=trk.GetHitT2(n).GetHitX();
	  signoftrkY+=trk.GetHitT2(n).GetHitY();
	  trk.GetHitT2(n).GetHitZ();
	  fabs(trk.GetHitT2(n).GetHitZ());
	}
    }   
  
  double armmult=0.;
  if((arm/double(trk.GetHitEntries()))>0.5)
    {
      arm=1.;
      armmult=-1.0;
    }
  else
    {
      arm=0.;
      armmult=1.0; 	
    }
  
  /*
  if(hemis!=arm)
    std::cout<<
  */

  if(signoftrkX<0)
    signoftrkX=-1.0;
  else
    signoftrkX=1.0;

  if(signoftrkY<0)
    signoftrkY=-1.0;
  else
    signoftrkY=1.0;
 
  //if(arm!=0.)
  // std::cout<<"WARN !! "<<armmult<<std::endl;
  if(usexytrack)
    {
      if(((armmult*trk.GetTx()*signoftrkX)>=0.0)&&((armmult*trk.GetTy()*signoftrkY)>=0.0))	
	flag=true;
      else
	{
//modification, 26 june
	  if((((armmult*trk.GetTx()*signoftrkX)<0.0))||(fabs(trk.GetTx())<IgnoredSmallAnglePrimarySlopeCut))	     
	    flag=true;


	   //if(((armmult*trk.GetTy()*signoftrkY)<0.0))
	     //if(fabs(trk.GetTy())<IgnoredSmallAnglePrimarySlopeCut)
	       //flag=true;

	  if((((armmult*trk.GetTy()*signoftrkY)<0.0))||(fabs(trk.GetTy())<IgnoredSmallAnglePrimarySlopeCut))
	    flag=true;
	     //if(fabs(trk.GetTy())<IgnoredSmallAnglePrimarySlopeCut)
	       //flag=true;


	  //	  std::cout<<"arm | "<<"Tx | "<<"SignX | "<<"Ty | "<<"SignY "<<" | Phi     :->   "<<arm<<" | "<<trk.GetTx()<<" | "<<signoftrkX<<" | "<<trk.GetTy()<<" | "<<signoftrkY<<" | "<<trk.Phi()*180.0/3.14159<<std::endl;
	}
    }
  else
    {
      if(trk.Theta()>=0)
	flag=true;
      else
	if(fabs(trk.Theta())<IgnoredSmallAnglePrimarySlopeCut)
	  flag=true;
    }


 

  // if(!flag)
  //std::cout<<"ThetaAsPrimary failed"<<std::endl;


  return flag;
}


bool T2SelectionCutUtils::DZCutPassed(T1T2Track trk)
{

  bool flag=false;
   double DZgood=1000.;
   DZgood=sigmaZfunctE10(fabs(trk.Eta()));
   DZgood=DZgood*T2_DZMultiplicator;

   double firsthitZ=trk.GetHitT2(0).GetHitZ();

   double C0=(trk.GetHitT2(0).GetHitX()*trk.GetHitT2(0).GetHitX()+trk.GetHitT2(0).GetHitY()*trk.GetHitT2(0).GetHitY())/(trk.GetTx()*trk.GetHitT2(0).GetHitX()+trk.GetTy()*trk.GetHitT2(0).GetHitY());
   double Z0impact=trk.GetHitT2(0).GetHitZ()-C0;

   


   if(trk.GetHitT2(0).GetHitClass()==9)
     firsthitZ=trk.GetHitT2(1).GetHitZ();


   //old strat
   /*
     if(fabs(trk.Z_at_Rmin())<DZgood)
     {
     flag=true;
     }
   


     if(((Z0impact*firsthitZ)>0)&&(fabs(Z0impact)<5000.))
     flag=true;
   */

   if(((Z0impact*firsthitZ)<0)||(fabs(Z0impact)<6500.))
      flag=true;

   /*
   if( (((trk.Z_at_Rmin()*firsthitZ)>0)&&(fabs(trk.Z_at_Rmin())<5000.))||
       (((trk.Z_at_Rmin()*firsthitZ)<0)&&(fabs(trk.bx_)<25.)&&(fabs(trk.by_)<25.))) //12/9/2010 Cut Put to 60-> gave Z Asym bad shape DN/deta!
     flag=true;
   */


   //if(!flag)
   //std::cout<<"DZCutPassed failed: "<<(trk.Z_at_Rmin())<<"Max:"<<DZgood<<std::endl;
   
   return flag;
}



bool T2SelectionCutUtils::ThetaCutPassed(T1T2Track trk)
{

  bool flag=false;
  if((fabs(trk.GetTx())<0.1)&&(fabs(trk.GetTy())<0.1))
    flag=true;
  /*  
  if(!flag)
   std::cout<<"ThetaCutPassed failed: "<<(trk.GetTx())<<"  "<<(trk.GetTy())<<std::endl;
  */
   return flag;
}



bool T2SelectionCutUtils::AcceptThisT2TrackWhatEverQuarter(T1T2Track trk)
{
  bool t2trkcutPassed=false;
  unsigned int cl1hitMult=0;
  for(unsigned int i=0;i<trk.GetHitEntries();i++)
    {
      if(trk.GetHitT2(i).GetHitClass()==1)
	cl1hitMult++;
    }


  if(ThetaCutPassed(trk)==true)   
    if(DZCutPassed(trk)==true)
      if(IsinT2(trk.Eta())==true)
	if(T2TrkMultCond(cl1hitMult)==true)    //if(T2TrkMultCond(trk.GetHitEntries())==true) 
	  if(ChiCutCond(trk)==true)  
	      if(ThetaAsPrimary(trk,XYFitUsed))
		t2trkcutPassed=true;

 return t2trkcutPassed;
}



bool T2SelectionCutUtils::AcceptThisT2Track(T1T2Track trk)
{
  bool t2trkcutPassed=false;
  unsigned int cl1hitMult=0; 
  //std::cout<<"AcceptThisT2Track-A"<<std::endl;
  for(unsigned int i=0;i<trk.GetHitEntries();i++)
    {
      if(trk.GetHitT2(i).GetHitClass()==1)
	cl1hitMult++;
    }
  //std::cout<<"AcceptThisT2Track-B"<<std::endl;

 if(DZCutPassed(trk)==true)
   if(ThetaCutPassed(trk)==true)   
      if(IsinT2(trk.Eta())==true)
	if(T2TrkMultCond(cl1hitMult)==true)    //	if(T2TrkMultCond(trk.GetHitEntries())==true)
	  if(ChiCutCond(trk)==true)  
	    if(TrkInQuarter(trk)==true)	      
	      if(ThetaAsPrimary(trk,XYFitUsed))
		t2trkcutPassed=true;
	
 /*
  if(!t2trkcutPassed)
    {
      std::cout<<std::endl;
      std::cout<<"AcceptThisT2Track failed because of:"<<std::endl;
      if(!ThetaCutPassed(trk))
	std::cout<<"ThetaCutPassed(trk)"<<std::endl;

      if(!DZCutPassed(trk))
	std::cout<<"DZCutPassed(trk)"<<std::endl;

      if(!IsinT2(trk.Eta()))
	std::cout<<"IsinT2(trk)"<<std::endl;

      if(!T2TrkMultCond(trk.GetHitEntries()))
	std::cout<<"T2TrkMultCond(trk) "<<trk.GetHitEntries()<<std::endl;

      if(!(ChiCutCond(trk)))
	std::cout<<"ChiCutCond(trk)"<<std::endl;

      if(!TrkInQuarter(trk))
	std::cout<<"TrkInQuarter(trk)"<<std::endl;

      if(!(ThetaAsPrimary(trk,XYFitUsed)))
	std::cout<<"ThetaAsPrimary(trk)"<<std::endl;
      std::cout<<std::endl;
      std::cout<<std::endl;
      }
 */
  
  return t2trkcutPassed;
}




long double T2SelectionCutUtils::CartesianDistance(T2SelectionCutUtils::point3d p1,T2SelectionCutUtils::point3d p2)
{
  long double distance;
  long double d2=(p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z);
  distance=sqrt(d2);  
  return distance;
}


T2SelectionCutUtils::TracksSeparation T2SelectionCutUtils::TrkZ0(T1T2Track t1, bool usesXYtracking){

   T2SelectionCutUtils::TracksSeparation separ;
  long double a1x ;
  long double b1x ;
  long double a1y;
  long double b1y ;

  long double a2x ;
  long double b2x ;
  long double a2y ;
  long double b2y ;


  if(usesXYtracking)
    {
      	a1x = t1.GetTx(); 
	b1x = t1.X0(); 
	a1y = t1.GetTy(); 
	b1y = t1.Y0();
	
	a2x = 0.; 
	b2x = 0.; 
	a2y = 0.; 
	b2y = 0.; 
    }
  else
    {
      a1x = t1.GetTy() * cos( t1.Phi() ); 
      b1x = t1.GetTx() * cos( t1.Phi() ); 
      a1y = t1.GetTy() * sin( t1.Phi() ); 
      b1y = t1.GetTx() * sin( t1.Phi() ); 
 
      a2x = 0.; 
      b2x = 0.; 
      a2y = 0.; 
      b2y = 0.; 
      
    }
  

  long double A=  a1x*a1x + a1y*a1y + 1;
  long double B=  a1x*a2x + a1y*a2y + 1;
  long double C=  a2x*a2x + a2y*a2y + 1;
  long double D=  a1x*(b1x-b2x)+ a1y*(b1y-b2y);
  long double E=  a2x*(b1x-b2x)+ a2y*(b1y-b2y);

  long double Zc1= (B*E-C*D)/(A*C-B*B);
  long double Zc2= (A*C-B*D)/(A*C-B*B);
 
  T2SelectionCutUtils::point3d P1;
  P1.x = b1x + a1x*Zc1;
  P1.y = b1y + a1y*Zc1;
  P1.z = Zc1;

  T2SelectionCutUtils::point3d P2;
  P2.x = b2x + a2x*Zc2;
  P2.y = b2y + a2y*Zc2;
  P2.z = Zc2;

//  T2SelectionCutUtils::point3d midpoint;
//  midpoint.x = (b1x + b2x + a1x*Zc1 + a2x*Zc2)/2.0;
//  midpoint.y = (b1y + b2y + a1y*Zc1 + a2y*Zc2)/2.0;
//  midpoint.z = (Zc1+Zc2)/2.0;

  long double distanza=CartesianDistance(P1,P2);

  separ.distance=distanza;
  separ.midPoint=P1;

 
  

  return separ;  


}







T2SelectionCutUtils::TracksSeparation T2SelectionCutUtils::TwoTracksSeparation(T1T2Track t1, T1T2Track t2, bool usesXYtracking){
  T2SelectionCutUtils::TracksSeparation separ;
  long double a1x ;
  long double b1x ;
  long double a1y;
  long double b1y ;

  long double a2x ;
  long double b2x ;
  long double a2y ;
  long double b2y ;


  if(usesXYtracking)
    {
      	a1x = t1.GetTx(); 
	b1x = t1.X0(); 
	a1y = t1.GetTy(); 
	b1y = t1.Y0();
	
	a2x = t2.GetTx(); 
	b2x = t2.X0(); 
	a2y = t2.GetTy(); 
	b2y = t2.Y0(); 
    }
  else
    {
      a1x = t1.GetTy() * cos( t1.Phi() ); 
      b1x = t1.GetTx() * cos( t1.Phi() ); 
      a1y = t1.GetTy() * sin( t1.Phi() ); 
      b1y = t1.GetTx() * sin( t1.Phi() ); 
 
      a2x = t2.GetTy() * cos( t2.Phi() ); 
      b2x = t2.GetTx() * cos( t2.Phi() ); 
      a2y = t2.GetTy() * sin( t2.Phi() ); 
      b2y = t2.GetTx() * sin( t2.Phi() );
    }
  

  long double A=  a1x*a1x + a1y*a1y + 1;
  long double B=  a1x*a2x + a1y*a2y + 1;
  long double C=  a2x*a2x + a2y*a2y + 1;
  long double D=  a1x*(b1x-b2x)+ a1y*(b1y-b2y);
  long double E=  a2x*(b1x-b2x)+ a2y*(b1y-b2y);

  long double Zc1= (B*E-C*D)/(A*C-B*B);
  long double Zc2= (A*E-B*D)/(A*C-B*B);
 
  T2SelectionCutUtils::point3d P1;
  P1.x = b1x + a1x*Zc1;
  P1.y = b1y + a1y*Zc1;
  P1.z = Zc1;

  T2SelectionCutUtils::point3d P2;
  P2.x = b2x + a2x*Zc2;
  P2.y = b2y + a2y*Zc2;
  P2.z = Zc2;

  T2SelectionCutUtils::point3d midpoint;
  midpoint.x = (b1x + b2x + a1x*Zc1 + a2x*Zc2)/2.0;
  midpoint.y = (b1y + b2y + a1y*Zc1 + a2y*Zc2)/2.0;
  midpoint.z = (Zc1+Zc2)/2.0;

  long double distanza=CartesianDistance(P1,P2);

  separ.distance=distanza;
  separ.midPoint=midpoint;
  /*
  std::cout<<"P1 xyz: "<<P1.x<<" | "<<P1.y<<" | "<<P1.z<<std::endl;
  std::cout<<"P2 xyz: "<<P2.x<<" | "<<P2.y<<" | "<<P2.z<<std::endl;
  std::cout<<"distanza = "<<distanza<<std::endl;
  std::cout<<"MidPoint xyz: "<<midpoint.x<<" | "<<midpoint.y<<" | "<<midpoint.z<<std::endl;
  */
  return separ;    
}




/*T2SelectionCutUtils::CandidateVtxWithTrks*/  
//void T2SelectionCutUtils::FindPrimaryVtx(std::vector<T1T2Track> startingtracks2,bool usexytrack, bool UseAnyQuarter, double Maxbxtrk=100.,double Maxbytrk=100.){
//void T2SelectionCutUtils::FindPrimaryVtx(std::vector<T1T2Track> startingtracks2,bool usexytrack, bool UseAnyQuarter=false, double Maxbxtrk=100.,double Maxbytrk=100.){  
void T2SelectionCutUtils::FindPrimaryVtx(std::vector<T1T2Track> startingtracks2,bool usexytrack, bool UseAnyQuarter, double Maxbxtrk, double Maxbytrk){  

  T2SelectionCutUtils::CandidateVtxWithTrks thevtx;
  T2SelectionCutUtils::TracksSeparation actualtrksep;
  std::vector<std::vector<T1T2Track> > allTracksclusters;
  std::vector<T1T2Track> startingtracks;
  T2GeometryUtil conv;
  
  storeVtx.VtxTrks.clear();
  storeVtx.VtxPos.z=-10000.;
  storeVtx.VtxPos.x=-10000.;
  storeVtx.VtxPos.y=-10000.;

  double bx_,by_;
  // std::cout<<"1-Vertex function will work with an input trk size of: "<<startingtracks2.size()<<std::endl;

  for(unsigned int j=0;j<startingtracks2.size();j++)
    {
      // bool t2trkcutPassed=false;
      T1T2Track trk=startingtracks2.at(j);
      bx_=trk.X0();
      by_=trk.Y0();	
      if(usexytrack==false)
	{
	  bx_ = trk.GetTx() * cos( trk.Phi() ); 
	  by_ = trk.GetTx() * sin( trk.Phi() );
	}
      if(((TrkInQuarter(trk))&&(UseAnyQuarter==false))||(UseAnyQuarter==true))
	if(fabs(bx_)<Maxbxtrk)
	  if(fabs(by_)<Maxbytrk)
	    if(ChiCutCond(trk, usexytrack,0.01, 0.01)) //added 13 09 2010
	      {
		//if(fabs(trk.Z_at_Rmin())<7000.)    
		startingtracks.push_back(trk);
		//   std::cout<<"From vtx1 "<<trk.Phi()*180.0/3.14159<<std::endl;
	      }
      
    }
  //std::cout<<"2-Vertex function will work with an input trk size of: "<<startingtracks.size()<<std::endl;

  for(unsigned int j=0; j<startingtracks.size();j++)  
    {
      if(allTracksclusters.size()==0)
	{
	  //	  std::cout<<"From vtx2 "<<startingtracks.at(j).Phi()*180.0/3.14159<<std::endl;
	  //  if(fabs(startingtracks.at(j).Z_at_Rmin())<7000)
	  // if(((startingtracks.at(j).Z_at_Rmin()*startingtracks.at(j).GetHitT2(0).GetHitZ()>0)&&(fabs(startingtracks.at(j).Z_at_Rmin())<5000.))||(startingtracks.at(j).Z_at_Rmin()*thetrk.GetHitT2(0).GetHitZ()<0))   //Mod on 22 june 010 

	  double firsthitZ=startingtracks.at(j).GetHitT2(0).GetHitZ();
	  if(startingtracks.at(j).GetHitT2(0).GetHitClass()==9)
	    firsthitZ=startingtracks.at(j).GetHitT2(1).GetHitZ();


	  if((((startingtracks.at(j).Z_at_Rmin()*firsthitZ)>0)&&(fabs(startingtracks.at(j).Z_at_Rmin())<5000.))||((startingtracks.at(j).Z_at_Rmin()*firsthitZ)<0))
	    {
	    if(ThetaAsPrimary(startingtracks.at(j),usexytrack))
	      {
		std::vector<T1T2Track> onetrk;
		onetrk.push_back(startingtracks.at(j));
		allTracksclusters.push_back(onetrk);
		
	      }
	    }
	}
      else
	{
	  bool CanTrkStayInClu=false;
	  unsigned int closercluster=10001; //you will never have 10001 clusters!
	  long double mincludistance=100000.;
	  long double minZc= 100000.;

	  for(unsigned int m=0;m<allTracksclusters.size();m++)
	    {
	      std::vector<T1T2Track> oneexistingcluster=allTracksclusters.at(m);
	      
	  
	      for(unsigned int u=0;u<oneexistingcluster.size();u++)
		{
		  actualtrksep=TwoTracksSeparation(oneexistingcluster.at(u),startingtracks.at(j),usexytrack);
		  // std::cout<<"Trk-Trk distance = "<<actualtrksep.distance<<std::endl;
		  if(actualtrksep.distance<60)//6cm= the beampipe
		    {
		      if(fabs(actualtrksep.midPoint.z)<5000)
			{
			  if(actualtrksep.distance<mincludistance)
			    if(actualtrksep.midPoint.z<minZc)
			      {
				CanTrkStayInClu=true;
				minZc=actualtrksep.midPoint.z;
				mincludistance=actualtrksep.distance;
				closercluster=m;
			      }
			}
		    }
		 
		}
	      

	    }

	  
	  if(CanTrkStayInClu==false) 
	    {
	      //create a new cluster
	      //if(fabs(startingtracks.at(j).Z_at_Rmin())<7000)
	       double firsthitZ=startingtracks.at(j).GetHitT2(0).GetHitZ();
	       if(startingtracks.at(j).GetHitT2(0).GetHitClass()==9)
		 firsthitZ=startingtracks.at(j).GetHitT2(1).GetHitZ();


	      if((((startingtracks.at(j).Z_at_Rmin()*firsthitZ)>0)&&(fabs(startingtracks.at(j).Z_at_Rmin())<5000.))||((startingtracks.at(j).Z_at_Rmin()*firsthitZ)<0))
		if(ThetaAsPrimary(startingtracks.at(j),usexytrack))
		{
		  std::vector<T1T2Track> onetrk;
		  onetrk.push_back(startingtracks.at(j));
		  allTracksclusters.push_back(onetrk);
		}
	    }
	  else
	    {
	      //add  track to an existing cluster
	      if(closercluster >= allTracksclusters.size())
		std::cout<<"Error in PrimaryVertexFinder: FindPrimaryVtx does not work properly!!";
	      else
		if(ThetaAsPrimary(startingtracks.at(j),usexytrack))
		  (allTracksclusters.at(closercluster)).push_back(startingtracks.at(j));	
	    }
	}
    } //end for on startingtracks
  

  // std::cout<<"Here1"<<std::endl;
  // std::cout<<"allTracksclusters.size(): "<<allTracksclusters.size()<<std::endl;
  //Now compute a Vtx position for each cluster 
  std::vector<T2SelectionCutUtils::point3d> allcandidatevtxPos;
  std::vector<T2SelectionCutUtils::point3d> allcandidatevtxPosXY;
  std::vector<T2SelectionCutUtils::point3d> allcandidatevtxPosRZ;
  std::vector<std::vector<T1T2Track> > Finalclusters;

  for(unsigned int f=0;f<allTracksclusters.size();f++)
    {
      std::vector<T1T2Track> onecluster=allTracksclusters.at(f);
      //compute the vertex position as the average of the mid-point distances between trks in cluster 
      if(onecluster.size()>=2)
	{
	  //std::cout<<"onecluster.size()"<<onecluster.size()<<std::endl;
	  T2SelectionCutUtils::point3d theaverage;
	  T2SelectionCutUtils::point3d theaverageRZ;
	  T2SelectionCutUtils::point3d theaverageXY;
	  theaverage.x=0.;
	  theaverage.y=0.;
	  theaverage.z=0.;
	  theaverageRZ=theaverage;
	  theaverageXY=theaverage;
	  long double howmanyadded=0.;
	  for(unsigned int s=0;s<onecluster.size()-1;s++)
	    for(unsigned int ss=s+1;ss<onecluster.size();ss++)
	      {
		//	std::cout<<"Trk-Trk z sep= "<<actualtrksep.midPoint.z<<std::endl;
		actualtrksep=TwoTracksSeparation(onecluster.at(s),onecluster.at(ss),usexytrack);
		theaverage.x=theaverage.x + actualtrksep.midPoint.x;
		theaverage.y=theaverage.y + actualtrksep.midPoint.y;
		theaverage.z=theaverage.z + actualtrksep.midPoint.z;
		std::vector<T2Hit> hitv1=HitsFromTrk(onecluster.at(s));
		std::vector<T2Hit> hitv2=HitsFromTrk(onecluster.at(ss));

		T1T2Track trk1;
		T1T2Track trk2; 
		 
		trk1=TrackFromHits(true,hitv1);
		trk2=TrackFromHits(true,hitv2);  
		actualtrksep=TwoTracksSeparation(trk1,trk2,true);
		theaverageXY.x=theaverageXY.x + actualtrksep.midPoint.x;
		theaverageXY.y=theaverageXY.y + actualtrksep.midPoint.y;
		theaverageXY.z=theaverageXY.z + actualtrksep.midPoint.z;

		trk1=TrackFromHits(false,hitv1);
		trk2=TrackFromHits(false,hitv2);  
		actualtrksep=TwoTracksSeparation(trk1,trk2,false);
		theaverageRZ.x=theaverageRZ.x + actualtrksep.midPoint.x;
		theaverageRZ.y=theaverageRZ.y + actualtrksep.midPoint.y;
		theaverageRZ.z=theaverageRZ.z + actualtrksep.midPoint.z;

		howmanyadded=howmanyadded+1.0;
	      }

	  theaverage.x=theaverage.x / howmanyadded;
	  theaverage.y=theaverage.y / howmanyadded;
	  theaverage.z=theaverage.z / howmanyadded;

	  theaverageXY.x=theaverageXY.x / howmanyadded;
	  theaverageXY.y=theaverageXY.y / howmanyadded;
	  theaverageXY.z=theaverageXY.z / howmanyadded;
	  
	  theaverageRZ.x=theaverageRZ.x / howmanyadded;
	  theaverageRZ.y=theaverageRZ.y / howmanyadded;
	  theaverageRZ.z=theaverageRZ.z / howmanyadded;

	  Finalclusters.push_back(onecluster);

	  allcandidatevtxPos.push_back(theaverage);
	  allcandidatevtxPosRZ.push_back(theaverageRZ);
	  allcandidatevtxPosXY.push_back(theaverageXY);
	}
    }
  // std::cout<<"Here2"<<std::endl;
  //Now select the Cluster you prefear: I take the one with min Z.
  double minz=100000.;
  unsigned int indexofselcluster=0;
  for(unsigned int f=0;f<allcandidatevtxPos.size();f++)
    {
      if(fabs(allcandidatevtxPos.at(f).z)<minz)
	{
	  minz=allcandidatevtxPos.at(f).z;
	  indexofselcluster=f;
	}
    }
 

  if(Finalclusters.size()>0)
    {
      storeVtx.VtxTrks= Finalclusters.at(indexofselcluster);
      storeVtx.VtxPos=allcandidatevtxPos.at(indexofselcluster);
      storeVtx.VtxPosXY=allcandidatevtxPosXY.at(indexofselcluster);
      storeVtx.VtxPosRZ=allcandidatevtxPosRZ.at(indexofselcluster);
      //   std::cout<<"Here3 minz= "<<minz<<std::endl;
    }
  
}




//Warning: this is a copy of the function above. Keep it updated.
std::vector<T2SelectionCutUtils::point3d> T2SelectionCutUtils::GivePrimaryVtx(std::vector<T1T2Track> startingtracks2, bool usexytrack, bool UseAnyQuarter){

  std::vector<T2SelectionCutUtils::point3d> vtxPosToRet;
  
  
  T2SelectionCutUtils::TracksSeparation actualtrksep;
  std::vector<std::vector<T1T2Track> > allTracksclusters;
  std::vector<T1T2Track> startingtracks;
  T2GeometryUtil conv;
  
 
  double bx_,by_;

  for(unsigned int j=0;j<startingtracks2.size();j++)
    {
      // bool t2trkcutPassed=false;
      T1T2Track trk=startingtracks2.at(j);
      bx_=trk.X0();
      by_=trk.Y0();	
      if(usexytrack==false)
	{
	  bx_ = trk.GetTx() * cos( trk.Phi() ); 
	  by_ = trk.GetTx() * sin( trk.Phi() );
	}
     if(((TrkInQuarter(trk))&&(UseAnyQuarter==false))||(UseAnyQuarter==true))
	if(fabs(bx_)<60)
	  if(fabs(by_)<60)
	    if((((trk.Z_at_Rmin()*trk.GetHitT2(0).GetHitZ())>0)&&(fabs(trk.Z_at_Rmin())<5000.))||((trk.Z_at_Rmin()*trk.GetHitT2(0).GetHitZ())<0))   //Mod on 22 june 010   
   	      startingtracks.push_back(trk);
     //(fabs(startingtracks.at(j).Z_at_Rmin())<5000) was before
    }
  // std::cout<<"Vertex trk size: "<<startingtracks.size()<<std::endl;

  for(unsigned int j=0; j<startingtracks.size();j++)  
    {
      if(allTracksclusters.size()==0)
	{
	  double firsthitZ=startingtracks.at(j).GetHitT2(0).GetHitZ();
	  if(startingtracks.at(j).GetHitT2(0).GetHitClass()==9)
	    firsthitZ=startingtracks.at(j).GetHitT2(1).GetHitZ();

	  // if(fabs(startingtracks.at(j).Z_at_Rmin())<5000)//Mod on 22 june 010   
	  if((((startingtracks.at(j).Z_at_Rmin()*firsthitZ)>0)&&(fabs(startingtracks.at(j).Z_at_Rmin())<5000.))||((startingtracks.at(j).Z_at_Rmin()*firsthitZ)<0))
	    if(ThetaAsPrimary(startingtracks.at(j),usexytrack))
	      {
		std::vector<T1T2Track> onetrk;
		onetrk.push_back(startingtracks.at(j));
		allTracksclusters.push_back(onetrk);	      
	      }
	 
	}
      else
	{
	  bool CanTrkStayInClu=false;
	  unsigned int closercluster=10001; //you will never have 10001 clusters!
	  long double mincludistance=100000.;
	  long double minZc= 100000.;

	  for(unsigned int m=0;m<allTracksclusters.size();m++)
	    {
	      std::vector<T1T2Track> oneexistingcluster=allTracksclusters.at(m);
	      
	  
	      for(unsigned int u=0;u<oneexistingcluster.size();u++)
		{
		  actualtrksep=TwoTracksSeparation(oneexistingcluster.at(u),startingtracks.at(j),usexytrack);
		  
		  if(actualtrksep.distance<60)//6cm= the beampipe
		    {
		      if(fabs(actualtrksep.midPoint.z)<5000)
			{
			  if(actualtrksep.distance<mincludistance)
			    if(actualtrksep.midPoint.z<minZc)
			      {
				CanTrkStayInClu=true;
				minZc=actualtrksep.midPoint.z;
				mincludistance=actualtrksep.distance;
				closercluster=m;
			      }
			}
		    }
		 
		}
	      

	    }

	  
	  if(CanTrkStayInClu==false) 
	    {
	      //create a new cluster
	      //if(fabs(startingtracks.at(j).Z_at_Rmin())<5000)

	       double firsthitZ=startingtracks.at(j).GetHitT2(0).GetHitZ();
	       if(startingtracks.at(j).GetHitT2(0).GetHitClass()==9)
		 firsthitZ=startingtracks.at(j).GetHitT2(1).GetHitZ();

	      

	      if((((startingtracks.at(j).Z_at_Rmin()*firsthitZ)>0)&&(fabs(startingtracks.at(j).Z_at_Rmin())<5000.))||((startingtracks.at(j).Z_at_Rmin()*firsthitZ)<0))
	      if(ThetaAsPrimary(startingtracks.at(j),usexytrack))
		{
		  std::vector<T1T2Track> onetrk;
		  onetrk.push_back(startingtracks.at(j));
		  allTracksclusters.push_back(onetrk);
		}
	    }
	  else
	    {
	      //add  track to an existing cluster
	      if(closercluster >= allTracksclusters.size())
		std::cout<<"Error in GivePrimaryVtx: FindPrimaryVtx does not work properly!!";
	      else
		if(ThetaAsPrimary(startingtracks.at(j),usexytrack))
		  (allTracksclusters.at(closercluster)).push_back(startingtracks.at(j));	
	    }
	}
    } //end for on startingtracks
  

 
  //Now compute a Vtx position for each cluster 
  std::vector<T2SelectionCutUtils::point3d> allcandidatevtxPos;
  std::vector<T2SelectionCutUtils::point3d> allcandidatevtxPosXY;
  std::vector<T2SelectionCutUtils::point3d> allcandidatevtxPosRZ;
  std::vector<std::vector<T1T2Track> > Finalclusters;

  for(unsigned int f=0;f<allTracksclusters.size();f++)
    {
      std::vector<T1T2Track> onecluster=allTracksclusters.at(f);
      //compute the vertex position as the average of the mid-point distances between trks in cluster 
      if(onecluster.size()>=2)
	{	  
	  T2SelectionCutUtils::point3d theaverage;
	  T2SelectionCutUtils::point3d theaverageRZ;
	  T2SelectionCutUtils::point3d theaverageXY;
	  theaverage.x=0.;
	  theaverage.y=0.;
	  theaverage.z=0.;
	  theaverageRZ=theaverage;
	  theaverageXY=theaverage;
	  long double howmanyadded=0.;
	  for(unsigned int s=0;s<onecluster.size()-1;s++)
	    for(unsigned int ss=s+1;ss<onecluster.size();ss++)
	      {		
		actualtrksep=TwoTracksSeparation(onecluster.at(s),onecluster.at(ss),usexytrack);
		theaverage.x=theaverage.x + actualtrksep.midPoint.x;
		theaverage.y=theaverage.y + actualtrksep.midPoint.y;
		theaverage.z=theaverage.z + actualtrksep.midPoint.z;
		std::vector<T2Hit> hitv1=HitsFromTrk(onecluster.at(s));
		std::vector<T2Hit> hitv2=HitsFromTrk(onecluster.at(ss));

		T1T2Track trk1;
		T1T2Track trk2; 
		 
		trk1=TrackFromHits(true,hitv1);
		trk2=TrackFromHits(true,hitv2);  
		actualtrksep=TwoTracksSeparation(trk1,trk2,true);
		theaverageXY.x=theaverageXY.x + actualtrksep.midPoint.x;
		theaverageXY.y=theaverageXY.y + actualtrksep.midPoint.y;
		theaverageXY.z=theaverageXY.z + actualtrksep.midPoint.z;

		trk1=TrackFromHits(false,hitv1);
		trk2=TrackFromHits(false,hitv2);  
		actualtrksep=TwoTracksSeparation(trk1,trk2,false);
		theaverageRZ.x=theaverageRZ.x + actualtrksep.midPoint.x;
		theaverageRZ.y=theaverageRZ.y + actualtrksep.midPoint.y;
		theaverageRZ.z=theaverageRZ.z + actualtrksep.midPoint.z;

		howmanyadded=howmanyadded+1.0;
	      }

	  theaverage.x=theaverage.x / howmanyadded;
	  theaverage.y=theaverage.y / howmanyadded;
	  theaverage.z=theaverage.z / howmanyadded;

	  theaverageXY.x=theaverageXY.x / howmanyadded;
	  theaverageXY.y=theaverageXY.y / howmanyadded;
	  theaverageXY.z=theaverageXY.z / howmanyadded;
	  
	  theaverageRZ.x=theaverageRZ.x / howmanyadded;
	  theaverageRZ.y=theaverageRZ.y / howmanyadded;
	  theaverageRZ.z=theaverageRZ.z / howmanyadded;

	  Finalclusters.push_back(onecluster);

	  allcandidatevtxPos.push_back(theaverage);	  //BUG: this line added 12-09-2010
	  allcandidatevtxPosRZ.push_back(theaverageRZ);
	  allcandidatevtxPosXY.push_back(theaverageXY);
	 
	  
	  
	}
     
      
    }


double minz=100000.;
  unsigned int indexofselcluster=0;
  for(unsigned int f=0;f<allcandidatevtxPos.size();f++)
    {
      if(fabs(allcandidatevtxPos.at(f).z)<minz)
	{
	  minz=allcandidatevtxPos.at(f).z;
	  indexofselcluster=f;
	}
    }
 

  if(Finalclusters.size()>0)
    {
      vtxPosToRet.push_back(allcandidatevtxPosXY.at(indexofselcluster));
      vtxPosToRet.push_back(allcandidatevtxPosRZ.at(indexofselcluster));      
    }

  return vtxPosToRet;

  
}




std::vector<T2SelectionCutUtils::CandidateVtxWithTrks> T2SelectionCutUtils::FindAllVtxs(std::vector<T1T2Track> startingtracks2, bool usexytrack, std::vector<int> Quartertouse, double ClusteringDistance_mm){

  std::vector<T2SelectionCutUtils::point3d> vtxPosToRet;
  std::vector<T2SelectionCutUtils::CandidateVtxWithTrks> vtxSToRet;
  T2SelectionCutUtils::TracksSeparation actualtrksep;
  std::vector<std::vector<T1T2Track> > allTracksclusters;
  std::vector<T1T2Track> startingtracks;
  T2GeometryUtil conv;
  


  for(unsigned int j=0;j<startingtracks2.size();j++)
    {
      // bool t2trkcutPassed=false;
      T1T2Track trk=startingtracks2.at(j);
      trk.X0();
      trk.Y0();
      if(usexytrack==false)
	{
	  trk.GetTx();
	  cos( trk.Phi() );
	  trk.GetTx();
	  sin( trk.Phi() );
	}

      //if(((TrkInQuarter(trk))&&(UseAnyQuarter==false))||(UseAnyQuarter==true))
      if(TrkAlsoInQuarter(trk,Quartertouse))
	 startingtracks.push_back(trk);
      
    }
  
  //std::cout<<"Vertex trk size: "<<startingtracks.size()<<std::endl;

  for(unsigned int j=0; j<startingtracks.size();j++)  
    {
      if(allTracksclusters.size()==0)
	{
	  startingtracks.at(j).GetHitT2(0).GetHitZ();
	  if(startingtracks.at(j).GetHitT2(0).GetHitClass()==9)
	    startingtracks.at(j).GetHitT2(1).GetHitZ();
	  
	  // if(fabs(startingtracks.at(j).Z_at_Rmin())<5000)//Mod on 22 june 010   
	  
	  std::vector<T1T2Track> onetrk;
	  onetrk.push_back(startingtracks.at(j));
	  allTracksclusters.push_back(onetrk);	      
	      
	  
	}
      else
	{
	  bool CanTrkStayInClu=false;
	  unsigned int closercluster=10001; //you will never have 10001 clusters!
	  long double mincludistance=100000.;
	  long double minZc= 100000.;


	  


	  for(unsigned int m=0;m<allTracksclusters.size();m++)
	    {
	      std::vector<T1T2Track> oneexistingcluster=allTracksclusters.at(m);
	      
	  
	      for(unsigned int u=0;u<oneexistingcluster.size();u++)
		{
		  actualtrksep=TwoTracksSeparation(oneexistingcluster.at(u),startingtracks.at(j),usexytrack);
		  
		  
		  if(actualtrksep.distance<mincludistance)
		    if(actualtrksep.midPoint.z<minZc)
		      {
			if(actualtrksep.distance<ClusteringDistance_mm)
			  CanTrkStayInClu=true;

			minZc=actualtrksep.midPoint.z;
			mincludistance=actualtrksep.distance;
			closercluster=m;
		      }
		
		  
		}
	      

	    }
	 

	  
	  if(CanTrkStayInClu==false) 
	    {
	      //create a new cluster
	      //if(fabs(startingtracks.at(j).Z_at_Rmin())<5000)

	      
	      std::vector<T1T2Track> onetrk;
	      onetrk.push_back(startingtracks.at(j));
	      allTracksclusters.push_back(onetrk);
		
	    }
	  else
	    {
	      //add  track to an existing cluster
	      if(closercluster >= allTracksclusters.size())
		std::cout<<"Error in GivePrimaryVtx: FindPrimaryVtx does not work properly!!";
	      else		
		(allTracksclusters.at(closercluster)).push_back(startingtracks.at(j));	
	    }
	}
    } //end for on startingtracks
  

 
  //Now compute a Vtx position for each cluster 
  // std::vector<T2SelectionCutUtils::point3d> allcandidatevtxPos;
  //std::vector<T2SelectionCutUtils::point3d> allcandidatevtxPosXY;
  //std::vector<T2SelectionCutUtils::point3d> allcandidatevtxPosRZ;
  //std::vector<std::vector<T1T2Track> > Finalclusters;

  for(unsigned int f=0;f<allTracksclusters.size();f++)
    {
      std::vector<T1T2Track> onecluster=allTracksclusters.at(f);
      //compute the vertex position as the average of the mid-point distances between trks in cluster 
      if(onecluster.size()>=2)
	{	  
	  T2SelectionCutUtils::point3d theaverage;
	  T2SelectionCutUtils::point3d theaverageRZ;
	  T2SelectionCutUtils::point3d theaverageXY;
	  theaverage.x=0.;
	  theaverage.y=0.;
	  theaverage.z=0.;
	  theaverageRZ=theaverage;
	  theaverageXY=theaverage;
	  long double howmanyadded=0.;
	  for(unsigned int s=0;s<onecluster.size()-1;s++)
	    for(unsigned int ss=s+1;ss<onecluster.size();ss++)
	      {		
		actualtrksep=TwoTracksSeparation(onecluster.at(s),onecluster.at(ss),usexytrack);
		theaverage.x=theaverage.x + actualtrksep.midPoint.x;
		theaverage.y=theaverage.y + actualtrksep.midPoint.y;
		theaverage.z=theaverage.z + actualtrksep.midPoint.z;
		std::vector<T2Hit> hitv1=HitsFromTrk(onecluster.at(s));
		std::vector<T2Hit> hitv2=HitsFromTrk(onecluster.at(ss));

		T1T2Track trk1;
		T1T2Track trk2; 
		 
		trk1=TrackFromHits(true,hitv1);
		trk2=TrackFromHits(true,hitv2);  
		actualtrksep=TwoTracksSeparation(trk1,trk2,true);
		theaverageXY.x=theaverageXY.x + actualtrksep.midPoint.x;
		theaverageXY.y=theaverageXY.y + actualtrksep.midPoint.y;
		theaverageXY.z=theaverageXY.z + actualtrksep.midPoint.z;

		trk1=TrackFromHits(false,hitv1);
		trk2=TrackFromHits(false,hitv2);  
		actualtrksep=TwoTracksSeparation(trk1,trk2,false);
		theaverageRZ.x=theaverageRZ.x + actualtrksep.midPoint.x;
		theaverageRZ.y=theaverageRZ.y + actualtrksep.midPoint.y;
		theaverageRZ.z=theaverageRZ.z + actualtrksep.midPoint.z;

		howmanyadded=howmanyadded+1.0;
	      }

	  theaverage.x=theaverage.x / howmanyadded;
	  theaverage.y=theaverage.y / howmanyadded;
	  theaverage.z=theaverage.z / howmanyadded;

	  theaverageXY.x=theaverageXY.x / howmanyadded;
	  theaverageXY.y=theaverageXY.y / howmanyadded;
	  theaverageXY.z=theaverageXY.z / howmanyadded;
	  
	  theaverageRZ.x=theaverageRZ.x / howmanyadded;
	  theaverageRZ.y=theaverageRZ.y / howmanyadded;
	  theaverageRZ.z=theaverageRZ.z / howmanyadded;

	  //Finalclusters.push_back(onecluster);
	  
	  T2SelectionCutUtils::CandidateVtxWithTrks thevtx;
	  thevtx.VtxPos= theaverageXY;
	  thevtx.VtxPosXY= theaverageXY;
	  thevtx.VtxPosRZ=  theaverageRZ;
	  thevtx.VtxTrks = onecluster;
	    
	  vtxSToRet.push_back(thevtx);
	  	  
	  
	}
    }
 
  //vtxSToRet ordering according to |Z| position
  
  std::vector<unsigned int> vtxpositionindex;
  
  double actualZmin=0.;
  unsigned int actualZminIndex=0.;
  for(unsigned int g=0;g<vtxSToRet.size();g++)
    {
      actualZmin=vtxSToRet.at(g).VtxPosXY.z;
      actualZminIndex=g;

      for(unsigned int gg=g;gg<vtxSToRet.size();gg++)
	{
	  if(fabs(vtxSToRet.at(gg).VtxPosXY.z)<actualZmin)
	    {
	      actualZmin=vtxSToRet.at(gg).VtxPosXY.z;
	      actualZminIndex=gg;
	    }
	}
      
      vtxpositionindex.push_back(actualZminIndex);
      
    }
  std::vector<T2SelectionCutUtils::CandidateVtxWithTrks> vtxSToRetOrdered;
  for(unsigned int h=0;h<vtxpositionindex.size();h++)
    {
      vtxSToRetOrdered.push_back(vtxSToRet.at(vtxpositionindex.at(h)));
    }

  return vtxSToRetOrdered;
  //return vtxSToRet;
  
}













//return x-y-z
std::vector<double> vtxPosToRet(std::vector<T2SelectionCutUtils::point3d> allvtxs) 
{
  std::vector<double> XYZ;
  double thex=0.;
  double they=0.;
  double thez=0.;
  
  for(unsigned int i=0;i<allvtxs.size();i++)
    {
      thex+=allvtxs.at(i).x;
      they+=allvtxs.at(i).y;
      thez+=allvtxs.at(i).z;
    }
  
  thex=thex/allvtxs.size();
  they=they/allvtxs.size();
  thez=thez/allvtxs.size();

  XYZ.push_back(thex);
  XYZ.push_back(they);
  XYZ.push_back(thez);

  return XYZ;
}


double T2SelectionCutUtils::FindAverageTrackZ(std::vector<std::vector<T2Hit> >* theTracks, char XorY){

  double Zavg=0.;
  unsigned int i=0;
  
  for(unsigned int f=0;f<theTracks->size();f++)
    {
      std::vector<double> corrFit=MyLinearfitXYSimple(theTracks->at(f));
      if(XorY=='X')
	{
	  Zavg+=(corrFit[1]/corrFit[0]);
	  i++;
	}
      else
	{
	  Zavg+=(corrFit[3]/corrFit[2]);
	  i++;
	}
    }

  Zavg=(-1.0)*Zavg/((double)i);
  

  // TVectorD vect(4);
  //vect[0] = corrFit[1];//b_xz;
  //vect[1] = corrFit[3];//b_yz;
  //vect[2] = corrFit[0];//a_xz;
  //vect[3] = corrFit[2];//a_yz;

  return Zavg;
}




std::vector<T2Hit> T2SelectionCutUtils::HitsFromTrk(T1T2Track trk){
  
  std::vector<T2Hit> hitv;
  unsigned int numhit=trk.GetHitEntries();
  for(unsigned int i=0;i<numhit;i++)
    hitv.push_back(trk.GetHitT2(i));
  
  return hitv;
}




std::vector<T2Hit> T2SelectionCutUtils::HitsFromTrkInOneQuarter(T1T2Track trk, unsigned int selectedquarter){
  

  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo;

  std::vector<T2Hit> hitv;
  unsigned int numhit=trk.GetHitEntries();
  unsigned int symbh0;
  for(unsigned int i=0;i<numhit;i++)
    {
      if(trk.GetHitT2(i).GetHitClass()!=9)
	{
	  planeinfo=conv.GetT2Info(trk.GetHitT2(i).GetHitDetRawId());
	  symbh0=planeinfo.arm*2+planeinfo.ht; 
	   if(symbh0==selectedquarter)
	     hitv.push_back(trk.GetHitT2(i));
	}
      else
	hitv.push_back(trk.GetHitT2(i));
    }
  return hitv;
}



std::vector<T2Hit> T2SelectionCutUtils::HitsFromTrkInOneQuarter( std::vector<T2Hit> hitv, unsigned int selectedquarter){ 

  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinfo;

  unsigned int numhit=hitv.size();
  unsigned int symbh0;
  for(unsigned int i=0;i<numhit;i++)
    {
      if(hitv.at(i).GetHitClass()!=9)
	{
	  planeinfo=conv.GetT2Info(hitv.at(i).GetHitDetRawId());
	  symbh0=planeinfo.arm*2+planeinfo.ht; 
	   if(symbh0==selectedquarter)
	     hitv.push_back(hitv.at(i));
	}
      else
	hitv.push_back(hitv.at(i));
    }
  return hitv;
}









// -1 no wrost hit found
// >= 0 worst hit index
int T2SelectionCutUtils::findWorstPoint(const T1T2Track &trk)//(T1T2Track* trk)
{
double  DropWorstRecHitChisquareRThreshold = 0.01;

  std::vector<double> RecoHitRErrorPerStripcount; //in the T2TrackProducer taken as input
  RecoHitRErrorPerStripcount.push_back(0.09);
  RecoHitRErrorPerStripcount.push_back(0.134);
  RecoHitRErrorPerStripcount.push_back(0.109);
  RecoHitRErrorPerStripcount.push_back(0.109);
  RecoHitRErrorPerStripcount.push_back(0.118);

  std::vector<T2Hit> hitvec;  
  std::vector<float> r;
  std::vector<float> z;
  std::vector<float> er;
  std::vector<float> ez;
  //std::cout<<"A"<<std::endl;
  unsigned int sizeHitv=trk.GetHitEntries();

  //Push the staff into vectors
  //IMPORTANT hitvecraw should be ordered in the producer!!
  for (unsigned int jj = 0;jj<sizeHitv;jj++){
    T2Hit hit = trk.GetHitT2(jj);

    r.push_back(hit.GetHitR());
    z.push_back(hit.GetHitZ());
    //std::cout<<"A1"<<std::endl;
    if (hit.GetHitClass() == 1){
      //Hit contains strips
      unsigned int stripCount = hit.GetHitNumStrip();
      if (stripCount < RecoHitRErrorPerStripcount.size()){	
	er.push_back(RecoHitRErrorPerStripcount[stripCount]);
      } else {
	er.push_back(RecoHitRErrorPerStripcount[RecoHitRErrorPerStripcount.size()-1]);
      }
    } else {
      //Only pad hit or Vtx
      er.push_back(hit.GetHitDR());
    }
    //std::cout<<"A4"<<std::endl;
    //ez.push_back(hitvec[jj].GetHitDZ());
    ez.push_back(hit.GetHitDZ());
  }
  // std::cout<<"B"<<std::endl;
  bool biggestContributionChiRSet = false;
  double biggestContributionChiR = 0;
  unsigned int worstHitPointIndexR = 0;
  
  float Sr=0;
  float Srz=0;
  float Szz_r=0;
  float Sz_r=0; 
  float S0_r=0; 

  //Calculate different sums
  for(unsigned int jj =0; jj<sizeHitv; jj++){      
    Srz += r[jj]*z[jj]/er[jj]/er[jj];
    Szz_r += z[jj]*z[jj]/er[jj]/er[jj];
    Sz_r += z[jj]/er[jj]/er[jj];
    Sr += r[jj]/er[jj]/er[jj];
    S0_r += 1.0/er[jj]/er[jj];
  }

  //Calculate the coefficients
  double a_rz = (Srz*S0_r - Sz_r*Sr) / (Szz_r*S0_r - Sz_r*Sz_r);   // angular coefficient
  double b_rz = (Sr*Szz_r - Sz_r*Srz) / (Szz_r*S0_r - Sz_r*Sz_r);  // intercept   R=(a_rz)Z + b_rz
  //double e_a_rz = sqrt( S0_r / (S0_r*Szz_r - Sz_r*Sz_r) );         
  //double e_b_rz = sqrt( Szz_r / (S0_r*Szz_r - Sz_r*Sz_r) );
  //std::cout<<"C"<<std::endl;
  //Calculate the chisquare and keep track which gives the worst
  double chi2r = 0.0;
  for(unsigned int jjj =0; jjj<sizeHitv; jjj++){    
    double chi2radd = (a_rz*z[jjj]+b_rz - r[jjj])*(a_rz*z[jjj]+b_rz - r[jjj])/er[jjj]/er[jjj];
    chi2r += chi2radd;
    
    if(biggestContributionChiRSet == false || biggestContributionChiR < chi2radd){     
      biggestContributionChiRSet = true;
      biggestContributionChiR = chi2radd;
      worstHitPointIndexR = jjj;
    }
  }
  
  //No this chi2r is calculated with all the hit. Suppose yuo found the worst hit and yor chi2 before was bad.
  //Using this condition you don't give the possibility to fit the track with the worst hit removed!!
  
  if (TMath::Prob(chi2r,sizeHitv-2) >= DropWorstRecHitChisquareRThreshold || biggestContributionChiRSet == false) {
    //Everything is OK. No worst point found.
    return -1;
  }
  /*
  if(verbosity)
    std::cout<<"Worst index in trk: "<<(int) worstHitPointIndexR<<std::endl;
  */
  //if(trk.GetHitT2(worstHitPointIndexR).GetHitClass==9)
    

  return (int) worstHitPointIndexR;
  
}





T1T2Track T2SelectionCutUtils::fitTrackspecialXY(std::vector<T2Hit> hitvec, bool  CheckWorstHit)
{
   
  
   unsigned int sizeHitv=hitvec.size();   

   int hemisphere=0; 

   if(sizeHitv>=2)
     {
       hemisphere = (int)(hitvec[1].GetHitZ()/fabs(hitvec[1].GetHitZ()));
     }
   else
     {
        hemisphere = (int)(hitvec[0].GetHitZ()/fabs(hitvec[0].GetHitZ()));
	std::cout<<"ERROR in T2SelectionCutUtils::fitTrackspecialXY: using vtx position to determine the hemisphere "<<std::endl;
     }


double Sx=0.;
double Sxz=0.;
double Szz_x=0.;
double Sz_x=0.; 
double S0_x=0.; 

double Sy=0.;
double Syz=0.;
double Szz_y=0.;
double Sz_y=0.; 
double S0_y=0.; 

std::vector<double> x;   
std::vector<double> z;
std::vector<double> y;
std::vector<double> ex;
std::vector<double> ez;
std::vector<double> ey;
 double a_xz, b_xz, a_yz, b_yz;

  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      //   std::cout<<hitvec[jj].GetHitX()<<" "<<hitvec[jj].GetHitPhi()<<" "<<hitvec[jj].GetHitZ()<<std::endl;
      //r->x   phi->y
      
      x.push_back(hitvec[jj].GetHitX());
      y.push_back(hitvec[jj].GetHitY());
      z.push_back(hitvec[jj].GetHitZ());    
 
      double phirad=hitvec[jj].GetHitPhi()*3.14159/180.0;
      double sigmay;
      double sigmax;
      //er.push_back(0.1);
      if(hitvec[jj].GetHitClass()==1)
	{
	  sigmax=cos(phirad)*cos(phirad)*0.12*0.12+sin(phirad)*sin(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
	  ex.push_back(sqrt(sigmax));
      
	  //ephi.push_back(0.1);
	  sigmay=sin(phirad)*sin(phirad)*0.12*0.12+cos(phirad)*cos(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
	  ey.push_back(sqrt(sigmay));
	}
      else   //To be tuned properly !!!!!!!!!!!!!
	{
	  double dr=hitvec[jj].GetHitDR();

	  if(hitvec[jj].GetHitClass()!=9)//9=THE VTX
	    {
	      sigmax=cos(phirad)*cos(phirad)*dr*dr+sin(phirad)*sin(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
	      sigmay=sin(phirad)*sin(phirad)*dr*dr+cos(phirad)*cos(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
	    }
	  else
	    {	      
	       sigmax=hitvec[jj].GetHitDX()*hitvec[jj].GetHitDX();
	       sigmay=hitvec[jj].GetHitDY()*hitvec[jj].GetHitDY();		
	    }
	 
	  ex.push_back(sqrt(sigmax));      
	  //ephi.push_back(0.1);
	  
	  ey.push_back(sqrt(sigmay));

	}
      if((sigmax==0)||(sigmay==0))
	std::cout<<"ERROR in T2TrackProducer2::fitTrackspecialXY: sigmaX=  "<<sigmax<<" sigmaY= "<<sigmay<<std::endl;
      //ez.push_back(hitvec[jj].GetHitDZ());
            
      Sxz += x[jj]*z[jj]/ex[jj]/ex[jj];
      Szz_x += z[jj]*z[jj]/ex[jj]/ex[jj];
      Sz_x += z[jj]/ex[jj]/ex[jj];
      Sx += x[jj]/ex[jj]/ex[jj];
      S0_x += 1.0/ex[jj]/ex[jj];



      Syz += y[jj]*z[jj]/ey[jj]/ey[jj];
      Szz_y += z[jj]*z[jj]/ey[jj]/ey[jj];
      Sz_y += z[jj]/ey[jj]/ey[jj];
      Sy += y[jj]/ey[jj]/ey[jj];
      S0_y += 1.0/ey[jj]/ey[jj];

  }


  //std::cout<<" Here 2 "<<std::endl;

a_xz = (Sxz*S0_x - Sz_x*Sx) / (Szz_x*S0_x - Sz_x*Sz_x);   // angular coefficient
b_xz = (Sx*Szz_x - Sz_x*Sxz) / (Szz_x*S0_x - Sz_x*Sz_x);  // intercept   X=(a_xz)Z + b_xz

a_yz = (Syz*S0_y - Sz_y*Sy) / (Szz_y*S0_y - Sz_y*Sz_y);   // angular coefficient
b_yz = (Sy*Szz_y - Sz_y*Syz) / (Szz_y*S0_y - Sz_y*Sz_y);  // intercept  Y=(a_yz)Z + b_yz 

double e_a_xz = sqrt( S0_x / (S0_x*Szz_x - Sz_x*Sz_x) );
double e_b_xz = sqrt( Szz_x / (S0_x*Szz_x - Sz_x*Sz_x) );
double e_a_yz = sqrt( S0_y / (S0_y*Szz_y - Sz_y*Sz_y) );
double e_b_yz = sqrt( Szz_y / (S0_y*Szz_y - Sz_y*Sz_y) );

//USE T1 Convenction
TVectorD vect(4);
 for(int oo=0; oo<4; oo++)
   vect[oo]=0.;

 vect[0] = b_xz;
 vect[1] = b_yz;
 vect[2] = a_xz;
 vect[3] = a_yz;
 
//double correlx=(-1.0)*Sz_x*(1.0/(S0_x*Szz_x-(Sz_x*Sz_x)));
//double correly=(-1.0)*Sz_y*(1.0/(S0_y*Szz_y-(Sz_y*Sz_y)));

 

 TMatrixD mat(4,4);
 for(unsigned int oo=0; oo<4; oo++)
   for(unsigned int ooo=0;ooo<4; ooo++)
     mat[oo][ooo]=0.;


 mat[0][0] = e_b_xz*e_b_xz;
    mat[1][1] = e_b_yz*e_b_yz;
    mat[2][2] = e_a_xz*e_a_xz;
    mat[3][3] = e_a_yz*e_a_yz;

    double chi2 = 0;
    double chi2X = 0;
    double chi2Y = 0;

for(unsigned int jjj =0; jjj<sizeHitv; jjj++)
    {
      chi2X += (a_xz*z[jjj]+b_xz - x[jjj])*(a_xz*z[jjj]+b_xz - x[jjj])/ex[jjj]/ex[jjj];
      chi2Y += (a_yz*z[jjj]+b_yz - y[jjj])*(a_yz*z[jjj]+b_yz - y[jjj])/ey[jjj]/ey[jjj];
    }

 chi2 = (chi2X+chi2Y)/2.0; 


 // std::cout<<"chi2X chi2Y "<<chi2X<<" "<<chi2Y<<std::endl;

T1T2Track fittedtrack(vect,mat,chi2,chi2X,chi2Y,hemisphere,2);
//std::cout<<"chi2X after: "<<fittedtrack.ChiSquaredX()<<"  chi2Y after: "<<fittedtrack.ChiSquaredY()<<std::endl;
//sizeHitv=trk.GetHitEntries();
 sizeHitv=hitvec.size();
for(unsigned int jj =0; jj<sizeHitv; jj++)
  {
    fittedtrack.AddHit(hitvec[jj]);
  }
 //fittedtrack.AddHit(trk.GetHitT2(jj));
  // fittedtrack.AddHit(hitvec[jj]);

//  std::cout<<" Here 4 "<<std::endl;

   return fittedtrack; 
}










/*

{
   
unsigned int sizeHitv=hitvec.size();

 int hemisphere = (int)(hitvec[0].GetHitZ()/fabs(hitvec[0].GetHitZ()));

double Sx=0.;
double Sxz=0.;
double Szz_x=0.;
double Sz_x=0.; 
double S0_x=0.; 

double Sy=0.;
double Syz=0.;
double Szz_y=0.;
double Sz_y=0.; 
double S0_y=0.; 

std::vector<double> x;   
std::vector<double> z;
std::vector<double> y;
std::vector<double> ex;
std::vector<double> ez;
std::vector<double> ey;
 double a_xz, b_xz, a_yz, b_yz;

  for(unsigned int jj =0; jj<sizeHitv; jj++)
    {
      //   std::cout<<hitvec[jj].GetHitX()<<" "<<hitvec[jj].GetHitPhi()<<" "<<hitvec[jj].GetHitZ()<<std::endl;
      //r->x   phi->y
      
      x.push_back(hitvec[jj].GetHitX());
      y.push_back(hitvec[jj].GetHitY());
      z.push_back(hitvec[jj].GetHitZ());    
 
      double phirad=hitvec[jj].GetHitPhi()*3.14159/180.0;
      double sigmay;
      double sigmax;
      //er.push_back(0.1);
      if(hitvec[jj].GetHitClass()==1)
	{
	  sigmax=cos(phirad)*cos(phirad)*0.12*0.12+sin(phirad)*sin(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
	  ex.push_back(sqrt(sigmax));
      
	  //ephi.push_back(0.1);
	  sigmay=sin(phirad)*sin(phirad)*0.12*0.12+cos(phirad)*cos(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
	  ey.push_back(sqrt(sigmay));
	}
      else   //To be tuned properly !!!!!!!!!!!!!
	{
	  double dr=hitvec[jj].GetHitDR();
	  sigmax=cos(phirad)*cos(phirad)*dr*dr+sin(phirad)*sin(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
	  ex.push_back(sqrt(sigmax));
      
	  //ephi.push_back(0.1);
	  sigmay=sin(phirad)*sin(phirad)*dr*dr+cos(phirad)*cos(phirad)*0.015*0.015*hitvec[jj].GetHitR()*hitvec[jj].GetHitR();
	  ey.push_back(sqrt(sigmay));

	}
      if((sigmax==0)||(sigmay==0))
	std::cout<<"ERROR in T2TrackProducer2::fitTrackspecialXY: sigmaX=  "<<sigmax<<" sigmaY= "<<sigmay<<std::endl;
      //ez.push_back(hitvec[jj].GetHitDZ());
            
      Sxz += x[jj]*z[jj]/ex[jj]/ex[jj];
      Szz_x += z[jj]*z[jj]/ex[jj]/ex[jj];
      Sz_x += z[jj]/ex[jj]/ex[jj];
      Sx += x[jj]/ex[jj]/ex[jj];
      S0_x += 1.0/ex[jj]/ex[jj];



      Syz += y[jj]*z[jj]/ey[jj]/ey[jj];
      Szz_y += z[jj]*z[jj]/ey[jj]/ey[jj];
      Sz_y += z[jj]/ey[jj]/ey[jj];
      Sy += y[jj]/ey[jj]/ey[jj];
      S0_y += 1.0/ey[jj]/ey[jj];

  }


  //std::cout<<" Here 2 "<<std::endl;

a_xz = (Sxz*S0_x - Sz_x*Sx) / (Szz_x*S0_x - Sz_x*Sz_x);   // angular coefficient
b_xz = (Sx*Szz_x - Sz_x*Sxz) / (Szz_x*S0_x - Sz_x*Sz_x);  // intercept   X=(a_xz)Z + b_xz

a_yz = (Syz*S0_y - Sz_y*Sy) / (Szz_y*S0_y - Sz_y*Sz_y);   // angular coefficient
b_yz = (Sy*Szz_y - Sz_y*Syz) / (Szz_y*S0_y - Sz_y*Sz_y);  // intercept  Y=(a_yz)Z + b_yz 

double e_a_xz = sqrt( S0_x / (S0_x*Szz_x - Sz_x*Sz_x) );
double e_b_xz = sqrt( Szz_x / (S0_x*Szz_x - Sz_x*Sz_x) );
double e_a_yz = sqrt( S0_y / (S0_y*Szz_y - Sz_y*Sz_y) );
double e_b_yz = sqrt( Szz_y / (S0_y*Szz_y - Sz_y*Sz_y) );

//USE T1 Convenction
TVectorD vect(4);
 for(int oo=0; oo<4; oo++)
   vect[oo]=0.;

 vect[0] = b_xz;
 vect[1] = b_yz;
 vect[2] = a_xz;
 vect[3] = a_yz;
 
double correlx=(-1.0)*Sz_x*(1.0/(S0_x*Szz_x-(Sz_x*Sz_x)));
double correly=(-1.0)*Sz_y*(1.0/(S0_y*Szz_y-(Sz_y*Sz_y)));

 

 TMatrixD mat(4,4);
 for(unsigned int oo=0; oo<4; oo++)
   for(unsigned int ooo=0;ooo<4; ooo++)
     mat[oo][ooo]=0.;


 mat[0][0] = e_b_xz*e_b_xz;
    mat[1][1] = e_b_yz*e_b_yz;
    mat[2][2] = e_a_xz*e_a_xz;
    mat[3][3] = e_a_yz*e_a_yz;

    double chi2 = 0;
    double chi2X = 0;
    double chi2Y = 0;

for(unsigned int jjj =0; jjj<sizeHitv; jjj++)
    {
      chi2X += (a_xz*z[jjj]+b_xz - x[jjj])*(a_xz*z[jjj]+b_xz - x[jjj])/ex[jjj]/ex[jjj];
      chi2Y += (a_yz*z[jjj]+b_yz - y[jjj])*(a_yz*z[jjj]+b_yz - y[jjj])/ey[jjj]/ey[jjj];
    }

 chi2 = (chi2X+chi2Y)/2.0;


 //std::cout<<"chi2X chi2Y "<<chi2X<<" "<<chi2Y<<std::endl;

T1T2Track fittedtrack(vect,mat,chi2,chi2X,chi2Y,hemisphere,2);
//std::cout<<"chi2X after: "<<fittedtrack.ChiSquaredX()<<"  chi2Y after: "<<fittedtrack.ChiSquaredY()<<std::endl;

 sizeHitv=hitvec.size();
 //std::cout<<"..."<<sizeHitv<<std::endl;
for(unsigned int jj =0; jj<sizeHitv; jj++)
  {
    fittedtrack.AddHit(hitvec[jj]);
  }
 
//  std::cout<<" Here 4 "<<std::endl;

   return fittedtrack; 
}


*/





/*
Modified by Fabrizio Ferro - INFN Genova for TOTEM
*/

#include <cmath>

#include "DataFormats/T1T2Track/interface/T1T2Track.h"


TMatrixD T1T2Track::TrackPointInterpolationCovariance(double z) const
{
  TMatrixD h(2,4);
  h(0,0)=1;
  h(1,1)=1;
  h(0,2)=z-z0_;
  h(1,3)=z-z0_;
  
  TMatrixD cov_matr(dimension,dimension);
  for(int i=0; i<dimension; ++i)
    for(int j=0; j<dimension; ++j)
      cov_matr(i,j)=par_covariance_matrix_[i][j];
  
  TMatrixD V_hT(cov_matr, TMatrixD::kMultTranspose, h);
  //h*=V_hT;
  //return h;
  return TMatrixD(h, TMatrixD::kMult, V_hT);
}


T1T2Track::T1T2Track(double z0, const TVectorD & track_params_vector, 
      const TMatrixD &par_covariance_matrix, double chiSquared, int hemi, int det) 
  :  _chi2R(0), _chi2Phi(0), z0_(z0), chiSquared_(chiSquared), valid_(true)
{
  _detector = det;
  _hemisphere = hemi;
  for(int i=0; i<dimension; ++i)
  {
    track_params_vector_[i]=track_params_vector[i];
    for(int j=0; j<dimension; ++j)
    {
   
      par_covariance_matrix_[i][j]=par_covariance_matrix(i,j);
    }
  }

      bx_ = track_params_vector_[0]; 
    
      by_ =track_params_vector_[1];       
 

}

T1T2Track::T1T2Track( const TVectorD & track_params_vector, 
      const TMatrixD &par_covariance_matrix, double chiSquared, int hemi, int det) 
      : _chi2R(0), _chi2Phi(0), z0_(0), chiSquared_(chiSquared), valid_(true)
{
  _hemisphere = hemi;
  _detector = det;
  for(int i=0; i<dimension; ++i)
  {
    track_params_vector_[i]=track_params_vector[i];
    for(int j=0; j<dimension; ++j)
    {
    par_covariance_matrix_[i][j]=par_covariance_matrix(i,j);
      if(i==j){
      par_covariance_vector_[i]=par_covariance_matrix[i][j];
//      std::cout << "MATRICIONA " << par_covariance_matrix[i][j] << " " <<par_covariance_vector_[i]<< std::endl;
      }
    }

  }

  double ax = track_params_vector_[2];
  double bx = track_params_vector_[0];
  double ay = track_params_vector_[3];
  double by = track_params_vector_[1];



  _Z_at_Rmin = -(ax*bx + ay*by)/(ax*ax + ay*ay);
  _Rmin = sqrt( (ax*_Z_at_Rmin + bx)*(ax*_Z_at_Rmin + bx) + (ay*_Z_at_Rmin + by)*(ay*_Z_at_Rmin + by) );
  _theta = acos(1./sqrt( (ax*ax)+(ay*ay)+1));
  _eta = fabs(-log(tan(_theta/2.))) * (double)_hemisphere;



  ax=ax*hemi;
  ay=ay*hemi;


  _phi = atan(ay/ax);
  
  if((ax<0)&&(ay>0))
    _phi= _phi+3.14159265;
  
  if((ax<0)&&(ay<0))
    _phi= _phi+3.14159265;
  
  if((ax>0)&&(ay<0))
    _phi= _phi+2*3.14159265;



  bx_ = track_params_vector_[0]; 
  by_ =track_params_vector_[1];   
}

T1T2Track::T1T2Track(double z0, const TVectorD & track_params_vector, 
		     const TMatrixD &par_covariance_matrix, double chiSquared,double chiSquaredX,double chiSquaredY, int hemi, int det) 
  : z0_(z0), chiSquared_(chiSquared),chiSquaredX_(chiSquaredX),chiSquaredY_(chiSquaredY), valid_(true)
{
  _detector = det;
  _hemisphere = hemi;
  for(int i=0; i<dimension; ++i)
    {
      track_params_vector_[i]=track_params_vector[i];
      for(int j=0; j<dimension; ++j)
	{
   
	  par_covariance_matrix_[i][j]=par_covariance_matrix(i,j);
	}
    }
}




//Constructor utilized for T2 XY fit
T1T2Track::T1T2Track( const TVectorD & track_params_vector, 
		      const TMatrixD &par_covariance_matrix, double chiSquared,double chiSquaredX,double chiSquaredY,  int hemi, int det) 
  : t2_roadID(0), z0_(0), chiSquared_(chiSquared),chiSquaredX_(chiSquaredX),chiSquaredY_(chiSquaredY),  valid_(true)
{
  _hemisphere = hemi;
  _detector = det;
  for(int i=0; i<dimension; ++i)
    {
      track_params_vector_[i]=track_params_vector[i];
      for(int j=0; j<dimension; ++j)
	{
	  par_covariance_matrix_[i][j]=par_covariance_matrix(i,j);
	  if(i==j){
	    par_covariance_vector_[i]=par_covariance_matrix[i][j];
	    //      std::cout << "MATRICIONA " << par_covariance_matrix[i][j] << " " <<par_covariance_vector_[i]<< std::endl;
	  }
	}

    }

  double ax = track_params_vector_[2];
  double bx = track_params_vector_[0];
  double ay = track_params_vector_[3];
  double by = track_params_vector_[1];
  //commented on 24/9/2010
  /*
  if(_detector == 2){
    ax=ax*hemi;
    ay=ay*hemi;
  }
  */
  _Z_at_Rmin = -(ax*bx + ay*by)/(ax*ax + ay*ay);
  _Rmin = sqrt( (ax*_Z_at_Rmin + bx)*(ax*_Z_at_Rmin + bx) + (ay*_Z_at_Rmin + by)*(ay*_Z_at_Rmin + by) );
  _theta = acos(1./sqrt( (ax*ax)+(ay*ay)+1));
  _eta = fabs(-log(tan(_theta/2.))) * (double)_hemisphere;


  ax=ax*hemi;
  ay=ay*hemi;
  

  _phi = atan(ay/ax);
  
  if((ax<0)&&(ay>0))
    _phi= _phi+3.14159265;
  
  if((ax<0)&&(ay<0))
    _phi= _phi+3.14159265;
  
  if((ax>0)&&(ay<0))
    _phi= _phi+2*3.14159265;

  bx_ = track_params_vector_[0]; 
  by_ =track_params_vector_[1];     

  _t2GhostId=0;

}






//Constructor utilized for T2 XY AND RZ fit
T1T2Track::T1T2Track( const TVectorD & track_params_vector, 
		      const TMatrixD &par_covariance_matrix, double chiSquared,double chiSquaredX,double chiSquaredY,  int hemi, int det,
		      double TanthetaRZ, double bRZ, double  phiRZ, double  e_TanthetaRZ, double e_bRZ, double e_phiRZ, double chi2R, double  chi2Phi, int theT2RoadID, unsigned int numHit_t2, unsigned int numCl1HitHit_t2, unsigned int numStripOnlyHit_t2, unsigned int numPadOnlyHit_t2

		      )   : t2_roadID(0),z0_(0), chiSquared_(chiSquared),chiSquaredX_(chiSquaredX),chiSquaredY_(chiSquaredY),  valid_(true)
{
  _hemisphere = hemi;
  _detector = det;
  for(int i=0; i<dimension; ++i)
    {
      track_params_vector_[i]=track_params_vector[i];
      for(int j=0; j<dimension; ++j)
	{
	  par_covariance_matrix_[i][j]=par_covariance_matrix(i,j);
	  if(i==j){
	    par_covariance_vector_[i]=par_covariance_matrix[i][j];
	    //      std::cout << "MATRICIONA " << par_covariance_matrix[i][j] << " " <<par_covariance_vector_[i]<< std::endl;
	  }
	}

    }
  /*

  inline double X0() const {return track_params_vector_[0];}
  inline double X0Sigma() const {return sqrt(par_covariance_vector_[0]);}
  inline double Y0() const {return track_params_vector_[1];}
  inline double Y0Sigma() const {return sqrt(par_covariance_vector_[1]);}
  inline double Z0() const {return z0_;}
  inline void Z0(double z0) {z0_=z0;}
  inline double GetTx() const {return track_params_vector_[2];}
  inline double GetTxSigma() const {return sqrt(par_covariance_vector_[2]);}
  inline double GetTy() const {return track_params_vector_[3];}
  inline double GetTySigma() const {return sqrt(par_covariance_vector_[3]);}

*/
  /*
   FittedParam(0)=b_xz;
   FittedParam(1)=b_yz;
   FittedParam(2)=a_xz;
   FittedParam(3)=a_yz;
  */

  //ax.bx.ay.by
  double ax = track_params_vector_[2];
  double bx = track_params_vector_[0];
  double ay = track_params_vector_[3];
  double by = track_params_vector_[1];
  //commented on 24/9/2010
  /*
  if(_detector == 2){
    ax=ax*hemi;
    ay=ay*hemi;
  }
  */
  _Z_at_Rmin = -(ax*bx + ay*by)/(ax*ax + ay*ay);
  _Rmin = sqrt( (ax*_Z_at_Rmin + bx)*(ax*_Z_at_Rmin + bx) + (ay*_Z_at_Rmin + by)*(ay*_Z_at_Rmin + by) );
  _theta = acos(1./sqrt( (ax*ax)+(ay*ay)+1));
  _eta = fabs(-log(tan(_theta/2.))) * (double)_hemisphere;

  ax=ax*hemi;
  ay=ay*hemi;
  _phi = atan(ay/ax);
  
  if((ax<0)&&(ay>0))
    _phi= _phi+3.14159265;
  
  if((ax<0)&&(ay<0))
    _phi= _phi+3.14159265;
  
  if((ax>0)&&(ay<0))
    _phi= _phi+2*3.14159265;

  bx_ = track_params_vector_[0]; 
  by_ =track_params_vector_[1]; 

 
  
  t2_roadID=  theT2RoadID;

  //T2-only field for the RZ fitting
  
  _TanthetaRZ = TanthetaRZ;
  _bRZ = bRZ;    
  _phiRZ= phiRZ;
  _e_TanthetaRZ=e_TanthetaRZ;
  _e_bRZ=e_bRZ;  
  _e_phiRZ =e_phiRZ;
  _chi2R = chi2R;
  _chi2Phi= chi2Phi;
  _etaRZ = fabs(-log(fabs(TanthetaRZ)/2.)) * (double)_hemisphere;  //approximated.
 
  _ax_rz = TanthetaRZ * cos( phiRZ );
  _ay_rz = TanthetaRZ * sin( phiRZ );
  _bxRZ = bRZ * cos( phiRZ ); 
  _byRZ = bRZ * sin( phiRZ );
  _Z_at_RminRZ = -bRZ/TanthetaRZ;
  z0_=-bRZ/TanthetaRZ;


  _numPadOnlyHit_t2 = numPadOnlyHit_t2;  
  _numStripOnlyHit_t2 = numStripOnlyHit_t2;
  _numCl1HitHit_t2 = numCl1HitHit_t2;
  _numHit_t2 = numHit_t2;
  
  _t2GhostId=0;
}

//Constructor from T2TrackProducer RZ fit
T1T2Track::T1T2Track(const TVectorD & track_params_vector, 
      const TMatrixD &par_covariance_matrix, double a_rz, double b_rz, double phim,double e_a_rz,
double e_b_rz,double e_phim,double chi2,double chi2R,double chi2Phi, double normchi2red, int hemisphere, int det): t2_roadID(0), _chi2R(chi2R), _chi2Phi(chi2Phi), z0_(0), chiSquared_(chi2),chiSquaredoverN_(normchi2red),   valid_(true)
{
  _detector=det;
  _hemisphere=hemisphere;  
  _theta =atan(fabs(a_rz));
  _eta = fabs(-log(tan(_theta/2.))) * (double)hemisphere;
  _phi = phim;
  // _Z_at_Rmin=-b_rz/tan(_theta);//until 21-12-09
  _Z_at_Rmin=-b_rz/atan(a_rz);
  _Rmin= 0.0;

  _TanthetaRZ=a_rz;
  _e_TanthetaRZ=e_a_rz;
 
  _bRZ  = b_rz;
  _e_bRZ = e_b_rz; 

  _phiRZ=phim;
  _e_phiRZ=e_phim;

  for(int i=0; i<dimension; ++i)
  {
    track_params_vector_[i]=track_params_vector[i];
    for(int j=0; j<dimension; ++j)
    {
    par_covariance_matrix_[i][j]=par_covariance_matrix(i,j);
      if(i==j){
      par_covariance_vector_[i]=par_covariance_matrix[i][j];
      }
    }
  }

  
  bx_ = a_rz * cos( phim ); 
    
  by_ = a_rz * sin( phim ); 
  
  _t2GhostId=0;
  // cout<<"Inside special constructor: "<<"  a_rz= "<< a_rz <<"  theta="<<_theta << "  _eta=" <<_eta <<"  _phi(rad)="<<_phi  <<endl;
}




TVectorD T1T2Track::ParameterVector() const 
{
  TVectorD v(dimension);
  
  for(int i=0; i<dimension; ++i)
    v[i]=track_params_vector_[i];
      
  return v;
}



void T1T2Track::ParameterVector(const TVectorD & track_params_vector)
{
  for(int i=0; i<dimension; ++i)
    track_params_vector_[i]=track_params_vector[i];
}




TMatrixD T1T2Track::CovarianceMatrix() const 
{
  TMatrixD m(dimension,dimension);
  
  for(int i=0; i<dimension; ++i)
    for(int j=0; j<dimension; ++j)
      m(i,j)=par_covariance_matrix_[i][j];
      
  return m;
}



void T1T2Track::CovarianceMatrix(const TMatrixD &par_covariance_matrix)
{
  for(int i=0; i<dimension; ++i)
    for(int j=0; j<dimension; ++j)
      par_covariance_matrix_[i][j]=par_covariance_matrix(i,j);
}


// Count number of common hits with other track, based on hit position
unsigned int T1T2Track::NCommonHits(T1T2Track& otherTrack)
{
  unsigned int nCommonHits = 0;
  if (_detector == 1) {  // currently implemented only for T1
    for (vector<T1RecHitGlobal>::const_iterator itH1 = track_hits_vector_1.begin(); itH1 != track_hits_vector_1.end(); itH1++)
      for (unsigned int iH2 = 0; iH2 < otherTrack.GetHitEntries(); iH2++) {
	GlobalPoint gp1 = (*itH1).GlobalPosition();
	GlobalPoint gp2 = otherTrack.GetHitT1(iH2).GlobalPosition();
	if (fabs(gp1.x() - gp2.x()) + fabs(gp1.y() - gp2.y()) + fabs(gp1.z() - gp2.z()) < 1.e-3)
	  nCommonHits++;
      }
  }
  return nCommonHits;
}

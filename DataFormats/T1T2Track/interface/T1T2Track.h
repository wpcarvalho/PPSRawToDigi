/*
  Modified by Fabrizio Ferro - INFN Genova for TOTEM
*/
#ifndef DataFormats_T1T2Track_h
#define DataFormats_T1T2Track_h

#include <Rtypes.h>
#include <TMath.h>
#include "TROOT.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include <DataFormats/T1Road/interface/T1RecHitGlobal.h>
#include <DataFormats/T2Road/interface/T2Road.h>
#include <DataFormats/T2Road/interface/T2RoadCollection.h>




class T1T2Track    //x = x0+tx*(z-z0) y = ...
{
 public:
  
  double bx_;
  double by_;  
  unsigned int t2_roadID;  //if assigned is strictly > 0.
  enum { dimension = 4 };
  // parameter vector size
  enum { covarianceSize = dimension * dimension };
  // covariance matrix size
  T1T2Track() : /*track_params_vector_(4),*/_chi2R(0), _chi2Phi(0),  z0_(0),
    /*par_covariance_matrix_(4,4),*/ chiSquared_(0), valid_(false) {
      
  }
  T1T2Track(int det) : /*track_params_vector_(4),*/  _chi2R(0), _chi2Phi(0), z0_(0),
    /*par_covariance_matrix_(4,4),*/ chiSquared_(0), valid_(true) {
    _detector = det;
  }

  T1T2Track(double z0, const TVectorD & track_params_vector, 
            const TMatrixD &par_covariance_matrix, double chiSquared, double chiSquaredX, double chiSquaredY,int hemi, int detector);
  T1T2Track( const TVectorD & track_params_vector, 
	     const TMatrixD &par_covariance_matrix, double chiSquared, int hemi,int detector);
  T1T2Track(double z0, const TVectorD & track_params_vector, 
	    const TMatrixD &par_covariance_matrix, double chiSquared, int hemi, int detector);



  //Constructor utilized for T2 XY fit  
  T1T2Track( const TVectorD & track_params_vector, 
	     const TMatrixD &par_covariance_matrix, double chiSquared, double chiSquaredX, double chiSquaredY,int hemi,int detector);


  //Constructor from T2TrackProducer RZ fit
  T1T2Track(const TVectorD & track_params_vector, const TMatrixD &par_covariance_matrix,double a_rz, double b_rz, double phim,double e_a_rz, 
	    double e_b_rz,double e_phim,double chi2,double chi2R,double chi2Phi,double normchi2red,int hemisphere, int detector); 

  
  //Constructor utilized for T2 XY AND RZ fit
  T1T2Track(const TVectorD & track_params_vector, const TMatrixD &par_covariance_matrix, 
	    double chiSquared,double chiSquaredX,double chiSquaredY,  int hemi, int det, double TanthetaRZ, double bRZ, double  phiRZ, double  e_TanthetaRZ, double e_bRZ, double e_phiRZ, double chi2R, double  chi2Phi, int theT2RoadID, unsigned int numHit_t2, unsigned int numCl1HitHit_t2, unsigned int numStripOnlyHit_t2, unsigned int numPadOnlyHit_t2);
		      



  virtual ~T1T2Track() {}
  inline unsigned int GetHitEntries() const {
    if(_detector == 1)
      return track_hits_vector_1.size();
    else
      return track_hits_vector_2.thisRoad.size();
  }

  inline const T1RecHitGlobal & GetHitT1(int i) const {
    return track_hits_vector_1[i];
  }

  inline const T2Hit & GetHitT2(int i) const {
    return track_hits_vector_2.thisRoad[i];
  }
  
  inline std::vector<T2Hit> GetT2TrackHits(){
    
    std::vector<T2Hit> t2hits;
    unsigned int roadsize=track_hits_vector_2.thisRoad.size();
    for(unsigned int i=0;i<roadsize;i++)
      t2hits.push_back(track_hits_vector_2.thisRoad[i]);

    return t2hits;
  }
  
  inline void AddHit(const T1RecHitGlobal &hit) {track_hits_vector_1.push_back(hit);}
  inline void AddHit(const T2Hit hit) {track_hits_vector_2.thisRoad.push_back(hit);}

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

  inline double Z_at_Rmin() const {return _Z_at_Rmin;}
  inline double Rmin() const {return _Rmin;}
  inline double Theta() const {return _theta;}
  inline double Eta() const {return _eta;}
  inline double Phi() const {return _phi;}

  inline TVector3 GetDirectionVector() const {TVector3 vect(GetTx(), GetTy(), 1); vect.SetMag(1.0); return vect;}
  TVectorD ParameterVector() const;
  void ParameterVector(const TVectorD & track_params_vector);
  TMatrixD CovarianceMatrix() const;
  void CovarianceMatrix(const TMatrixD &par_covariance_matrix);
  inline double ChiSquared() const {return chiSquared_;}
  inline void ChiSquared(double & chiSquared) {chiSquared_=chiSquared;}
  inline double ChiSquaredR() const {return _chi2R;}
  inline double ChiSquaredPhi() const {return _chi2Phi;}
  inline double ChiSquaredX() const {return chiSquaredX_;}
  inline void ChiSquaredX(double & chiSquared) {chiSquaredX_=chiSquared;}

  inline double ChiSquaredY() const {return chiSquaredY_;}
  inline void ChiSquaredY(double & chiSquared) {chiSquaredY_=chiSquared;}
  inline double ChiSquaredOverN() const {
    if(_detector == 1)
      return chiSquared_/(float)(2*(track_hits_vector_1.size()-2)); 
    else
      return chiSquaredoverN_;
  }
  inline double ChiSquaredXOverN() const {
    if(_detector == 1)
      return chiSquaredX_/(float)((track_hits_vector_1.size()-2));
    else
      return -1;
  }
    
  inline double ChiSquaredYOverN() const {
    if(_detector == 1)
      return chiSquaredY_/(float)((track_hits_vector_1.size()-2));
    else
      return -1;
  }
    
  inline TVector2 GetTrackPoint(double z) const //returns x,y vector
    {
      double delta_z = z-z0_;
      return TVector2(track_params_vector_[0]+track_params_vector_[2]*delta_z, 
		      track_params_vector_[1]+track_params_vector_[3]*delta_z);
    }
  inline TVector3 TrackCentrePoint() {return TVector3(track_params_vector_[0], 
						      track_params_vector_[1], z0_);}
    
  TMatrixD TrackPointInterpolationCovariance(double z) const;
  inline bool IsValid() {return valid_;}
  inline void IsValid(bool valid) {valid_=valid;}
    
  inline void Reset() {track_hits_vector_1.clear();
  track_hits_vector_2.thisRoad.clear();
  }
    
  inline void SetDet(int det){_detector = det;}
  inline int GetDet() const {return _detector;}

  unsigned int NCommonHits(T1T2Track& otherTrack);

  //T2-only field:  this parameter are calculated using fits in RZ plane 
  //while the T2 track parameter fitted in the XY plane can be taken with 
  //the same method utilized by T1.

  double _TanthetaRZ;
  double _bRZ;               
  double _bxRZ; 
  double _byRZ; 
  double _phiRZ;
  double _ax_rz;
  double _ay_rz;
  double _bxRZ_; 
  double _byRZ_;
  double _e_TanthetaRZ;
  double _e_bRZ;  
  double _e_phiRZ;  
  double _etaRZ;    
  double _chi2R;
  double _chi2Phi; 
  double _Z_at_RminRZ; 
  unsigned int _numHit_t2, _numPadOnlyHit_t2, _numStripOnlyHit_t2, _numCl1HitHit_t2;
  int _t2GhostId;          //0: no ghost found, otherwise track with the same index can be the ghost of some other 
                           //trk having the same index. Warning: if 4 trks have are given by 2 reals and 2 ghosts I will write   
                           //I will put the same index for the same column. But inside the same column I will put positive  
                           //index for the track with more cl1Hit (more likely to be the real one) and negative for the ghost-candidate.

  
 private:
  std::vector<T1RecHitGlobal> track_hits_vector_1;
  T2Road track_hits_vector_2;

  //TVectorD track_params_vector_;  //x0, y0, tx, ty; x = x0+tx*(z-z0) ...
  double track_params_vector_[dimension];
  double z0_;
  //TMatrixD par_covariance_matrix_;  //x0, y0, tx, ty...
  double par_covariance_matrix_[dimension][dimension];
  double par_covariance_vector_[dimension];
  double chiSquared_;

  double chiSquaredoverN_;  //Need for T2
       
  double chiSquaredX_;
  double chiSquaredY_;

  double _Z_at_Rmin;
  double _Rmin;
  double _theta;
  double _phi;
  double _eta;
  int _hemisphere;
  int _detector; // 1 = T1, 2 = T2

  
  bool valid_;
};



#endif

/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Hubert Niewiadomski
 *   Jan Kašpar (jan.kaspar@gmail.com)
 *
 ****************************************************************************/

#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrack.h"

//----------------------------------------------------------------------------------------------------

TMatrixD RPFittedTrack::TrackPointInterpolationCovariance(double z) const
{
  TMatrixD h(2,4);
  h(0,0)=1;
  h(1,1)=1;
  h(0,2)=z-z0_;
  h(1,3)=z-z0_;
  
  TMatrixD cov_matr(dimension,dimension);
  for(int i=0; i<dimension; ++i)
    for(int j=0; j<dimension; ++j)
      cov_matr(i,j)=CovarianceMatrixElement(i,j);
  
  TMatrixD V_hT(cov_matr, TMatrixD::kMultTranspose, h);
  //h*=V_hT;
  //return h;
  return TMatrixD(h, TMatrixD::kMult, V_hT);
}

//----------------------------------------------------------------------------------------------------

RPFittedTrack::RPFittedTrack(double z0, const TVectorD & track_params_vector, 
      const TMatrixD &par_covariance_matrix, double chiSquared) 
      : z0_(z0), chiSquared_(chiSquared), valid_(true), u_id_(0), v_id_(0), sourceTrackCandidateValid(false)
{
  for(int i=0; i<dimension; ++i)
  {
    track_params_vector_[i]=track_params_vector[i];
    for(int j=0; j<dimension; ++j)
    {
      CovarianceMatrixElement(i,j)=par_covariance_matrix(i,j);
    }
  }
}

//----------------------------------------------------------------------------------------------------

TVectorD RPFittedTrack::ParameterVector() const 
{
  TVectorD v(dimension);
  
  for(int i=0; i<dimension; ++i)
    v[i]=track_params_vector_[i];
      
  return v;
}

//----------------------------------------------------------------------------------------------------

void RPFittedTrack::ParameterVector(const TVectorD & track_params_vector)
{
  for(int i=0; i<dimension; ++i)
    track_params_vector_[i]=track_params_vector[i];
}

//----------------------------------------------------------------------------------------------------

TMatrixD RPFittedTrack::CovarianceMatrix() const 
{
  TMatrixD m(dimension,dimension);
  
  for(int i=0; i<dimension; ++i)
    for(int j=0; j<dimension; ++j)
      m(i,j)=CovarianceMatrixElement(i,j);
      
  return m;
}

//----------------------------------------------------------------------------------------------------

void RPFittedTrack::CovarianceMatrix(const TMatrixD &par_covariance_matrix)
{
  for(int i=0; i<dimension; ++i)
    for(int j=0; j<dimension; ++j)
      CovarianceMatrixElement(i,j)=par_covariance_matrix(i,j);
}

//----------------------------------------------------------------------------------------------------

bool operator< (const RPFittedTrack &l, const RPFittedTrack &r)
{
  if (l.z0_ < r.z0_)
    return true;
  if (l.z0_ > r.z0_)
    return true;

  for (int i = 0; i < RPFittedTrack::dimension; ++i)
  {
    if (l.track_params_vector_[i] < r.track_params_vector_[i])
      return true;
    if (l.track_params_vector_[i] > r.track_params_vector_[i])
      return true;
  }

  return false;
}

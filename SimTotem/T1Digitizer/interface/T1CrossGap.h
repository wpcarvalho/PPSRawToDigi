#ifndef T1_CROSS_GAP_H
#define T1_CROSS_GAP_H

/** \class T1CrossGap
 * Class used to represent one T1 gas gap crossing by a charged track.
 *
 * \author Tim Cox
 *
 * This is used by T1T1GasCollisions in the digitization of the T1s.
 * Actually this may NOT model the whole gas gap, but just the path length
 * corresponding to a PSimHit in the gap (i.e. from the PSimHit 'entry point'
 * to its 'exit point'). PSimHit's can have a pretty evanescent existence
 * and I can't pretend to have a full understanding of their full range
 * of characteristics. Yet another thing 'to be studied'.
 *
 */

#include <vector>
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
class HepParticleData;

class T1CrossGap {
 public:

  /**
   * iam = particle type in PDG code.
   * mom = momentum of particle
   * gap = space std::vector representing the hit entry and exit
   */
  T1CrossGap(int iam, float mom, LocalVector gap);
  ~T1CrossGap(){};

  std::vector<LocalPoint> ionClusters() const { return clusters; }
  int noOfClusters() const { return clusters.size(); }
  std::vector<int> electrons() const { return electronsInClusters; }
  int noOfElectrons() const { return electronsInClusters.size(); }
  std::vector<double> stepLengths() const { return steps; }
  int noOfSteps() const { return steps.size(); }
  std::vector<float> eLossPerStep() const { return elosses; }
  int noOfElosses() const { return elosses.size(); }

  void addCluster(LocalPoint here) { clusters.push_back( here ); }
  void addElectrons(int nelec = 1) { electronsInClusters.push_back( nelec ); }
  void addElectronToBack() { ++electronsInClusters.back(); }

  void addStep( double step ) { steps.push_back( step ); }
  void addEloss( float eloss ) { elosses.push_back( eloss ); }
  
  double logGamma( double mass, float momentum );
  double logGamma(){ return loggam; }
  double beta2() const { return theBeta2; }
  double gamma() const { return theGamma; }
  LocalVector gapVector() const { return theGap; }
  LocalVector unitVector() const { return theGap.unit(); }
  float length() const { return theGap.mag(); }

 private:
  int setParticle(int);

  double theBeta2; // Lorentz beta^2
  double theGamma; // Lorentz gamma
  double loggam; 
  LocalVector theGap;

  const HepParticleData * theParticleData; // buffer for incident charged particle

  std::vector<LocalPoint> clusters;
  std::vector<int> electronsInClusters;
  std::vector<double> steps;
  std::vector<float> elosses;

};

#endif

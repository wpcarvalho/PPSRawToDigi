/**
\mainpage
\section Links Links
  - \ref RPSimuReco
  - \ref TotemReadout
**/


/**\defgroup RPSimuReco RP simulation and reconstruction chain
\brief Description of RP simulation and reconstruction classes.

\section s1 Step 1: Event generation
- classes: e.g. ElegentSource
- output: HepMCProduct

\section s2 Step 2: Smearing
- classes: e.g. EnergySmearing, VertexSmearing
- output: HepMCProduct

\section s3 Step 3: Geant simulation
- classes: 
- output: PSimHit mix collection

\section s4 Step 4: Digitization
- classes: see \ref TotemDigiProduction
- output: DetSetVector<RPStripDigi>, DetSetVector<RPDetTrigger> (see RPStripDigi and RPDetTrigger)

\section s5 Step 5: Clusterization
- classes: see \ref RPClusterization
- output: DetSetVector<RPDigCluster> (see RPDigCluster)

\section s6 Step 6: Reco production
- classes: see \ref RPRecoHitProduction
- output: DetSetVector<RPRecoHit> (see RPRecoHit)

\section s7 Step 7: Track candidates search
- classes: RPSingleCandidateTrackFinder or RPNonParallelTrackCandidateFinder
- output: RPTrackCandidateCollection

\section s8 Step 8: Track fitting
- classes: RPTrackCandidateCollectionFitter
- output: RPFittedTrackCollection

\section s9 Step 9: Elastic or inelastic reconstruction
- classes (elastic): \ref ElasticReconstruction
- output (elastic): RPRecoElasticEvent
 */

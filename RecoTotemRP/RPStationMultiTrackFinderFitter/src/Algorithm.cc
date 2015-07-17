/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Jakub Sawicki (jakub.kuba.sawicki@gmail.com)
 *
 ****************************************************************************/

#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/Algorithm.h"

namespace RPStationMultiTrackFinderFitter {

Algorithm::Algorithm(
        const edm::ParameterSet& ps,
        GeometryHelper &geometry,
        IClusterBuilder &clusterBuilder) :
    verbosity(ps.getUntrackedParameter<unsigned int>("verbosity")),
    geometry(geometry),
    clusterBuilder(clusterBuilder),
    overlapsRemoval(ps),
    filter_dist_cut(ps.getParameter<double>("filter_dist_cut")),
    track_pattern_assoc_cut(ps.getParameter<double>("track_pattern_assoc_cut")),
    z_h(ps.getParameter<double>("z_h")),
    minRPsPerStation(ps.getParameter<unsigned int>("minRPsPerStation")),
    verboseEnabled(true)
{ 
}

//----------------------------------------------------------------------------------------------------

Algorithm::~Algorithm()
{
}

//----------------------------------------------------------------------------------------------------

void Algorithm::reconstruct(
        vector<Cluster> *output,
        const RPRecognizedPatternsCollection *input)
{
    // clusters
    vector<Cluster> clusters;
    clusterBuilder.build(clusters, *input);

    // TODO: redundant ??
    // rough selection of clusters
    vector<Cluster> selectedClusters;
    selectClusters(selectedClusters, clusters);

    // filter clusters
    vector<Cluster> filteredClusters;
    filterClusters(filteredClusters, selectedClusters, *input);

    // remove overlaps
    removeOverlaps(*output, filteredClusters, *input);
}

//----------------------------------------------------------------------------------------------------

void Algorithm::selectClusters(vector<Cluster> &out, const vector<Cluster> &in)
{
    if (verbose(2))
        printf(">> Algorithm::selectClusters\n");

    for (unsigned int ci = 0; ci < in.size(); ci++)
    {
        const Cluster &c = in[ci];
        
        bool selected = selectCluster(c);

        if (selected)
            out.push_back(c);

        if (verbose(2))
        {
            printf("cluster #%u: RPs = %lu | ", ci, c.rpIdCount.size());
            for (map<RPId, unsigned int>::const_iterator rpIt = c.rpIdCount.begin(); rpIt != c.rpIdCount.end(); ++rpIt)
                printf("%u -> %u, ", rpIt->first, rpIt->second);
            printf(" | selected = %u [%u]\n", selected, minRPsPerStation);
        }

    }

    if (verbose(1))
        printf("\t--> selected clusters: %lu\n", out.size());
}

//----------------------------------------------------------------------------------------------------

bool Algorithm::selectCluster(const Cluster &cluster)
{
    if (cluster.rpIdCount.size() < minRPsPerStation)
        return false;

    for (map<RPId, unsigned int>::const_iterator rpIt = cluster.rpIdCount.begin(); rpIt != cluster.rpIdCount.end(); ++rpIt)
    {
        if (rpIt->second < minRPsPerStation - 1)
            return false;
    }

    return true;
}

//----------------------------------------------------------------------------------------------------

void Algorithm::filterClusters(
        vector<Cluster> &out, const vector<Cluster> &in,
        const RPRecognizedPatternsCollection &rpPatColl)
{
    if (verbose(2))
        printf(">> Algorithm::filterClusters\n");

    // discard clusters that can not be associated with input hits (= 1-RP U/V patterns) 
    for (unsigned int ci = 0; ci < in.size(); ci++)
    {
        const Cluster &cluster = in[ci];

        if (verbose(1))
        {
            printf("* selected cluster #%u: ", ci);
            Track t = cluster.GetCenterTrack();
            t.z0 = z_h;
            t.Print();
        }

        if (clusterMatchesPatterns(cluster, rpPatColl))
        {
            if (verbose(2))
                printf("\t\tAccepted: #%u -> #%lu\n", ci, out.size());

            out.push_back(cluster);
        }
    }

    if (verbose(1))
        printf("\t--> filtered clusters: %lu\n", out.size());
}

//----------------------------------------------------------------------------------------------------

bool Algorithm::clusterMatchesPatterns(
        const Cluster &cluster,
        const RPRecognizedPatternsCollection &rpPatColl)
{
    unsigned int compatibleUPatterns = 0;
    unsigned int compatibleVPatterns = 0;
    unsigned int totalRPsCompatible = 0;

    for (RPRecognizedPatternsCollection::const_iterator rpIt = rpPatColl.begin(); rpIt != rpPatColl.end(); ++rpIt)
    {
        RPId rpId = rpIt->first;
        const vector<RPRecognizedPatterns::Line> &uPatterns = rpIt->second.uLines;
        const vector<RPRecognizedPatterns::Line> &vPatterns = rpIt->second.vLines;

        RPInfo rpInfo = geometry.getRPInfo(rpId);

        if (verbose(2))
            printf("\tRP: %u\n", rpId);

        // check if track candidate should be visible in the RP
        if (!trackHitsRP(cluster.GetCenterTrack(), rpId))
        {
            if (verbose(2))
                printf("\t\tpoint outside RP anyway, skip\n");

            continue;
        }

        // candidate visible
        totalRPsCompatible++;
        
        // U patterns
        if (verbose(2))
            printf("\t\tU\n");
        pair<unsigned int, double> assoc = assocClusterWithPattern(cluster, rpId, uPatterns, U);
        if (assoc.first != uPatterns.size() && assoc.second < filter_dist_cut)
            compatibleUPatterns++;

        // U patterns
        if (verbose(2))
            printf("\t\tV\n");
        assoc = assocClusterWithPattern(cluster, rpId, vPatterns, V);
        if (assoc.first != vPatterns.size() && assoc.second < filter_dist_cut)
            compatibleVPatterns++;
    }

    return (totalRPsCompatible >= minRPsPerStation && compatibleVPatterns == totalRPsCompatible
      && compatibleUPatterns == totalRPsCompatible);
}

//----------------------------------------------------------------------------------------------------

void Algorithm::removeOverlaps(
        vector<Cluster> &out, const vector<Cluster> &in,
        const RPRecognizedPatternsCollection &rpPatColl)
{
    if (verbosity > 2)
        printf(">> Algorithm::removeOverlaps\n");

    vector<unsigned int> hitCount;
    vector<vector<unsigned int> > cluHits;
    //for (unsigned int ci = 0; ci < in.size(); ci++) {
        //vector<unsigned int> vec;
        //cluHits.push_back(vec);
    //}
    associateClustersWithHits(cluHits, in, rpPatColl);

    // must be in the order associateclusterswithhits uses
    for (RPRecognizedPatternsCollection::const_iterator rpIt = rpPatColl.begin(); rpIt != rpPatColl.end(); ++rpIt)
    {
        hitCount.push_back(rpIt->second.uLines.size());
        hitCount.push_back(rpIt->second.vLines.size());
    }

    overlapsRemoval.run(out, in, hitCount, cluHits);
}

//----------------------------------------------------------------------------------------------------

bool Algorithm::trackHitsRP(Track track, RPId rpID)
{
    RPInfo rpInfo = geometry.getRPInfo(rpID);

    double z_p = rpInfo.cz,
           x = track.ax * (z_p - z_h) + track.bx,
           y = track.ay * (z_p - z_h) + track.by;
    return !geometry.isPointOutside(rpID, CLHEP::Hep3Vector(x, y, z_p));
}

//----------------------------------------------------------------------------------------------------

pair<unsigned int, double> Algorithm::assocClusterWithPattern(
        const Cluster &cluster,
        RPId rpID,
        const vector<RPRecognizedPatterns::Line> &patterns,
        Direction direction)
{
    RPInfo rpInfo = geometry.getRPInfo(rpID);

    double d2_best = 1E100;
    unsigned int d2_best_pi = patterns.size();

    double dx = 0, dy = 0;

    switch (direction)
    {
        case U:
            dx = rpInfo.u_dx;
            dy = rpInfo.u_dy;
            break;
        case V:
            dx = rpInfo.v_dx;
            dy = rpInfo.v_dy;
            break;
    }

    double z_p = rpInfo.cz;
    TMatrixD P(1, 4);
    P(0, 0) = (z_p - z_h) * dx;
    P(0, 1) = dx;
    P(0, 2) = (z_p - z_h) * dy;
    P(0, 3) = dy;

    TMatrixD PT(TMatrixD::kTransposed, P);

    TVectorD cpv = P * cluster.center;
    TMatrixD cpV = P * cluster.sum_W_i * PT;

    double cp = cpv(0);
    double cp_u = sqrt(cpV(0, 0));

    // patterns
    for (unsigned int pi = 0; pi < patterns.size(); pi++)
    {
        double p = patterns[pi].b;
        double p_u = 20E-3;

        double diff = p - cp;
        double d2 = diff*diff / (p_u*p_u + cp_u*cp_u);

        if (d2 < d2_best)
        {
            d2_best = d2;
            d2_best_pi = pi;
        }

        if (verbose(2))
            printf("\t\t\tp=%.3f +- %.3f, cp=%.3f +- %.3f | d2 = %.2f\n", p, p_u, cp, cp_u, d2);
    }

    return make_pair(d2_best_pi, d2_best);
}

//----------------------------------------------------------------------------------------------------

vector<unsigned int> Algorithm::associateClusterWithHits(
        const Cluster &cluster,
        const RPRecognizedPatternsCollection &rpPatColl)
{
    vector<unsigned int> cluHits;

    for (RPRecognizedPatternsCollection::const_iterator rpIt = rpPatColl.begin(); rpIt != rpPatColl.end(); ++rpIt)
    {
        unsigned int rpId = rpIt->first;
        const vector<RPRecognizedPatterns::Line> &uPatterns = rpIt->second.uLines;
        const vector<RPRecognizedPatterns::Line> &vPatterns = rpIt->second.vLines;

        RPInfo rpInfo = geometry.getRPInfo(rpId);

        unsigned int uAssoc = uPatterns.size(),
                     vAssoc = vPatterns.size();

        // check if track candidate should be visible in the RP
        if (trackHitsRP(cluster.GetCenterTrack(), rpId))
        {
            verboseEnable(false);
            pair<unsigned int, double> assoc = assocClusterWithPattern(cluster, rpId, uPatterns, U);
            if (assoc.second < track_pattern_assoc_cut)
                uAssoc = assoc.first;

            assoc = assocClusterWithPattern(cluster, rpId, vPatterns, V);
            if (assoc.second < track_pattern_assoc_cut)
                vAssoc = assoc.first;

            verboseEnable(true);
        }

        cluHits.push_back(uAssoc);
        cluHits.push_back(vAssoc);

        if (verbose(2))
        {
            if (uAssoc != uPatterns.size())
                printf("%u ", uAssoc);
            else
                printf("X ");

            if (vAssoc != vPatterns.size()) 
                printf("%u ", vAssoc);
            else
                printf("X ");
        }
    }
    
    if (verbose(2))
        printf("\n");

    return cluHits;
}

//----------------------------------------------------------------------------------------------------

void Algorithm::associateClustersWithHits(
        vector<vector<unsigned int> > &cluHits, const vector<Cluster> &clusters,
        const RPRecognizedPatternsCollection &rpPatColl)
{
    for (unsigned int ci = 0; ci < clusters.size(); ci++)
    {
        const Cluster &cluster = clusters[ci];

        if (verbose(2))
            printf("cluster #%u: ", ci);

        cluHits.push_back(associateClusterWithHits(cluster, rpPatColl));
    }
}

//----------------------------------------------------------------------------------------------------

bool Algorithm::verbose(unsigned int level)
{
    return verboseEnabled && (verbosity > level);
}

//----------------------------------------------------------------------------------------------------

void Algorithm::verboseEnable(bool enable)
{
    verboseEnabled = enable;
}

} // namespace RPStationMultiTrackFinderFitter

/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Jakub Sawicki (jakub.kuba.sawicki@gmail.com)
 *
 ****************************************************************************/

#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/ClusterBuilderSubsets.h"
#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/Cluster.h"

#include <vector>

namespace RPStationMultiTrackFinderFitter {

ClusterBuilderSubsets::ClusterBuilderSubsets(const edm::ParameterSet &ps, GeometryHelper &geometry) :
    IClusterBuilder(ps, geometry),
    ax_mean(ps.getParameter<double>("ax_mean")),
    ax_cut(ps.getParameter<double>("ax_cut")),
    ay_mean(ps.getParameter<double>("ay_mean")),
    ay_cut(ps.getParameter<double>("ay_cut")),
    z_h(ps.getParameter<double>("z_h")),
    uv_unc(ps.getParameter<double>("uv_unc")),
    edge_create_th(ps.getParameter<double>("edge_create_th")),
    minRPsPerStation(ps.getParameter<unsigned int>("minRPsPerStation"))
{
}

//----------------------------------------------------------------------------------------------------

void ClusterBuilderSubsets::build(
        vector<Cluster> &clusters,
        const RPRecognizedPatternsCollection &rpPatColl)
{
    vector<Candidate> candidates;
    set<Edge> edges;
    populateGraph(candidates, edges, rpPatColl);

    set<set<CandidateId> > subsets; // = list of cliques (= clusters)
    findLargestNonInclusiveCliques(subsets, candidates, edges);
    
    if (verbosity > 2)
    {
        printf("Subsets (%lu):\n", subsets.size());

        unsigned int ss_id = 0;
        for (set<set<CandidateId> >::const_iterator ss_i = subsets.begin(); ss_i != subsets.end(); ss_i++, ss_id++) {
            printf("> subset %u:", ss_id);
            for (set<CandidateId>::const_iterator s_i = ss_i->begin(); s_i != ss_i->end(); s_i++)
                printf(" %u", *s_i);

            printf("\n");
        }
    }

    mergeSubsetsIntoClusters(clusters, subsets, candidates);
}

//----------------------------------------------------------------------------------------------------

void ClusterBuilderSubsets::populateGraph(vector<Candidate> &candidates, set<Edge> &edges,
        const RPRecognizedPatternsCollection &rpPatColl)
{
    // loop over RP pairs
    for (RPRecognizedPatternsCollection::const_iterator rit1 = rpPatColl.begin(); rit1 != rpPatColl.end(); ++rit1)
    {
        for (RPRecognizedPatternsCollection::const_iterator rit2 = rit1; rit2 != rpPatColl.end(); ++rit2)
        {
            if (rit2 == rit1)
                continue;

            if (verbosity > 2)
                printf("%u -- %u\n", rit1->first, rit2->first);

            vector<Candidate> newCandidates = buildCandidatesFrom2RPs(rit1->first, rit2->first, rit1->second, rit2->second);

            for (vector<Candidate>::const_iterator ncIt = newCandidates.begin(); ncIt != newCandidates.end(); ncIt++)
            {
                for (unsigned int c_i = 0; c_i < candidates.size(); c_i++)
                {
                    if (candidates[c_i].dist(*ncIt) < edge_create_th)
                        edges.insert(Edge(c_i, candidates.size()));
                }
                candidates.push_back(*ncIt);
            }
        }
    }

    if (verbosity > 2)
    {
        for (set<Edge>::const_iterator edge_i = edges.begin(); edge_i != edges.end(); edge_i++)
            printf("edge: %u-%u\n", edge_i->first, edge_i->second);
    }
}

//----------------------------------------------------------------------------------------------------

vector<ClusterBuilderSubsets::Candidate> ClusterBuilderSubsets::buildCandidatesFrom2RPs(
    RPId rpID1, RPId rpID2,
    const RPRecognizedPatterns &rpPat1, const RPRecognizedPatterns &rpPat2)
{

    RPInfo rpInfo1 = geometry.getRPInfo(rpID1);
    RPInfo rpInfo2 = geometry.getRPInfo(rpID2);
    //printf("%u: %f %f %f %f %f\n", rpId1, rpInfo1.cz,
            //rpInfo1.u_dx, rpInfo1.u_dy, rpInfo1.v_dx, rpInfo1.v_dy);

    const vector<RPRecognizedPatterns::Line> &pat_1u = rpPat1.uLines;
    const vector<RPRecognizedPatterns::Line> &pat_1v = rpPat1.vLines;
    const vector<RPRecognizedPatterns::Line> &pat_2u = rpPat2.uLines;
    const vector<RPRecognizedPatterns::Line> &pat_2v = rpPat2.vLines;

    vector<Candidate> candidates;

    // loop over track candidates in first RP

    for (unsigned int p1ui = 0; p1ui < pat_1u.size(); p1ui++)
    {
        for (unsigned int p1vi = 0; p1vi < pat_1v.size(); p1vi++)
        {
            for (unsigned int p2ui = 0; p2ui < pat_2u.size(); p2ui++)
            {
                for (unsigned int p2vi = 0; p2vi < pat_2v.size(); p2vi++)
                {
                    if (verbosity > 2)
                        printf("\tp1ui=%u, p1vi=%u | p2ui=%u, p2vi=%u\n", p1ui, p1vi, p2ui, p2vi);

                    Candidate cand = buildCandidate(rpID1, rpID2,
                            pat_1u[p1ui], pat_1v[p1vi],
                            pat_2u[p2ui], pat_2v[p2vi]);

                    if (candidateOk(rpID1, rpID2, cand))
                    {
                        if (verbosity > 2)
                        {
                            double ax = cand.th(0), ax_u = sqrt(cand.V_th(0, 0)),
                                   bx = cand.th(1), bx_u = sqrt(cand.V_th(1, 1)),
                                   ay = cand.th(2), ay_u = sqrt(cand.V_th(2, 2)),
                                   by = cand.th(3), by_u = sqrt(cand.V_th(3, 3));

                            printf("\tcandidate #%lu: ax=%.1f +- %.1f urad, ay=%.1f +- %.1f urad, bx=%.3f +- %.3f mm, by=%.3f +- %.3f mm\n",
                                candidates.size(), ax*1E6, ax_u*1E6, ay*1E6, ay_u*1E6, bx, bx_u, by, by_u);
                        }

                        candidates.push_back(cand);
                    } else {
                        if (verbosity > 2)
                            printf("\t\t outside the rp's edge or too inclined\n");
                    }
                }
            }
        }
    }

    return candidates;
}

//----------------------------------------------------------------------------------------------------

ClusterBuilderSubsets::Candidate ClusterBuilderSubsets::buildCandidate(
        RPId rpID1, RPId rpID2,
        const Line &u1, const Line &v1,
        const Line &u2, const Line &v2)
{

    RPInfo rpInfo1 = geometry.getRPInfo(rpID1),
           rpInfo2 = geometry.getRPInfo(rpID2);

    double u1_p  = u1.b,
           u1_z  = rpInfo1.cz,
           u1_dx = rpInfo1.u_dx,
           u1_dy = rpInfo1.u_dy,

           v1_p  = v1.b,
           v1_z  = rpInfo1.cz,
           v1_dx = rpInfo1.v_dx,
           v1_dy = rpInfo1.v_dy,

           u2_p  = u2.b,
           u2_z  = rpInfo2.cz,
           u2_dx = rpInfo2.u_dx,
           u2_dy = rpInfo2.u_dy,

           v2_p  = v2.b,
           v2_z  = rpInfo2.cz,
           v2_dx = rpInfo2.v_dx,
           v2_dy = rpInfo2.v_dy;

    TVectorD mv(4);
    mv(0) = u1_p;
    mv(1) = v1_p;
    mv(2) = u2_p;
    mv(3) = v2_p;

    TMatrixDSym V_mv(4);
    V_mv(0, 0) = uv_unc * uv_unc;
    V_mv(1, 1) = uv_unc * uv_unc;
    V_mv(2, 2) = uv_unc * uv_unc;
    V_mv(3, 3) = uv_unc * uv_unc;

    double M_data[] = {
        u1_dx*(u1_z - z_h), u1_dx, u1_dy*(u1_z - z_h), u1_dy,
        v1_dx*(v1_z - z_h), v1_dx, v1_dy*(v1_z - z_h), v1_dy,
        u2_dx*(u2_z - z_h), u2_dx, u2_dy*(u2_z - z_h), u2_dy,
        v2_dx*(v2_z - z_h), v2_dx, v2_dy*(v2_z - z_h), v2_dy
    };

    TMatrixD M(4, 4);
    M.SetMatrixArray(M_data);

    TMatrixD Mi(TMatrixD::kInverted, M);
    TMatrixD MiT(TMatrixD::kTransposed, Mi);

    TVectorD th = Mi * mv;
    TMatrixD V_th = Mi * V_mv * MiT;

    TMatrixD V_th_i(TMatrixD::kInverted, V_th);

    Candidate cand;
    cand.th = th;
    cand.V_th = V_th;
    cand.V_th_i = V_th_i;
    cand.rpId1 = rpID1;
    cand.rpId2 = rpID2;

    copy(u1.hits.begin(), u1.hits.end(),
            inserter(cand.hits, cand.hits.begin()));
    copy(v1.hits.begin(), v1.hits.end(),
            inserter(cand.hits, cand.hits.begin()));
    copy(u2.hits.begin(), u2.hits.end(),
            inserter(cand.hits, cand.hits.begin()));
    copy(v2.hits.begin(), v2.hits.end(),
            inserter(cand.hits, cand.hits.begin()));

    return cand;
}

//----------------------------------------------------------------------------------------------------

bool ClusterBuilderSubsets::candidateOk(RPId rpID1, RPId rpID2, const Candidate &candidate)
{
    if (geometry.isPointOutside(rpID1, coordsAtRP(rpID1, candidate)) ||
          geometry.isPointOutside(rpID2, coordsAtRP(rpID2, candidate)))
        return false;

    double ax = candidate.th(0),
           ay = candidate.th(2);
    if (fabs(ax - ax_mean) > ax_cut || fabs(ay - ay_mean) > ay_cut)
        return false;

    return true;
}

//----------------------------------------------------------------------------------------------------

CLHEP::Hep3Vector ClusterBuilderSubsets::coordsAtRP(RPId rpID, const Candidate &candidate) {
    double z = geometry.getRPInfo(rpID).cz;

    double ax = candidate.th(0),
           bx = candidate.th(1),
           ay = candidate.th(2),
           by = candidate.th(3);

    return CLHEP::Hep3Vector(ax*(z-z_h)+bx, ay*(z-z_h)+by, z);
}

//----------------------------------------------------------------------------------------------------

void ClusterBuilderSubsets::findLargestNonInclusiveCliques(set<set<CandidateId> > &subsets,
        const vector<Candidate> &candidates, const set<Edge> &edges)
{
    // TODO: understand

    CandidateId cand_id = 0;
    vector<CandidateId> cand_id_v;
    // two algorithms combined, first look for all cliques composed of vertices (track
    // candidates) with rising ids (no duplicates thanks to this) then remove cliques which
    // can be expanded (are not largest cliques possible) thus having left only the non
    // overlapping largest cliques
    while (true) {
        if (cand_id != candidates.size()) {
            // check if track candidate cand_id is connected to all
            // candidates in cand_id_v, if so attach the cand_id to the cand_id_v
            if (allCandidateIdsConnectedTo(cand_id_v, cand_id, edges)) {
                cand_id_v.push_back(cand_id);
            }
        } else if (cand_id_v.empty()) {
            // end of iteration
            break;
        } else {
            // check if cand_id_v has an interesting clique
            bool interesting = true;
            map<RPId, unsigned int> RP_cand;
            for (vector<CandidateId>::const_iterator cand_id_v_it = cand_id_v.begin();
                    cand_id_v_it != cand_id_v.end(); cand_id_v_it++) {
                const Candidate &cand = candidates[*cand_id_v_it];
                RP_cand[cand.rpId1]++;
                RP_cand[cand.rpId2]++;
            }
            // if less than set number of pots are involved, reject
            if (RP_cand.size() < minRPsPerStation*(minRPsPerStation-1)/2) {
                interesting = false;
            }
            // if any pot appears less than number of RPs needed minus one times
            // in the combinations then reject
            for (map<RPId, unsigned int>::const_iterator RP_cand_it = RP_cand.begin();
                    interesting && RP_cand_it != RP_cand.end(); RP_cand_it++) {
                if (RP_cand_it->second < minRPsPerStation - 1) {
                    interesting = false;
                }
            }
            // check if the clique is actually largest
            // if possible to attach any other candidate, then reject as
            // there must be other clique already found covering this one
            for (cand_id = 0; interesting && cand_id < candidates.size(); cand_id++) {
                // ok not to exclude case where cand_id is in cand_id_v as there are no
                // edge cand_id-cand_id so it is automatically ignored
                if (allCandidateIdsConnectedTo(cand_id_v, cand_id, edges)) {
                    interesting = false;
                }
            }
            // cand_id invalidated

            // if still interested, include in the resulting set
            if (interesting) {
                set<CandidateId> subset(cand_id_v.begin(), cand_id_v.end());
                subsets.insert(subset);
            }

            // resume iteration with last element popped from cand_id_v
            cand_id = cand_id_v.back();
            cand_id_v.pop_back();
        }

        cand_id++;
    }
}

//----------------------------------------------------------------------------------------------------

bool ClusterBuilderSubsets::allCandidateIdsConnectedTo(const vector<CandidateId> &cand_v,
        CandidateId cand_id, const set<Edge> &edges)
{
    for (vector<CandidateId>::const_iterator cand_it = cand_v.begin(); cand_it != cand_v.end(); cand_it++)
    {
        if (edges.find(Edge(*cand_it, cand_id)) == edges.end())
            return false;
    }

    return true;
}

//----------------------------------------------------------------------------------------------------

void ClusterBuilderSubsets::mergeSubsetsIntoClusters(
        vector<Cluster> &clusters,
        const set<set<CandidateId> > &subsets,
        const vector<Candidate> &candidates)
{
    for (set<set<CandidateId> >::const_iterator set_i = subsets.begin(); set_i != subsets.end(); set_i++)
    {
        Cluster cluster;
        for (set<CandidateId>::const_iterator cand_i = set_i->begin(); cand_i != set_i->end(); cand_i++)
        {
            const Candidate &cand = candidates[*cand_i];
            cluster.Add(cand.rpId1, cand.rpId2, cand.th, cand.V_th_i);
            cluster.hits.insert(cand.hits.begin(), cand.hits.end());
        }

        clusters.push_back(cluster);
    }
}

//----------------------------------------------------------------------------------------------------

double ClusterBuilderSubsets::Candidate::dist(const Candidate &other)
{
    TVectorD D = th - other.th;
    TMatrixD D_V_i = V_th + other.V_th;
    D_V_i.SetTol(1e-22); // in case small numbers are contained in the matrix
    D_V_i.Invert();

    double d2 = 0.;
    unsigned int dim = D_V_i.GetNrows();
    for (unsigned int i = 0; i < dim; i++)
        for (unsigned int j = 0; j < dim; j++)
            d2 += D(i) * D_V_i(i, j) * D(j);

    return d2;
}

} // namespace RPStationMultiTrackFinderFitter

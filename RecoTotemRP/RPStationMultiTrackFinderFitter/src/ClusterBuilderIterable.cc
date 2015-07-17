/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Jakub Sawicki (jakub.kuba.sawicki@gmail.com)
 *
 ****************************************************************************/

#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/ClusterBuilderIterable.h"
#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/Cluster.h"

#include "TMatrixD.h"
#include "TVectorD.h"

#include <vector>

namespace RPStationMultiTrackFinderFitter {

ClusterBuilderIterable::ClusterBuilderIterable(const edm::ParameterSet &ps, GeometryHelper &geometry) :
    IClusterBuilder(ps, geometry),
    ax_cut(ps.getParameter<double>("ax_cut")),
    ay_cut(ps.getParameter<double>("ay_cut")),
    z_h(ps.getParameter<double>("z_h")),
    cluster_merge_si_cut(ps.getParameter<double>("cluster_merge_si_cut"))
{
}

//----------------------------------------------------------------------------------------------------

void ClusterBuilderIterable::build(
        vector<Cluster> &clusters,
        const RPRecognizedPatternsCollection &rpPatColl)
{
    // first combination of two RPs?
    bool first_RP_combination = true;

    // loop over RP pairs
    for (RPRecognizedPatternsCollection::const_iterator rit1 = rpPatColl.begin(); rit1 != rpPatColl.end(); ++rit1)
    {
        for (RPRecognizedPatternsCollection::const_iterator rit2 = rit1; rit2 != rpPatColl.end(); ++rit2)
        {
            if (rit2 == rit1)
                continue;

            // updated collection of clusters, after processing this RP pair
            vector<Cluster> newClusters;

            RPId rpId1 = rit1->first;
            RPId rpId2 = rit2->first;

            RPInfo rpInfo1 = geometry.getRPInfo(rpId1);
            RPInfo rpInfo2 = geometry.getRPInfo(rpId2);
            //printf("%u: %f %f %f %f %f\n", rpId1, rpInfo1.cz,
                    //rpInfo1.u_dx, rpInfo1.u_dy, rpInfo1.v_dx, rpInfo1.v_dy);

            if (verbosity > 2)
                printf("%u -- %u\n", rit1->first, rit2->first);

            const vector<RPRecognizedPatterns::Line> &pat_1u = rit1->second.uLines;
            const vector<RPRecognizedPatterns::Line> &pat_1v = rit1->second.vLines;
            const vector<RPRecognizedPatterns::Line> &pat_2u = rit2->second.uLines;
            const vector<RPRecognizedPatterns::Line> &pat_2v = rit2->second.vLines;

            // loop over hit candidates in first RP
            for (unsigned int p1ui = 0; p1ui < pat_1u.size(); p1ui++)
            {
                for (unsigned int p1vi = 0; p1vi < pat_1v.size(); p1vi++)
                {
                    // loop over hit candidates in second RP
                    for (unsigned int p2ui = 0; p2ui < pat_2u.size(); p2ui++)
                    {
                        for (unsigned int p2vi = 0; p2vi < pat_2v.size(); p2vi++)
                        {
                            if (verbosity > 2)
                                printf("\tp1ui=%u, p1vi=%u | p2ui=%u, p2vi=%u\n", p1ui, p1vi, p2ui, p2vi);

                            double u1_p  = pat_1u[p1ui].b,
                                   u1_z  = rpInfo1.cz,
                                   u1_dx = rpInfo1.u_dx,
                                   u1_dy = rpInfo1.u_dy;

                            double v1_p  = pat_1v[p1vi].b,
                                   v1_z  = rpInfo1.cz,
                                   v1_dx = rpInfo1.v_dx,
                                   v1_dy = rpInfo1.v_dy;

                            double u2_p  = pat_2u[p2ui].b,
                                   u2_z  = rpInfo2.cz,
                                   u2_dx = rpInfo2.u_dx,
                                   u2_dy = rpInfo2.u_dy;

                            double v2_p  = pat_2v[p2vi].b,
                                   v2_z  = rpInfo2.cz,
                                   v2_dx = rpInfo2.v_dx,
                                   v2_dy = rpInfo2.v_dy;

                            // check if the two hit candidates aren't outside the sensor's edge
                            double x, y, z;
                            z = u1_z;
                            x = u1_p * u1_dx + v1_p * v1_dx;
                            y = u1_p * u1_dy + v1_p * v1_dy;
                            
                            if (geometry.isPointOutside(rpId1, CLHEP::Hep3Vector(x,y,z)))
                            {
                                if (verbosity > 2)
                                    printf("\t\trejected: p1 outside sensor\n");

                                continue;
                            }

                            z = u2_z;
                            x = u2_p * u2_dx + v2_p * v2_dx;
                            y = u2_p * u2_dy + v2_p * v2_dy;

                            if (geometry.isPointOutside(rpId2, CLHEP::Hep3Vector(x,y,z)))
                            {
                                if (verbosity > 2)
                                    printf("\t\trejected: p2 outside sensor\n");

                                continue;
                            }

                            // determine track-candidate parameters
                            TVectorD mv(4);
                            mv(0) = u1_p;
                            mv(1) = v1_p;
                            mv(2) = u2_p;
                            mv(3) = v2_p;

                            double uv_unc = 20E-3;    // mm
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

                            double ax = th(0), ax_u = sqrt(V_th(0, 0));
                            double bx = th(1), bx_u = sqrt(V_th(1, 1));
                            double ay = th(2), ay_u = sqrt(V_th(2, 2));
                            double by = th(3), by_u = sqrt(V_th(3, 3));

                            if (verbosity > 2)
                                printf("\t\tax=%.1f +- %.1f urad, ay=%.1f +- %.1f urad, bx=%.3f +- %.3f mm, by=%.3f +- %.3f mm\n",
                                    ax*1E6, ax_u*1E6, ay*1E6, ay_u*1E6, bx, bx_u, by, by_u);

                            // TODO: ax_mean ??
                            // reasonability cuts
                            if (fabs(ax) > ax_cut || fabs(ay) > ay_cut)
                                continue;

                            if (first_RP_combination)
                            {
                                // for first RP combination, all track candidates are stored as cluster seeds
                                Cluster c;
                                c.Add(rpId1, rpId2, th, V_th_i);

                                if (verbosity > 2)
                                    printf("\t\t\tnew cluster #%lu\n", newClusters.size());

                                newClusters.push_back(c);

                                if (verbosity > 2)
                                    c.GetCenterTrack().Print("\t\t\t\tcenter: ");
                            } else {
                                // for second and further RP combinations, add the track candidate to
                                // all (already existing) compatible clusters

                                // find indeces of clusters compatible with the track candidate
                                vector<unsigned int> mergeable = vector<unsigned int>();
                                for (unsigned int ci = 0; ci < clusters.size(); ci++)
                                {
                                    if (verbosity > 3)
                                        printf("\t\t\tchecking cluster #%u\n", ci);
    
                                    double dist = clusters[ci].Distance(th, V_th, V_th_i);
      
                                    if (verbosity > 3)
                                    {
                                        const Track &center = clusters[ci].GetCenterTrack();
                                        center.Print("\t\t\t\tcenter: ");
                                        printf("\t\t\t\tdistance: %f\n", dist);
                                    }
      
                                    if (dist < cluster_merge_si_cut)
                                    {
                                        if (verbosity > 2)
                                            printf("\t\t\t**** CREATING cluster #%lu from #%u (last_dist = %.2f)\n", newClusters.size(), ci, dist);
    
                                        mergeable.push_back(ci);
                                    }
                                }

                                // add the track candidate to the previously selected clusters
                                for (unsigned int ci_i = 0; ci_i < mergeable.size(); ci_i++)
                                {
                                    Cluster newCluster;
                                    clusters[mergeable[ci_i]].Copy(newCluster);
                                    newCluster.Add(rpId1, rpId2, th, V_th_i);
                                    newClusters.push_back(newCluster);
                                }
                            }
                        }
                    }
                }
            }

            // replace cluster collection with the updated one
            // this implicitly discards clusters with no track candidate from this RP pair
            clusters.clear();
            clusters.swap(newClusters);

            first_RP_combination = false;
        }
    }

    if (verbosity > 2)
        printf(">> Clusters: %lu\n", clusters.size());
}


} // namespace RPStationMultiTrackFinderFitter

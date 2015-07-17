/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Jakub Sawicki (jakub.sawicki@cern.ch, jakub.kuba.sawicki@gmail.com)
 *
 ****************************************************************************/

#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/OverlapsRemoval.h"

namespace RPStationMultiTrackFinderFitter {

OverlapsRemoval::OverlapsRemoval(const edm::ParameterSet &ps) :
    verbosity(ps.getUntrackedParameter<unsigned int>("verbosity"))
{
}

//----------------------------------------------------------------------------------------------------

OverlapsRemoval::~OverlapsRemoval()
{
}

//----------------------------------------------------------------------------------------------------

void OverlapsRemoval::run(
        vector<Cluster> &out,
        const vector<Cluster> &in,
        const vector<unsigned int> &hitCount,
        const vector<vector<unsigned int> > &cluHits)
{
    if (verbosity > 2)
        printf("\n>> OverlapsRemoval::run\n");     

    this->hitCountOriginal = &hitCount;
    this->cluHits = &cluHits;

    calculatePossibleHitCount();

    vector<bool> chosen = vector<bool>(in.size(), false);
    v_chosen = vector<vector<bool> >();

    goDown(chosen, 0);

    // if more than 1 possibilities then reject all
    for (unsigned int c_i = 0; c_i < chosen.size(); c_i++)
    {
        if (v_chosen.size() == 1 && v_chosen[0][c_i])
        {
            out.push_back(in[c_i]);

            if (verbosity > 2)
                printf("%u true\n", c_i);
        } else {
            if (verbosity > 2)
                printf("%u false\n", c_i);
        }
    }
}

//----------------------------------------------------------------------------------------------------

void OverlapsRemoval::calculatePossibleHitCount()
{
    hitCountPossible = vector<unsigned int>(*hitCountOriginal);

    if (verbosity > 2)
        printf("Possible hits count:");

    for (unsigned int plane_i = 0; plane_i < hitCountOriginal->size(); plane_i++)
    {
        // how many cluster hits are associated with each U/V pattern
        vector<unsigned int> clu_per_hit = vector<unsigned int>((*hitCountOriginal)[plane_i]+1, 0);

        // iterate over all candidates to see which candidate belongs to which pattern
        for (unsigned int clu_i = 0; clu_i < cluHits->size(); clu_i++)
            clu_per_hit[(*cluHits)[clu_i][plane_i]]++;
        clu_per_hit.pop_back();

        // count patterns which have at least one candidate associated with them
        unsigned int maxHits = 0;
        for (unsigned int hit_i = 0; hit_i < clu_per_hit.size(); hit_i++)
        {
            if (clu_per_hit[hit_i] > 0)
                maxHits++;
        }

        hitCountPossible[plane_i] = maxHits;

        if (verbosity > 2)
            printf(" %u", maxHits);
    }

    if (verbosity > 2)
        printf("\n");
}

//----------------------------------------------------------------------------------------------------

void OverlapsRemoval::goDown(vector<bool> &chosen, unsigned int pos)
{
    if (pos >= cluHits->size())
        return;

    // too many possible combinations so reject all
    if (v_chosen.size() > 1)
        return;

    bool feasible, complete;

    chosen[pos] = true;
    feasibilityCheck(&feasible, &complete, chosen);
    if (verbosity > 2)
    {
        printf("case: ");
        for (unsigned int c_i = 0; c_i < chosen.size(); c_i++)
            printf("%u ",(chosen[c_i]) ? 1 : 0);
        printf("feas: %u complete: %u\n", feasible ? 1 : 0, complete ? 1 : 0);
    }

    if (complete)
    {
        // don't go down, as more tracks won't fit the set of chosen ones anyway
        vector<bool> copied_chosen;
        copied_chosen.reserve(chosen.size());
        copy(chosen.begin(), chosen.end(), back_inserter(copied_chosen));
        v_chosen.push_back(copied_chosen);
    } else {
        // it is acceptable but not enough tracks has been associated so go down
        if (feasible)
          goDown(chosen, pos+1);
    }

    chosen[pos] = false;

    // go down without check
    // * won't complete with no track chosen
    // * if not feasible path then it would've been rejected sooner
    goDown(chosen, pos+1);

    // leaving the method the chosen[pos] must be set to false to be consistent
}

//----------------------------------------------------------------------------------------------------

void OverlapsRemoval::feasibilityCheck(bool *feasible, bool *complete,
        const vector<bool> &chosen)
{
    bool _feasible = true, _complete = true;

    // get the maximum number of hits per plane
    // it should equal the number of real tracks
    // we must've had lots of bad luck for overlaps
    // to occur in all planes simultaneously
    unsigned int maxHits = 0;
    for (unsigned int i = 0; i < hitCountPossible.size(); i++)
        maxHits = max(maxHits, hitCountPossible[i]);

    // analyse each plane by itself
    for (unsigned int plane_i = 0; plane_i < hitCountPossible.size() && _feasible; plane_i++)
    {
        // if some hits are missing it probably is due to overlaps
        unsigned int maxOverlaps = maxHits - hitCountPossible[plane_i];

        // how many tracks per hit we have
        vector<unsigned int> deps = vector<unsigned int>((*hitCountOriginal)[plane_i]+1, 0);
        for (unsigned int c_i = 0; c_i < cluHits->size(); c_i++) if (chosen[c_i])
            deps[(*cluHits)[c_i][plane_i]]++;

        unsigned int overlapping = 0, missing_associations = 0;
        for (unsigned int hit_i = 0; hit_i < (*hitCountOriginal)[plane_i]; hit_i++)
        {
            // how many hits (patterns) doesn't have a candidate associated
            if (deps[hit_i] < 1)
                missing_associations++;

            // count how many overlaps we have in this plane
            if (deps[hit_i] > 1)
                overlapping += deps[hit_i] - 1;
        }

        // if there are more missing than the difference between the original and possible
        // number of hit counts than we are not complete yet
        if (missing_associations > (*hitCountOriginal)[plane_i] - hitCountPossible[plane_i])
            _complete = false;

        // too many overlaps, reject
        if (overlapping > maxOverlaps)
        {
            _feasible = false;
            _complete = false;
        }
    }

    if (feasible != 0)
      *feasible = _feasible;

    if (complete != 0)
      *complete = _complete;
}

} // namespace RPStationMultiTrackFinderFitter

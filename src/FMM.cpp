/** 
 * @file
 * @brief Implementation of Fast Multipole Method for force calculation.
 */

#include <memory>
#include <cassert>
#include <utility>
#include "FMM.hpp"

namespace sim {

/**
 * @class FMM
 * @brief Fast Multipole Method implementation.
 */

FMM::FMM(int p,
        double theta,
        int maxPerCell,
        int maxPairwiseLimit,
        Kernels* kernels,
        Force const* forceLaw):
    Interaction(forceLaw),
    mP(p),
    mTheta(theta),
    mMaxPerCell(maxPerCell),
    mMaxPairwiseLimit(maxPairwiseLimit),
    mKernels(kernels) {};

/**
 * @brief Multipole acceptance criteria between two octree nodes
 *
 * @return Whether it is allowed to multipole expand between the two nodes.
 */
bool FMM::MAC(Octree const* node1, Octree const* node2) const {
    Vec dr = node1->coc - node2->coc;
    double dist = dr.GetNorm();
    // TODO: does this really work well for general rect boxes?
    double sum_lengths = node1->GetMaxLength() + node2->GetMaxLength();
    return dist * mTheta >= sum_lengths;
}

/**
 * @brief Dual tree traversal for M2L expansions.
 *
 * This implementation is similar to that described in Dehnen 2002, except for a
 * simpler boundary condition. With a max limit on leaf particle number at tree
 * construction time, we are not allowed to further subdivide when we encounter
 * leaf nodes.
 *
 * The method is destructive, and will modify grid and nodes as it traverses the
 * tree to save multipole coefficients and partial forces
 *
 * @param[in, out] node1 First node in interaction pair
 * @param[in, out] node2 Second node in interaction pair
 * @param[in, out] grid Grid to which we save F coefficients and partial forces.
 */
void FMM::Interact(Octree* node1, Octree* node2, Grid& grid) const {
    assert(node1);
    assert(node2);
    // a size of 5e4 could easily overflow the pairwise terms, and cause all
    // interactoins to go through the brute-force method... Very funny.
    int64_t n1 = node1->mSouls.size();
    int64_t n2 = node2->mSouls.size();
    int64_t pairwise_terms = n1 * n2;
    if (node1 == node2) {
        // self-interaction, pair number different.
        // n1 >= 1 for all created nodes.
        assert(n1 >= 1);
        pairwise_terms = (n1 * (n1 - 1)) / 2;
    }

    if ((node1->IsLeaf() && node2->IsLeaf())
            || pairwise_terms <= mMaxPairwiseLimit)
    {
        // pairwise_terms = 0 is fine, this means there should be no effect from
        // one of the nodes, and this (node1, node2) pair may be safely ignored.
        // lower than pairwise limit or both leafs (irreducible) so carry out
        // direct mutual interaction
        if (node1 == node2) {
            // self-interacting nodes
            for (int i = 0; i < n1; i++) {
                Particle& p1 = grid[node1->mSouls[i]];
                for (int j = i + 1; j < n1; j++) {
                    Particle& p2 = grid[node1->mSouls[j]];
                    const Vec& force = GetForce(p1, p2);
                    const double pot = GetPot(p1, p2);

                    p1.accel += force / p1.GetMass();
                    p2.accel -= force / p2.GetMass();

                    p1.pot += pot;
                    p2.pot += pot;
                }
            } 
        } else {
            // inter-node direct interaction
            for (int s1 : node1->mSouls) {
                Particle& p1 = grid[s1];
                for (int s2 : node2->mSouls) {
                    Particle& p2 = grid[s2];  
                    const Vec& force = GetForce(p1, p2);
                    const double pot = GetPot(p1, p2);

                    p1.accel += force / p1.GetMass();
                    p2.accel -= force / p2.GetMass();

                    p1.pot += pot;
                    p2.pot += pot;
                }
            }
        }
    } else if (MAC(node1, node2)) {
        // number of particles exceeded, but well-separated
        // should cover parent-neightbour child-non-neighbour nodes
        // We leave F tensors as lazy labels for now
        mKernels->M2L(node1, node2);
        mKernels->M2L(node2, node1);
    } else {
        // "interaction cannot be performed" per Dehnen 2002
        // can and must recurse downwards for all pairs,
        // including pair of the same cell for box self-interactions
        if (node1 == node2) {
            // recurse all pairs of children
            for (int i = 0; i < node1->mBoxes; i++) {
                Octree* child_node1 = node1->GetChild(i);
                if (!child_node1) continue;
                for (int j = i; j < node1->mBoxes; j++) {
                    Octree* child_node2 = node1->GetChild(j);
                    if (!child_node2) continue;
                    Interact(child_node1, child_node2, grid);
                }
            }
        } else {
            // IsLeaf():    node1   node2
            // 1            false   false   larger node -> node1;
            // 2            false   true    keep node1
            // 3            true    false   swap to node2
            // 4            true    true    don't care, have directly interacted
            if (node1->IsLeaf()) {
                // cases 3 & 4
                std::swap(node1, node2);
            } else if (!node2->IsLeaf()
                    && node1->GetMaxLength() < node2->GetMaxLength()) {
                //case 1
                std::swap(node1, node2);
            } // else case 2, do nothing.
            for (int i = 0; i < node1->mBoxes; i++) {
                Octree* child_node = node1->GetChild(i);
                if (!child_node) continue;
                Interact(child_node, node2, grid);
            }
        }
    }
}

/**
 * @brief Evaluate acceleration by using L2L push-down, then L2P.
 * 
 * The method is destructive as both node and grid will be modified. This is
 * similar to algorithm 2 (evaluate gravity) of Dehnen 2002.
 *
 * @param[in, out] node Current node of DFS
 * @param[in, out] grid Grid on which new \f$F_n^m\f$, accel and pot are saved.
 */
void FMM::EvaluateAccel(Octree* node, Grid& grid) const {
    if (node->IsLeaf()) {
        for (int soul : node->mSouls) {
            Particle& par = grid[soul];
            mKernels->L2P(node, par);
        }
    } else {
        for (int i = 0; i < node->mBoxes; i++) {
            Octree* child_node = node->GetChild(i);
            if (!child_node) continue;
            // push-down
            mKernels->L2L(node, child_node);
            EvaluateAccel(child_node, grid);
        }
    }
}

/**
 * @brief Calculate acceleration and potential of all particles in grid g1
 *
 * @param[in] g1 Initial grid for which acceleration and potential is calculated
 * @return Grid with acceleration and potential information
 */
Grid FMM::Calculate(const Grid& g1) const {
    Grid grid = g1;
    std::unique_ptr<Octree> root = Octree::BuildTree(g1, mMaxPerCell, mP);
    mKernels->CalculateM(root.get());
    Interact(root.get(), root.get(), grid);
    EvaluateAccel(root.get(), grid); 
    return grid;
}

}

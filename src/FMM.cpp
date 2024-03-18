#include "FMM.hpp"

namespace sim {

FMM::FMM(int p,
        double theta,
        int maxPerCell,
        int maxPairwiseLimit,
        std::unique_ptr<const Kernels> kernels,
        std::unique_ptr<const Force> forceLaw):
    Interaction(std::move(forceLaw)),
    mP(p),
    mTheta(theta),
    mMaxPerCell(maxPerCell),
    mMaxPairwiseLimit(maxPairwiseLimit),
    mKernels(std::move(kernels)) {};

bool FMM::MAC(const Octree* node1, const Octree* node2) const {
    Vec dr = node1->com - node2->com;
    double dist = dr.GetNorm();
    // TODO: does this really work well for general rect boxes?
    double sum_lengths = node1->GetMaxLength() + node2->GetMaxLength();
    return dist * mTheta > sum_lengths;
}

// dual tree traversal similar to Dehnen 2002
// with simplified boundary cases: leaf boxes have a maximum size,
// and is not allowed to further subdivide.
void FMM::Interact(Octree* node1, Octree* node2, Grid& grid) const {
    int n1 = node1->mSouls.size();
    int n2 = node2->mSouls.size();
    int pairwise_terms = n1 * n2;

    // impose order
    if (node1->GetMaxLength() < node2->GetMaxLength()) {
        std::swap(node1, node2);
    }

    if (node1 == node2 && node1->IsLeaf()) {
        // self-interacting leaf nodes
        for (int i = 0; i < n1; i++) {
            Particle& p1 = grid[node1->mSouls[i]];
            for (int j = i + 1; j < n1; j++) {
                Particle& p2 = grid[node1->mSouls[j]];
                const Vec& force = GetForce(p1, p2);
                p1.accel += force / p1.GetMass();
                p2.accel -= force / p2.GetMass();
            }
        } 
    } else if (pairwise_terms < mMaxPairwiseLimit ||
              (node1->IsLeaf() && node2->IsLeaf())) {
        // small interaction count, or different leaf nodes 
        // so we can't reduce further
        for (int s1 : node1->mSouls) {
            Particle& p1 = grid[s1];
            for (int s2 : node2->mSouls) {
                Particle& p2 = grid[s2];  
                const Vec& force = GetForce(p1, p2);
                p1.accel += force / p1.GetMass();
                p2.accel -= force / p2.GetMass();
            }
        }
    } else if (MAC(node1, node2)) {
        // number of particles exceeded, but well-separated
        // should cover parent-neightbour child-non-neighbour nodes
        // We leave F tensors as lazy labels for now
        mKernels->M2L(node1, node2);
        // TODO: symmetry argument here?
        mKernels->M2L(node2, node1);
    } else {
        // "interaction cannot be performed" per Dehnen 2002
        // can and must recurse downwards (neighbours)
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
            // r_max(1) >= r_max(2) as guaranteed by std::swap
            // break down and recurse larger box 1.
            for (int i = 0; i < node1->mBoxes; i++) {
                Octree* child_node = node1->GetChild(i);
                if (!child_node) continue;
                Interact(child_node, node2, grid);
            }
        }
    }
}

// algorithm 2 (evaluate gravity) as per Dehnen 2002
void FMM::EvaluateAccel(Octree* node, Grid& grid) const {
    if (node->IsLeaf()) {
        for (int soul : node->mSouls) {
            Particle& par = grid[soul];
            par = mKernels->L2P(node, par);
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

Grid FMM::Calculate(const Grid& g1) const {
    Grid grid = g1;
    std::unique_ptr<Octree> root = Octree::BuildTree(g1, mMaxPerCell, mP);
    mKernels->CalculateM(root.get());
    Interact(root.get(), root.get(), grid);
    EvaluateAccel(root.get(), grid); 
    return grid;
}

}

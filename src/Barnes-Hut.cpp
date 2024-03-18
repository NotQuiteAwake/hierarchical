#include <cassert>
#include "Barnes-Hut.hpp"

namespace sim {

bool BarnesHut::IsFarAway(Octree const* node, const Particle& par) const {
    const Vec dr = par.pos - node->com;
    const double dist = dr.GetNorm();
    double l = node->GetMaxLength();
    return l < dist * mTheta; // for better double precision
}

Particle& BarnesHut::EvaluateAccel(
        int soul,
        Particle& par,
        Octree const* node
) const {
    // this node does not exist.
    if (!node) { return par; }
    //DEBUG
    //std::cout << node->com << std::endl;

    if (node->IsLeaf()) {
        assert(node->mSouls.size() == 1);
        if (node->mSouls[0] == soul) {
            // this node is the one with the par
            return par;
        }
        else {
            // different single particle node
            const Particle& par2 = node->GetParticle(node->mSouls[0]);
            par.accel += GetForce(par, par2) / par.GetMass();
            return par;
        }
    } else if (IsFarAway(node, par)) {
        // employ multiple M2P kernel
        par = mKernels->M2P(node, par);
    } else {
        // not far enough away, descend into children
        for (int i = 0; i < node->mBoxes; i++) {
            const Octree* child_node = node->GetChild(i);
            par = EvaluateAccel(soul, par, child_node);
        }
    }
    return par;
}

BarnesHut::BarnesHut(int p,
                     double theta,
                     std::unique_ptr<const Kernels> kernels,
                     std::unique_ptr<const Force> forceLaw) :
    Interaction(std::move(forceLaw)),
    mP(p),
    mTheta(theta),
    mKernels(std::move(kernels)) {};

int BarnesHut::GetP() const {
    return mP;
}

double BarnesHut::GetTheta() const {
    return mTheta;
}

Grid BarnesHut::Calculate(const Grid& g1) const {
    Grid grid = g1;
    std::unique_ptr<Octree> root = Octree::BuildTree(g1, 1, mP);
    mKernels->CalculateM(root.get());

    for (int i = 0; i < g1.GetSize(); i++) {
        grid[i] = EvaluateAccel(i, grid[i], root.get());
    }
    return grid;
}

}

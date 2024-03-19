#include <cassert>
#include "Kernels.hpp"

namespace sim {

Kernels::Kernels(int p): mP(p) {};

void Kernels::M2L(Octree const* source, Octree* sink) const {
    assert(source && sink);
    
    sink->F += M2X(source, sink->com);
}

void Kernels::M2P(Octree const* source, Particle& sinkPar) const {
    assert(source);

    ComplexMatrix psi = M2X(source, sinkPar.pos);
    AddAccel(sinkPar, psi);
}

void Kernels::L2L(Octree const* parent, Octree* child) const {
    assert(parent && child);
    assert(parent == child->GetParent());

    const Vec& sp = child->com;
    child->F += L2X(parent, sp);
}

// TODO: store the psi matrix onto the particles
// particle MUST be contained in leaf's octant.
// assertion is NOT made due to performance cost.
void Kernels::L2P(Octree const* leaf, Particle& containedPar) const {
    assert(leaf);
    assert(leaf->IsLeaf());

    ComplexMatrix psi(mP, mP);
    const Vec& x = containedPar.pos;

    psi = L2X(leaf, x);
    AddAccel(containedPar, psi);
}

void Kernels::CalculateM(Octree* node) const {
    assert(node); // in this implementation must not trip this

    if (node->IsLeaf()) {
        P2M(node);
    } else {
        for (int i = 0; i < node->mBoxes; i++) {
            Octree* child_node = node->GetChild(i);
            if (!child_node) continue;

            CalculateM(child_node);
            M2M(child_node, node);
        }
    }
}

}

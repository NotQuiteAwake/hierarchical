/**
 * @file
 * @brief Definition of the general Kernels interface.
 */

#include <cassert>
#include "Kernels.hpp"

namespace sim {

/**
 * @class Kernels
 *
 * @brief Definition of the general Kernels interface.
 *
 * Non-const reference or pointer parameters in method signature are usually
 * subject to direct modification, as they are the the "target" of the
 * operation.
 */

/**
 * @brief Initialises Kernels
 *
 * @param[in] p order of multipole expansion
 */
Kernels::Kernels(int p): mP(p) {};

void Kernels::M2L(Octree const* source, Octree* sink) {
    assert(source && sink);
    
    sink->F += M2X(source, sink->coc);
}

void Kernels::M2P(Octree const* source, Particle& sinkPar) {
    assert(source);

    ComplexMatrix psi = M2X(source, sinkPar.pos);
    AddAccel(sinkPar, psi);
}

void Kernels::L2L(Octree const* parent, Octree* child) {
    assert(parent && child);
    assert(parent == child->GetParent());

    const Vec& sp = child->coc;
    child->F += L2X(parent, sp);
}

/**
 * @brief Apply L2P kernel to particle contained in leaf node.
 *
 * There is a reason we don't enumerate particles on leaf directly. We can make
 * leaf non-const, but Octree stores mGrid as const to prevent side effects. We
 * can make that non-const as well, but when building the Octree sometimes we
 * have only available to us the original grid (We shouldn't need to first
 * create the result grid just to build the octree anyway.) Creating a copy of
 * the grid also won't work, since modifications to it won't go back to the
 * output result grid anyway. All things considered this awkward signature is
 * the best compromise and makes its effects explicit. (What a headache!)
 *
 * @param[in] leaf Leaf node
 * @param[in, out] containedPar Particle to apply L2P to
 */

void Kernels::L2P(Octree const* leaf, Particle& containedPar) {
    assert(leaf);
    assert(leaf->IsLeaf());
    assert(leaf->GetOctant().Within(containedPar.pos));

    ComplexMatrix psi(mP, mP);
    const Vec& x = containedPar.pos;

    psi = L2X(leaf, x);
    AddAccel(containedPar, psi);
}

/**
 * @brief Calculate M coefficients on a fully-built octree.
 *
 * @param[in, out] node Current position in DFS. Start from root.
 */
void Kernels::CalculateM(Octree* node) {
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

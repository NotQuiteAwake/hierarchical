/**
 * @file
 * @brief Octree implementation
 */

#include "Octree.hpp"
#include <cassert>

namespace sim {

/**
 * @class Octree
 * @brief Octree (node) implementation.
 */

/**
 * @brief Initialise an Octree node
 *
 * @param[in] parent pointer to parent node. Take nullptr for root.
 * @param[in] grid reference to grid the Octree is representing.
 * @param[in] new_oct Octant of space this node represents.
 * @param[in] maxParticles Max number of particles in the leaf node.
 * @param[in] p Order of muiltipole expansion.
 */
Octree::Octree(Octree* parent,
               const Grid& grid,
               const Octant& new_oct,
               int maxParticles,
               int p) :
    mMaxParticles(maxParticles),
    mP(p),
    mParent(parent),
    mGrid(grid),
    mOctant(new_oct),
    M(ComplexMatrix(p, p)),
    F(ComplexMatrix(p, p)) {};

/**
 * @brief Set provided node as child at (index), taking ownership.
 *
 * The child is saved to a std::unique_ptr and so should not and cannot be
 * deleted elsewhere. Each node has ownership over all their direct children,
 * not the parent. Outside access of nodes should be with raw pointers as they
 * don't take ownership. This way the ownership structure is clear and the tree
 * can be recursively destructed upon destruction of root.
 *
 * @param[in] index index (octantNumber) of the child node @param[in] octree
 * Pointer to node to be set as child
 */
void Octree::SetChild(int index, Octree* octree) {
    mLeafCount += bool(octree) - bool(mChildren[index]);
    mChildren[index].reset(octree);
}

/**
 * @brief Get raw pointer to child at (index)
 */
Octree* Octree::GetChild(int index) {
    return mChildren[index].get();
}

/**
 * @brief Get raw pointer to constant child at (index)
 */
Octree const* Octree::GetChild(int index) const {
    return mChildren[index].get();
}

/**
 * @brief Get raw pointer to parent.
 */
Octree* Octree::GetParent() {
    return mParent;
}

/**
 * @brief Get raw pointer to constant parent.
 */
Octree const* Octree::GetParent() const {
    return mParent;
}

/**
 * @brief Get length of longest side of the octant represented by this node.
 */
double Octree::GetMaxLength() const {
    return mOctant.GetMaxLength();
}

/**
 * @brief Test if current node is a leaf.
 */
bool Octree::IsLeaf() const {
    return mLeafCount == 0;
}

/**
 * @brief Create an empty child node with its octant number as its index
 *
 * @param[in] octantNumber octant number of the child.
 */
void Octree::GenChildNode(int octantNumber) {
    assert(!mChildren[octantNumber]);
    const Octant new_oct = mOctant.GetOctant(octantNumber);
    // pointer handed to unique_ptr
    SetChild(octantNumber, new Octree(this, mGrid, new_oct, mMaxParticles, mP));
}

/**
 * @brief Split the leaf node and distribute particles to sub-octants
 *
 * A child node is only created if there will be non-zero number of particles
 * that end up in it.
 */
void Octree::Split() {
    assert(mSouls.size()); // why else do you split?
    assert(IsLeaf()); // must not repeat splitting
    
    for (int i = 0; i < mBoxes; i++) {
        // empty octants are not instantiated
        if (mOctantSouls[i].size()) {
            GenChildNode(i);

            for (int soul : mOctantSouls[i]) {
                mChildren[i]->AddParticle(soul);
            }
            
            mOctantSouls[i].clear();    // information pushed down
        }
    }
}

/**
 * @brief Get the 3D box/octant that this node represents.
 */
Octant Octree::GetOctant() const {
    return mOctant;
}

/**
 * @brief Add a particle to this node.
 *
 * At a leaf node if limit maxParticles is exceeded then the node is
 * automatically split. This should always be called at the root node, for as
 * the particle is pushed down the relevant COCs are also updated.
 *
 * @param[in] soul Index of the particle in grid
 */
void Octree::AddParticle(int soul) {
    mSouls.push_back(soul);
    const Particle& par = GetParticle(soul);
    const Vec& par_pos = par.pos;
    assert(mOctant.Within(par_pos));
    const double par_charge = par.GetCharge();
    const int octant_num = mOctant.GetOctantNumber(par_pos);

    coc = (coc * charge + par_pos * par_charge) / (charge + par_charge);
    charge += par_charge;
    if (IsLeaf()) {
        mOctantSouls[octant_num].push_back(soul);
        // soul is pushed down via mOctantSouls
        if (mSouls.size() > mMaxParticles) { Split(); } 
    } else { 
        // soul is pushed down directly
        if (!GetChild(octant_num)) {
            // this child node does not exist 
            GenChildNode(octant_num); 
        }
        GetChild(octant_num)->AddParticle(soul);
    }
}

/**
 * @brief Overwrite the octant that the current node represnets.
 *
 * This should not be called outside of sim::IO, where this is used to load back
 * the octants dumped to files.
 */
void Octree::SetOctant(const Octant& octant) {
    assert(!mOctant.IsInitialised());
    mOctant = octant;
} 

/**
 * @brief Build octree from a grid.
 *
 * @param[in] grid Provides list of particles to construct the octree from
 * @param[in] maxParticles Max number of particles at leaf node
 * @param[in] p order of multipole expansion
 */
std::unique_ptr<Octree> Octree::BuildTree(const Grid& grid,
                                          int maxParticles,
                                          int p) {
    std::unique_ptr<Octree> root(
            new Octree(nullptr, grid, grid.GetOctant(), maxParticles, p)
            );
    for (int i = 0; i < grid.GetSize(); i++) {
        root->AddParticle(i);
    }
    return root;
}

/**
 * @brief Get a copy of a particle from the grid by its soul number.
 *
 */
Particle Octree::GetParticle(int soul) const {
    return mGrid[soul];
}

int Octree::GetP() const {
    return mP;
}

int Octree::GetMaxParticles() const {
    return mMaxParticles;
}

}

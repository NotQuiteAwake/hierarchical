#include "Octree.hpp"
#include <cassert>

namespace sim {

Octree::Octree(Octree* parent,
               const Grid& grid,
               const Octant& new_oct,
               int maxParticles,
               int p) :
    mMaxParticles(maxParticles),
    mP(p),
    mParent(parent),
    mGrid(grid),
    octant(new_oct),
    M(ComplexMatrix(p, p)),
    F(ComplexMatrix(p, p)) {};

void Octree::SetChild(int index, Octree* octree) {
    mLeafCount += bool(octree) - bool(mChildren[index]);
    mChildren[index].reset(octree);
}

Octree* Octree::GetChild(int index) {
    return mChildren[index].get();
}

Octree const* Octree::GetChild(int index) const {
    return mChildren[index].get();
}

Octree* Octree::GetParent() {
    return mParent;
}

Octree const* Octree::GetParent() const {
    return mParent;
}

double Octree::GetMaxLength() const {
    return octant.GetMaxLength();
}

bool Octree::IsLeaf() const {
    return mLeafCount == 0;
}

void Octree::GenChildNode(int octantNumber) {
    assert(!mChildren[octantNumber]);
    const Octant new_oct = octant.GetOctant(octantNumber);
    // pointer handed to unique_ptr
    SetChild(octantNumber, new Octree(this, mGrid, new_oct, mMaxParticles, mP));
}

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

void Octree::AddParticle(int soul) {
    mSouls.push_back(soul);
    const Particle& par = GetParticle(soul);
    const Vec& par_pos = par.pos;
    assert(octant.Within(par_pos));
    const double par_mass = par.GetMass();
    const int octant_num = octant.GetOctantNumber(par_pos);

    com = (com * mass + par_pos * par_mass) / (mass + par_mass);
    mass += par_mass;
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

std::unique_ptr<Octree> Octree::BuildTree(const Grid& grid,
                                          int maxParticles,
                                          int p) {
    std::unique_ptr<Octree> root(
            new Octree(nullptr, grid, grid.GetLimits(), maxParticles, p)
            );
    for (int i = 0; i < grid.GetSize(); i++) {
        root->AddParticle(i);
    }
    return root;
}

Particle Octree::GetParticle(int soul) const {
    return mGrid[soul];
}

int Octree::GetMaxParticles() const {
    return mMaxParticles;
}

}

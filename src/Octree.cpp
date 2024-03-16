#include "Octree.hpp"
#include <cassert>
#include <iostream>

namespace sim {

Octree::Octree(Octree* parent, const Grid& grid, const Octant& new_oct,
        int maxParticles)
    : mParent(parent), mGrid(grid), octant(new_oct),
    mMaxParticles(maxParticles) {};

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

double Octree::GetMaxLength() const {
    return octant.GetMaxLength();
}

bool Octree::IsLeaf() const {
    return mLeafCount == 0;
}

void Octree::GenChildNode(int octantNumber) {
    assert(!mChildren[octantNumber]);
    const Octant new_oct = octant.GetOctant(octantNumber);
    SetChild(octantNumber, new Octree(this, mGrid, new_oct, mMaxParticles));
}

void Octree::Split() {
    assert(mSouls.size()); // why else do you split?
    assert(IsLeaf()); // must not repeat splitting
    
    for (int i = 0; i < mBoxes; i++) {
        // empty octants are not instantiated
        if (mOctantSouls[i].size()) {
            GenChildNode(i);
            Octree& child = *mChildren[i];

            for (int soul : mOctantSouls[i]) {
                child.AddParticle(soul);
            }
            
            mOctantSouls[i].clear();    // information pushed down
        }
    }
}

void Octree::AddParticle(int soul) {
    mSouls.push_back(soul);
    const Vec& pos = GetParticle(soul).pos;
    const int octant_num = octant.GetOctantNumber(pos);
    if (IsLeaf()) {
        mOctantSouls[octant_num].push_back(soul);
        // soul is pushed down via mOctantSouls
        if (mSouls.size() > mMaxParticles) { Split(); } 
    } else { // need to push down soul directly
        // this child node does not exist 
        if (!GetChild(octant_num)) {
            GenChildNode(octant_num); 
        }
        GetChild(octant_num)->AddParticle(soul);
    }
    
}

void Octree::BuildAsRoot() {
    assert(mParent == nullptr); // must be root
    for (int i = 0; i < mGrid.GetSize(); i++) {
        AddParticle(i);
    }
}

Particle Octree::GetParticle(int soul) const {
    return mGrid[soul];
}

int Octree::GetMaxParticles() const {
    return mMaxParticles;
}

}

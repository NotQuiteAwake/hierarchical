#include "Octree.hpp"
#include <cassert>
#include <algorithm>

namespace sim {

Octree::Octree(Octree* parent, const Grid& grid, const Octant& new_oct,
        int maxParticles)
    : mParent(parent), mGrid(grid), octant(new_oct),
    mMaxParticles(maxParticles) {};

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
    return mIsLeaf;
}

void Octree::Split() {
    assert(mSouls.size()); // why else do you split?
    assert(mIsLeaf); // must not repeat splitting
    
    for (int i = 0; i < mBoxes; i++) {
        // empty octants are not instantiated
        if (mOctantSouls[i].size()) {
            const Octant new_oct = octant.GetOctant(i);
            
            mChildren[i].reset(new Octree(this, mGrid, new_oct, mMaxParticles));
            Octree& child = *mChildren[i];

            for (int soul : mOctantSouls[i]) {
                const Particle& particle = mGrid[soul];
                child.AddParticle(particle, soul);
            }
            
            mOctantSouls[i].clear();    // information pushed down
        }
    }
    
    mIsLeaf = false; 
}

void Octree::AddParticle(const Particle& par, int soul) {
    mSouls.push_back(soul);
    const Vec& pos = par.GetPos();
    const int octantNumber = octant.GetOctantNumber(pos);
    mOctantSouls[octantNumber].push_back(soul);
    
    if (mSouls.size() > mMaxParticles) { Split(); }
}

Particle Octree::GetParticle(int soul) const {
    return mGrid[soul];
}

}

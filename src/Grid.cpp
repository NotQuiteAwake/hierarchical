#include <cassert>
#include "Grid.hpp"

namespace sim {

Grid::Grid(Octant maxLim): mMaxLim(maxLim), mOctant(Octant()) {}

Grid::Grid(int size, Octant maxLim): 
    mMaxLim(maxLim),
    mOctant(Octant()) { Reserve(size); }
    
int Grid::GetSize() const {
    return mParticles.size();
}

const Octant Grid::GetLimits() const {
    return mMaxLim;
}

const Octant Grid::GetOctant() const {
    return mOctant;
}

Particle& Grid::operator[](int index) {
    return mParticles[index];
}

Particle Grid::operator[](int index) const {
    return mParticles[index];
}

void Grid::AddParticle(const Particle& par) {
    const Vec& pos = par.pos;
    mOctant.Relax(pos);
    // otherwise impose no limit on dynamic box size
    if (mMaxLim.IsInitialised()) {
        if (mMaxLim.Within(pos)) {
            mParticles.push_back(par);
        }
    } else {
        mParticles.push_back(par);
    }
}

void Grid::SetOctant(const Octant& octant) {
    assert(!mOctant.IsInitialised());
    mOctant = octant;
}

void Grid::Reserve(int size) {
    mParticles.reserve(size);
}

}

#include "Grid.hpp"

namespace sim {

Grid::Grid(Octant octant) : mOctant(octant) {}
    
int Grid::GetSize() const {
    return mParticles.size();
}

const Octant Grid::GetLimits() const {
    return mOctant;
}

Particle& Grid::operator[](int index) {
    return mParticles[index];
}

Particle Grid::operator[](int index) const {
    return mParticles[index];
}

const std::vector<Particle>& Grid::GetParticles() const {
    return mParticles;
}

void Grid::AddParticle(const Particle& par) {
    mParticles.push_back(par);
}

}

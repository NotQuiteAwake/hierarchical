#include <cassert>
#include "Grid.hpp"

namespace sim {

Grid::Grid(Octant maxLim): mMaxLim(maxLim), mOctant(Octant()) {}

Grid::Grid(int size, Octant maxLim): 
    mMaxLim(maxLim),
    mOctant(Octant()) { Reserve(size); }
    
Grid::Grid(const std::vector<Particle>& pars, Octant maxLim):
    mMaxLim(maxLim),
    mOctant(Octant()) {
        for (const Particle& par : pars) {
            // necessary, or else mOctant.Relax is not called.
            AddParticle(par);
        }
    }

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

void Grid::AddParticles(const std::vector<Particle>& par_list) {
    for (const Particle& par : par_list) {
        AddParticle(par);
    }
}

double Grid::GetPE() const {
    double PE = 0;
    for (const Particle& par : mParticles) {
        PE += par.GetPE();
    }
    return PE / 2; // we've double counted
}

double Grid::GetKE() const {
    double KE = 0;
    for (const Particle& par : mParticles) {
        KE += par.GetKE();
    }
    return KE;
}

double Grid::GetE() const {
    return GetKE() + GetPE();
}

Vec Grid::GetCOM() const {
    Vec com({0, 0, 0});
    double mass = 0; 
    for (const Particle& par : mParticles) {
        const double par_mass = par.GetMass();
        com += par.pos * par_mass;
        mass += par_mass;
    }

    assert(mass > 0);
    return com / mass;
}

Vec Grid::GetL(const Vec centre) const {
    Vec L({0, 0, 0});
    for (const Particle& par : mParticles) {
        L += CrossProduct(par.pos - centre, par.vel) * par.GetMass();
    }
    return L;
}

Vec Grid::GetP() const {
    Vec P({0, 0, 0});
    for (const Particle& par : mParticles) {
        P += par.GetP();
    }
    return P;
}

void Grid::SetOctant(const Octant& octant) {
    assert(!mOctant.IsInitialised());
    mOctant = octant;
}

void Grid::Reserve(int size) {
    mParticles.reserve(size);
}

}

/**
 * @file
 * @brief Implement Grid to describe simulation space and particles within.
 */

#include <cassert>
#include "Grid.hpp"

namespace sim {

/**
 * @class Grid
 * @brief Data structure to describe simulation space and particles within.
 *
 * If maxLim is not initialised, grid expands its simulation space dynamically
 * as particles are appended. if maxLim is initialised however, space can grow
 * only up to maxLim, and particles lying outside are rejected on appending.
 */

/**
 * @brief Initialise setting maxLim
 *
 * @param[in] maxLim Maximum limits for the simulation space
 */
Grid::Grid(Octant maxLim): mMaxLim(maxLim), mOctant(Octant()) {}

/**
 * @brief Initialise grid, preallocating space for (size) particles
 *
 * @param[in] size Expected number of particles
 * @param[in] maxLim Maximum limits for simulation space
 */
Grid::Grid(int size, Octant maxLim): 
    mMaxLim(maxLim),
    mOctant(Octant()) { Reserve(size); }
    
/**
 * @brief Initialise grid with a list of particles
 *
 * @param[in] pars List of particles to add to grid
 * @param[in] maxLim maximum limits for simulation space
 */
Grid::Grid(const std::vector<Particle>& pars, Octant maxLim):
    mMaxLim(maxLim),
    mOctant(Octant()) {
        for (const Particle& par : pars) {
            // necessary, or else mOctant.Relax is not called.
            AddParticle(par);
        }
    }

/**
 * @brief Get number of particles in grid
 *
 * @return Number of particles in grid.
 */
int Grid::GetSize() const {
    return mParticles.size();
}

/**
 * @brief Get Maximum limits imposed on simulation space
 *
 * @return An Octant object. If uninitialised, grid does not have a size limit.
 */
const Octant Grid::GetLimits() const {
    return mMaxLim;
}

/**
 * @brief Get current simulation space limits.
 *
 * @return Current simulation space limits.
 */
const Octant Grid::GetOctant() const {
    return mOctant;
}

/**
 * @brief Get particle by index (soul) in grid
 *
 * @return Particle at provided index
 */
Particle& Grid::operator[](int index) {
    return mParticles[index];
}

/**
 * @brief Get Particle by index (soul) in grid
 *
 * @return particle at provided index
 */
Particle Grid::operator[](int index) const {
    return mParticles[index];
}

/**
 * @brief Append particle to grid, resizing simulation space if needed.
 *
 * If maxLim is initialised for the grid, and the particle is outside that
 * limit, then it is merely ignored and not added to the grid.
 *
 * @param[in] par Particle to be added
 */
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

/**
 * @brief Append a list of particles to grid
 *
 * @param[in] par_list List of particles to be added.
 */
void Grid::AddParticles(const std::vector<Particle>& par_list) {
    for (const Particle& par : par_list) {
        AddParticle(par);
    }
}

/**
 * @brief Get total potential energy in grid
 *
 * This assumes the interaction is symmetric and pairwise.
 *
 * @return Total potential energy.
 */
double Grid::GetPE() const {
    double PE = 0;
    for (const Particle& par : mParticles) {
        PE += par.GetPE();
    }
    return PE / 2; // we've double counted
}

/**
 * @brief Get total kinetic energy in grid
 *
 * @return Total kinetic energy
 */
double Grid::GetKE() const {
    double KE = 0;
    for (const Particle& par : mParticles) {
        KE += par.GetKE();
    }
    return KE;
}

/**
 * @brief Get total energy in grid
 *
 * @return Total energy.
 */
double Grid::GetE() const {
    return GetKE() + GetPE();
}

/**
 * @brief Get centre of mass of grid
 *
 * @return Centre of mass of all particles in grid.
 */
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

/**
 * @brief Get angular momentum of grid
 *
 * @param[in] centre Centre for L calculation. Defaults to the origin.
 * @return Angular momentum about set centre
 */
Vec Grid::GetL(const Vec centre) const {
    Vec L({0, 0, 0});
    for (const Particle& par : mParticles) {
        L += CrossProduct(par.pos - centre, par.vel) * par.GetMass();
    }
    return L;
}

/**
 * @brief Get total momentum of grid
 *
 * @return Sum of momentum of all particles
 */
Vec Grid::GetP() const {
    Vec P({0, 0, 0});
    for (const Particle& par : mParticles) {
        P += par.GetP();
    }
    return P;
}

/**
 * @brief Force set grid dimensions
 *
 * This is only possible when the grid dimensions are uninitialised, that is
 * when there are no particles in grid.
 *
 * @param[in] Grid dimensions
 */
void Grid::SetOctant(const Octant& octant) {
    assert(GetSize() == 0);
    assert(!mOctant.IsInitialised());
    mOctant = octant;
}

/**
 * @brief Reserve space for (size) particles
 *
 * @param[in] size Number of particles to reserve space for.
 */
void Grid::Reserve(int size) {
    mParticles.reserve(size);
}

}

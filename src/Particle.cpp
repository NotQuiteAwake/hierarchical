/** 
 * @file
 * @brief Data structure representing a particle.
 */

#include <cassert>
#include "Particle.hpp"

namespace sim {

/** 
 * @class Particle
 * @brief Data structure representing a particle.
 */

// C++ constructs other members nontheless. Default values specified in header.
Particle::Particle(double mass, double charge):
    mMass(mass),
    mCharge(charge) {
        assert(mass > 0);
    }

double Particle::GetMass() const {
    return mMass;
}

double Particle::GetCharge() const {
    return mCharge;
}

/**
 * @brief Get potential energy on the particle
 *
 * Note the risk of double counting, since for pairwise forces potential is
 * really a value shared between a pair, while we store this whole value in a
 * particle.
 *
 */
double Particle::GetPE() const {
    return pot;
}

/**
 * @brief Get kinetic energy on the particle
 */
double Particle::GetKE() const {
    double speed = vel.GetNorm();
    return mMass * speed * speed / 2;
}

/**
 * @brief Get momentum on the particle
 */
Vec Particle::GetP() const {
    return vel * mMass;
}

}

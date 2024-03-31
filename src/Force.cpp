/**
 * @file
 * @brief Virtual class for general pairwise force and its implementations
 */

#include <cassert>
#include "Force.hpp"

namespace sim {

/**
 * @class Force
 * @brief Virtual class for a general pairwise force.
 */

/**
 * @brief Get force of p2 on p1.
 *
 * This checks for particle overlap. Derived classes implement instead ForceLaw,
 * which carries out no such checks.
 *
 * @return If two particles do not overlap return force. Otherwise, null vector.
 */
Vec Force::GetForce(const Particle& p1, const Particle& p2) const {
    if (CheckDistinct(p1, p2)) {
        return ForceLaw(p1, p2);
    } else {
        // in debug version CheckDistinct trips asserts.
        // in production as last ditch effort don't return nan's...
        return Vec();
    }
}

/**
 * @brief Get potential between two particles.
 *
 * This checks for particle overlap. Derived classes implement instead PotLaw,
 * which does not carry out such checks.
 *
 * @return If two particles do not overlap return potential. Otherwise, 0.
 */
double Force::GetPot(const Particle& p1, const Particle& p2) const {
    if (CheckDistinct(p1, p2)) {
        return PotLaw(p1, p2);
    } else {
        return 0;
    }
}

/**
 * @brief Check and assert that two particles are distinct and not overlapping
 *
 */
bool Force::CheckDistinct(const Particle& p1, const Particle& p2) const {
    bool same_par_addr = (&p1 == &p2);
    bool same_pos = (p1.pos == p2.pos);
    assert(!same_par_addr);
    assert(!same_pos);
    return !(same_par_addr || same_pos);
}

/**
 * @class DummyForce
 * @brief Dummy implementation of Force that does nothing.
 */
Vec DummyForce::ForceLaw(const Particle& p1, const Particle& p2) const {
    return Vec();
}

double DummyForce::PotLaw(const Particle& p1, const Particle& p2) const {
    return 0;
}

/**
 * @class InvSqForce
 * @brief Implementation of Force that provides the inverse-square law.
 * 
 */

/**
 * @brief Initialiser for InvSqForce
 *
 * @param[in] G Coupling constant. In our convention G < 0 for gravity.
 */
InvSqForce::InvSqForce(double G): mG(G) {};

Vec InvSqForce::ForceLaw(const Particle& p1, const Particle& p2) const {
    Vec dr = p2.pos - p1.pos;
    double dist = dr.GetNorm();
    return -dr * mG * p1.GetCharge() * p2.GetCharge() / (dist * dist * dist);
} 

double InvSqForce::PotLaw(const Particle& p1, const Particle& p2) const {
    Vec dr = p2.pos - p1.pos;
    double dist = dr.GetNorm();
    return mG * p1.GetCharge() * p2.GetCharge() / dist;
}

}

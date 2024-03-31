/**
 * @file Interaction
 * @brief Define the general Interaction interface.
 */

#include "Interaction.hpp"

namespace sim {

/**
 * @class Interaction
 * @brief Class defining the general Interaction interface.
 */

/**
 * @brief Get force due to p2 on p1.
 *
 * @return force.
 */
Vec Interaction::GetForce(const Particle& p1, const Particle& p2) const {
    return mForceLaw->GetForce(p1, p2);
}

double Interaction::GetPot(const Particle& p1, const Particle& p2) const {
    return mForceLaw->GetPot(p1, p2);
}

Interaction::Interaction(Force const* forceLaw): mForceLaw(forceLaw) {}

}

#include "Interaction.hpp"

namespace sim {

Vec Interaction::GetForce(const Particle& p1, const Particle& p2) const {
    return mForceLaw -> GetForce(p1, p2);
}

Interaction::Interaction(std::unique_ptr<const Force> forceLaw):
    mForceLaw(std::move(forceLaw)) {}

}

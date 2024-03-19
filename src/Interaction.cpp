#include "Interaction.hpp"

namespace sim {

Vec Interaction::GetForce(const Particle& p1, const Particle& p2) const {
    return mForceLaw->GetForce(p1, p2);
}

Interaction::Interaction(Force const* forceLaw): mForceLaw(forceLaw) {}

}

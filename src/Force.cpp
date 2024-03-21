#include <cassert>
#include "Force.hpp"

namespace sim {

Vec Force::GetForce(const Particle& p1, const Particle& p2) const {
    if (CheckDistinct(p1, p2)) {
        return ForceLaw(p1, p2);
    } else {
        // in debug version CheckDistinct trips asserts.
        // in production as last ditch effort don't return nan's...
        return Vec();
    }
}

bool Force::CheckDistinct(const Particle& p1, const Particle& p2) const {
    bool same_par_addr = (&p1 == &p2);
    bool same_pos = (p1.pos == p2.pos);
    assert(!same_par_addr);
    assert(!same_pos);
    return !(same_par_addr || same_pos);
}

Gravity::Gravity(double G): mG(G) {};

Vec Gravity::ForceLaw(const Particle& p1, const Particle& p2) const {
    Vec dr = p2.pos - p1.pos;
    double dist = dr.GetNorm();

    // no 'self-interaction'!
    assert(&p1 != &p2);
    assert(dist);

    return  dr * mG * p1.GetMass() * p2.GetMass() / (dist * dist * dist);
} 

Vec DummyForce::ForceLaw(const Particle& p1, const Particle& p2) const {
    return Vec();
}

}

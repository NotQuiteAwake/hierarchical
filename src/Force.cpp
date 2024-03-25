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

double Force::GetPot(const Particle& p1, const Particle& p2) const {
    if (CheckDistinct(p1, p2)) {
        return PotLaw(p1, p2);
    } else {
        return 0;
    }
}

bool Force::CheckDistinct(const Particle& p1, const Particle& p2) const {
    bool same_par_addr = (&p1 == &p2);
    bool same_pos = (p1.pos == p2.pos);
    assert(!same_par_addr);
    assert(!same_pos);
    return !(same_par_addr || same_pos);
}

Vec DummyForce::ForceLaw(const Particle& p1, const Particle& p2) const {
    return Vec();
}

double DummyForce::PotLaw(const Particle& p1, const Particle& p2) const {
    return 0;
}

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

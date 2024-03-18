#include "Particle.hpp"

namespace sim {

// C++ constructs other members nontheless.
Particle::Particle(double mass, double charge):
    mMass(mass),
    mCharge(charge) {}

double Particle::GetMass() const {
    return mMass;
}

double Particle::GetCharge() const {
    return mCharge;
}

}

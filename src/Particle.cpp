#include "Particle.hpp"

namespace sim {

// C++ constructs members nontheless.
Particle::Particle() {}

Particle::Particle(double mass, double charge) {
    mMass = mass;
    mCharge = charge;
}

double Particle::GetMass() const {
    return mMass;
}

double Particle::GetCharge() const {
    return mCharge;
}

Vec Particle::GetPos() const {
    return mPos;
}

Vec Particle::GetVel() const {
    return mVel;
}

Vec Particle::GetAccel() const {
    return mAccel;
}

void Particle::SetPos(const Vec& pos) {
    mPos = pos;
}

void Particle::SetVel(const Vec& vel) {
    mVel = vel;
}

void Particle::SetAccel(const Vec& accel) {
    mAccel = accel;
}

}

#include "Force.hpp"
#include <cmath>

namespace sim {

Gravity::Gravity(double G): mG(G) {};

Vec Gravity::GetForce(const Particle& p1, const Particle& p2) const {
    Vec dr = p2.pos - p1.pos;
    double dist = dr.GetNorm();
    return  dr * mG * p1.GetMass() * p2.GetMass() / pow(dist, 3);
} 

Vec DummyForce::GetForce(const Particle& p1, const Particle& p2) const {
    return Vec();
}

}

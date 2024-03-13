#include "Force.hpp"
#include <cmath>

namespace sim {

Vec Gravity::GetForce(const Particle& p1, const Particle& p2) const {
    Vec dr = p2.GetPos() - p1.GetPos();
    double dist = dr.CalculateNorm();
    // assumes G = 1
    return  dr * p1.GetMass() * p2.GetMass() / pow(dist, 3);
} 

}

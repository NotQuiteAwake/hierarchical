#include "Brute.hpp"

namespace sim {

Brute::Brute(std::unique_ptr<const Force> forceLaw):
    Interaction(std::move(forceLaw)) {};

Grid Brute::Calculate(const Grid& g1) const {
    int size = g1.GetSize();
    Grid g2 = g1;
    for (int i = 0; i < size; i++) {
        const Particle& pi = g1[i];
        for (int j = i + 1; j < size; j++) {
            const Particle& pj = g1[j];
            Vec force = GetForce(pi, pj); 
            // update acceleration pair-wise
            g2[i].accel = g2[i].accel + force / g2[i].GetMass();
            g2[j].accel = g2[j].accel - force / g2[j].GetMass();
        }
    }

    return g2;
}

}

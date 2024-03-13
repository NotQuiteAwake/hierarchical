#include "Brute.hpp"

namespace sim {

Grid Brute::Calculate(const Grid& g1) const {
    int size = g1.GetSize();
    Grid g2 = g1;
    for (int i = 0; i < size; i++) {
        const Particle& pi = g1[i];
        for (int j = i + 1; j < size; j++) {
            const Particle& pj = g1[j];
            Vec force = GetForce(pi, pj); 
            // update acceleration pair-wise
            g2[i].SetAccel(g2[i].GetAccel() + force / g2[i].GetMass());
            g2[j].SetAccel(g2[j].GetAccel() - force / g2[j].GetMass());
        }
    }

    return g2;
}

}

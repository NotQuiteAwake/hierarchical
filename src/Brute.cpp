/**
 * @file
 * @brief Implementation of brute-force method for force calculation
 */

#include "Brute.hpp"

namespace sim {

/**
 * @class Brute
 * @brief Brute-force implementation of force calculation for whole grid
 */

Brute::Brute(Force const* forceLaw):
    Interaction(forceLaw) {};

/**
 * @brief Calculate potential and acceleration for each charge in grid.
 *
 * @param[in] g1 input grid
 * @return Grid containing particles whose accel and potential are updated.
 */
Grid Brute::Calculate(const Grid& g1) const {
    int size = g1.GetSize();
    Grid g2 = g1;
    for (int i = 0; i < size; i++) {
        const Particle& pi = g1[i];
        for (int j = i + 1; j < size; j++) {
            const Particle& pj = g1[j];
            Vec force = GetForce(pi, pj); 
            double pot = GetPot(pi, pj);
            // update acceleration pair-wise
            g2[i].accel += force / g2[i].GetMass();
            g2[j].accel -= force / g2[j].GetMass();

            g2[i].pot += pot;
            g2[j].pot += pot;
        }
    }

    return g2;
}

}

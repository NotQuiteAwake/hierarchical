/**
 * @file
 * @brief Define general integrator and implement Euler and leapfrog.
 */

#include "Integrator.hpp"
#include "Vec.hpp"

namespace sim {

/**
 * @class Integrator
 * @brief Defines the general integrator interface.
 */

/**
 * @brief Initialise integrator
 *
 * @param[in] step Step size of integrator.
 */
Integrator::Integrator(double step): mStep(step) {}

double Integrator::GetStep() const { return mStep; }

/**
 * @brief Evolve a particle forward by internal step size
 *
 * @param[in] p1 Particle with acceleration precalculated
 * @return Particle after an internal step size.
 */
Particle Integrator::Evolve(const Particle& p1) const {
    return Evolve(p1, mStep);
}

/**
 * @brief Evolve full grid forward in time by step.
 *
 * Derived classes must implement Evolve(Particle, double).
 *
 * @param[in] g1 Grid to be evolved in time
 * @param[in] step step size
 * @return Grid evolved in time
 */
Grid Integrator::Evolve(const Grid& g1, double step) const {
    Grid g2 = Grid(g1.GetLimits());
    g2.Reserve(g1.GetSize()); // reduce allocation on the fly
    for (int i = 0; i < g1.GetSize(); i++) {
        g2.AddParticle(Evolve(g1[i], step));
    }
    return g2;
}

/**
 * @brief Evolve full grid forward in time by internal step size
 *
 * @param[in] g1 Grid to be evolved in time
 * @return Grid evolved in time
 */
Grid Integrator::Evolve(const Grid& g1) const {
    return Evolve(g1, mStep);
}

/**
 * @class Euler
 * @brief Euler integrator implementing the Integrator interface
 */

Euler::Euler(double step): Integrator(step) {}

Particle Euler::Evolve(const Particle& p1, const double step) const {
    Vec new_vel = p1.vel + p1.accel * step;
    Vec new_pos = p1.pos + p1.vel * step;

    Particle p2(p1.GetMass(), p1.GetCharge());
    p2.vel = new_vel;
    p2.pos = new_pos;

    return p2;
}

/**
 * @class LeapFrog
 * @brief Leap frog implementation of the Integrator interface.
 */

LeapFrog::LeapFrog(double step): Integrator(step) {}

Particle LeapFrog::Evolve(const Particle& p1, const double step) const {
    Vec new_vel = p1.vel + p1.accel * step;
    Vec new_pos = p1.pos + new_vel * step;

    Particle p2(p1.GetMass(), p1.GetCharge());
    p2.vel = new_vel;
    p2.pos = new_pos;

    return p2;
}

// TODO: implement RK4?

}

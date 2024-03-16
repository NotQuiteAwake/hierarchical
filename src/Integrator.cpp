#include "Integrator.hpp"
#include "Vec.hpp"

namespace sim {

Integrator::Integrator(double step): mStep(step) {}

double Integrator::GetStep() const { return mStep; }

Particle Integrator::Evolve(const Particle& p1) const {
    return Evolve(p1, mStep);
}

Grid Integrator::Evolve(const Grid& g1, double step) const {
    Grid g2 = Grid(g1.GetLimits());
    g2.Reserve(g1.GetSize()); // reduce allocation on the fly
    for (int i = 0; i < g1.GetSize(); i++) {
        g2.AddParticle(Evolve(g1[i], step));
    }
    return g2;
}

Grid Integrator::Evolve(const Grid& g1) const {
    return Evolve(g1, mStep);
}

Euler::Euler(double step): Integrator(step) {}

Particle Euler::Evolve(const Particle& p1, const double step) const {
    Vec new_vel = p1.vel + p1.accel * step;
    Vec new_pos = p1.pos + p1.vel * step;

    Particle p2(p1.GetMass(), p1.GetCharge());
    p2.vel = new_vel;
    p2.pos = new_pos;

    return p2;
}

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

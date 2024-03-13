#include "Integrator.hpp"
#include "Vec.hpp"

namespace sim {

Integrator::Integrator(double step): mStep(step) {}

double Integrator::GetStep() const { return mStep; }

Particle Integrator::Euler(const Particle& p1) const {
    Vec new_vel = p1.GetVel() + p1.GetAccel() * mStep;
    Vec new_pos = p1.GetPos() + p1.GetVel() * mStep;

    Particle p2;
    p2.SetVel(new_vel);
    p2.SetPos(new_pos);

    return p2;
}

Grid Integrator::Euler(const Grid& g1) const {
    Grid g2 = Grid(g1.GetLimits());
    for (int i = 0; i < g1.GetSize(); i++) {
        g2.AddParticle(Euler(g1[i]));
    }
    return g2;
}

Particle Integrator::LeapFrog(const Particle& p1) const {
    Vec new_vel = p1.GetVel() + p1.GetAccel() * mStep;
    Vec new_pos = p1.GetPos() + new_vel * mStep;

    Particle p2;
    p2.SetVel(new_vel);
    p2.SetPos(new_pos);

    return p2;
}

Grid Integrator::LeapFrog(const Grid& g1) const {
    Grid g2 = Grid(g1.GetLimits());
    for (int i = 0; i < g1.GetSize(); i++) {
        g2.AddParticle(LeapFrog(g1[i]));
    }
    return g2;
}

// TODO: implement RK4?

}

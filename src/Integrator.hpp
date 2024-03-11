#ifndef INTEGRATORHEADERDEF
#define INTEGRATORHEADERDEF

#include "Grid.hpp"
#include "Particle.hpp"

namespace sim {

class Integrator {
    private:
        const double mStep;
         
    public:
        Integrator(double step);
        double GetStep() const;

        Particle Euler(const Particle& p1) const;
        Grid Euler(const Grid& g1) const;

        Particle LeapFrog(const Particle& p1) const;
        Grid LeapFrog(const Grid& g1) const;

        Particle RK4(const Particle& p1) const;
        Grid RK4(const Grid& g1) const;
};

}

#endif

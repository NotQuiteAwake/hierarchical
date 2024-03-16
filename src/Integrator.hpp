#ifndef INTEGRATORHEADERDEF
#define INTEGRATORHEADERDEF

#include "Grid.hpp"
#include "Particle.hpp"

namespace sim {

class Integrator {
    protected:
        const double mStep;
         
    public:
        Integrator(double step);
        double GetStep() const;

        virtual Particle Evolve(const Particle& p1, const double step) const = 0;
        Particle Evolve(const Particle& p1) const;
        Grid Evolve(const Grid& g1, const double step) const;
        Grid Evolve(const Grid& g1) const;
};

class Euler : public Integrator {
    public:
        using Integrator::Evolve;
        Euler(double step);
        Particle Evolve(const Particle& p1, const double step) const override;
};

class LeapFrog : public Integrator {
    public:
        using Integrator::Evolve;
        LeapFrog(double step);
        Particle Evolve(const Particle& p1, const double step) const override;
};

}

#endif

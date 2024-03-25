#ifndef FORCEHEADERDEF
#define FORCEHEADERDEF

#include "Vec.hpp"
#include "Particle.hpp"

namespace sim {

// this is much better than passing in a *funptr...
class Force {
    public:
        // else leads to undefined behaviour
        virtual ~Force() = default;
        // We can only do pairwise interactions...
        // Further require that forces obey Newton III.
        // employ template method pattern.
        Vec GetForce(const Particle& p1, const Particle& p2) const;
        double GetPot(const Particle& p1, const Particle& p2) const;

    protected:
        virtual Vec ForceLaw(const Particle& p1, const Particle& p2) const = 0;
        virtual double PotLaw(const Particle& p1, const Particle& p2) const = 0;

    private:
        bool CheckDistinct(const Particle& p1, const Particle& p2) const;
};

class DummyForce : public Force {
    public:
        Vec ForceLaw(const Particle& p1, const Particle& p2) const override;
        double PotLaw(const Particle& p1, const Particle& p2) const override;
};

class InvSqForce : public Force {
    protected:
        double mG;

    public:
        // defaults to attractive coupling (gravity)
        InvSqForce(double G = -1);
        Vec ForceLaw(const Particle& p1, const Particle& p2) const override;
        double PotLaw(const Particle& p1, const Particle& p2) const override;
};

}

#endif

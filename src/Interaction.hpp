#ifndef INTERACTIONHEADERDEF
#define INTERACTIONHEADERDEF

#include "Force.hpp"
#include "Vec.hpp"
#include "Particle.hpp"
#include "Grid.hpp"

namespace sim {

class Interaction {
    private:
        Force const* mForceLaw;

    protected:
        Vec GetForce(const Particle& p1, const Particle& p2) const;
        double GetPot(const Particle& p1, const Particle& p2) const;

    public:
        // allows for polymorphism.
        Interaction(Force const* forceLaw);
        virtual ~Interaction() = default;
        virtual Grid Calculate(const Grid& g1) const = 0;
};

}

#endif

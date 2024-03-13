#ifndef FORCEHEADERDEF
#define FORCEHEADERDEF

#include "Vec.hpp"
#include "Particle.hpp"

namespace sim {

// this is much better than passing in a *funptr...
class Force {
    public:
        // We can only do pairwise interactions...
        // Further require that forces obey Newton III.
        virtual Vec GetForce(const Particle& p1, const Particle& p2) const = 0;
};

class Gravity : public Force {
    public:
        Vec GetForce(const Particle& p1, const Particle& p2) const override;
};

}

#endif

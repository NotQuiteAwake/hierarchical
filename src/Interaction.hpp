#ifndef INTERACTIONHEADERDEF
#define INTERACTIONHEADERDEF

#include <memory>
#include "Force.hpp"
#include "Vec.hpp"
#include "Particle.hpp"
#include "Grid.hpp"

namespace sim {

class Interaction {
    private:
        const std::shared_ptr<const Force> mForceLaw;

    protected:
        Vec GetForce(const Particle& p1, const Particle& p2) const;

    public:
        Interaction(const std::shared_ptr<const Force> forceLaw); // allows for polymorphism
        virtual Grid Calculate(const Grid& g1) const = 0;
};

}

#endif

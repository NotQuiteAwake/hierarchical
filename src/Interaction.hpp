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
        const std::unique_ptr<const Force> mForceLaw;

    protected:
        Vec GetForce(const Particle& p1, const Particle& p2) const;

    public:
        // allows for polymorphism.
        // default to dummy force as some calculations are force specific
        // and ignore the forceLaw
        Interaction(const std::unique_ptr<const Force> forceLaw
                = std::unique_ptr<const Force>(new DummyForce())); 
        virtual ~Interaction() = default;
        virtual Grid Calculate(const Grid& g1) const = 0;
};

}

#endif

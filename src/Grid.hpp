#ifndef GRIDHEADERDEF
#define GRIDHEADERDEF

#include <vector>
#include "Particle.hpp"
#include "Octant.hpp"

namespace sim {

class Grid {
    private:
        const Octant mOctant;
        std::vector<Particle> mParticles;

    public:
        Grid(Octant octant);

        int GetSize() const;
        const Octant GetLimits() const;

        Particle& operator[](int index);
        Particle operator[](int index) const;

        const std::vector<Particle>& GetParticles() const;
        void AddParticle(const Particle& par);
};

}

#endif


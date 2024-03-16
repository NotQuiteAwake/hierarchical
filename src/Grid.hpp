#ifndef GRIDHEADERDEF
#define GRIDHEADERDEF

#include <vector>
#include "Particle.hpp"
#include "Octant.hpp"

namespace sim {

class Grid {
    private:
        Octant mOctant;
        // note, that vector will allocate its memory on the heap. GREAT!
        // THANK GOD!
        std::vector<Particle> mParticles;

    public:
        Grid(Octant octant);

        int GetSize() const;
        const Octant GetLimits() const;

        Particle& operator[](int index);
        Particle operator[](int index) const;

        void AddParticle(const Particle& par);
        void Reserve(int size);
};

}

#endif


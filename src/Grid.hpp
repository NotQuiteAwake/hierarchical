#ifndef GRIDHEADERDEF
#define GRIDHEADERDEF

#include <vector>
#include "Particle.hpp"
#include "Octant.hpp"

namespace sim {

class Grid {
    private:
        Octant mMaxLim;
        Octant mOctant;

        std::vector<Particle> mParticles; // heap allocated

    public:
        Grid(Octant maxLim = Octant());
        Grid(int size, Octant maxLim = Octant());
        Grid(const std::vector<Particle>& pars, Octant maxLim = Octant());

        int GetSize() const;
        const Octant GetLimits() const;
        const Octant GetOctant() const;

        Particle& operator[](int index);
        Particle operator[](int index) const;

        void AddParticle(const Particle& par);
        void AddParticles(const std::vector<Particle>& par_list);
        double GetPE() const;
        double GetKE() const;
        double GetE() const;
        Vec GetCOM() const;
        Vec GetL(const Vec centre = Vec({0, 0, 0})) const;
        Vec GetP() const;

        void SetOctant(const Octant& octant);
        void Reserve(int size);
};

}

#endif


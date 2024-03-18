#ifndef KERNELSHEADERDEF
#define KERNELSHEADERDEF

#include "Octree.hpp"
#include "Matrix.hpp"

namespace sim {

class Kernels {
    protected:
        int mP;

    public:
        Kernels(int p);
        virtual ~Kernels() = default;

        // virtual particles... like photons in QED...
        virtual Particle AddAccel(const Particle& par,
                const ComplexMatrix& psi) const = 0;

        virtual void P2M(Octree* leaf) const = 0;
        virtual void M2M(Octree const* child, Octree* parent) const = 0;

        virtual ComplexMatrix M2X(Octree const* source, const Vec& s) const = 0;
        void M2L(Octree const* source, Octree* sink) const;
        Particle M2P(Octree const* source, const Particle& sinkPar) const;

        virtual ComplexMatrix L2X(Octree const* previous,
                const Vec& sp) const = 0;
        void L2L(Octree const* parent, Octree* child) const;
        Particle L2P(Octree* leaf, const Particle& containedPar) const;
        
        void CalculateM(Octree* node) const;
};

}

#endif

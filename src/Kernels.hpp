#ifndef KERNELSHEADERDEF
#define KERNELSHEADERDEF

#include "Octree.hpp"
#include "Matrix.hpp"

namespace sim {

class Kernels {
    protected:
        const int mP;

    public:
        Kernels(int p);
        virtual ~Kernels() = default;

        virtual void AddAccel(Particle& par,
                const ComplexMatrix& psi) const = 0;

        virtual void P2M(Octree* leaf) = 0;
        virtual void M2M(Octree const* child, Octree* parent) = 0;

        virtual ComplexMatrix M2X(Octree const* source, const Vec& s) = 0;
        void M2L(Octree const* source, Octree* sink);
        void M2P(Octree const* source, Particle& sinkPar);

        virtual ComplexMatrix L2X(Octree const* previous, const Vec& sp) = 0;
        void L2L(Octree const* parent, Octree* child);
        void L2P(Octree const* leaf, Particle& containedPar);
        
        void CalculateM(Octree* node);
};

}

#endif

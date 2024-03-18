#ifndef INVSQKERNELSHEADERDEF
#define INVSQKERNELSHEADERDEF

#include <complex>
#include "Kernels.hpp"

namespace sim {

class InvSqKernels : public Kernels {
    public:
        typedef std::complex<double> cdouble;

        InvSqKernels(int p);

        double Prefactor(int n, int m) const;
        // surface spherical harmonics via Boost
        cdouble Y(const Vec& v, int n, int m) const;
        // solid spherical harmonics as per Dehnen 2014
        cdouble Gamma(const Vec& v, int n, int m) const;
        cdouble Theta(const Vec& v, int n, int m) const;

        Particle AddAccel(const Particle& par,
                const ComplexMatrix& F) const override;
        void P2M(Octree* leaf) const override;
        void M2M(Octree const* child, Octree* parent) const override;
        ComplexMatrix M2X(Octree const* source, const Vec& s) const override;
        ComplexMatrix L2X(Octree const* previous, const Vec& sp) const override;
};

}
#endif

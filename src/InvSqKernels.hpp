#ifndef INVSQKERNELSHEADERDEF
#define INVSQKERNELSHEADERDEF

#include <complex>
#include "Kernels.hpp"

namespace sim {

class InvSqKernels : public Kernels {
    private:
        const int useBoostLimit = -1;
        const double mG;
        ComplexMatrix mTempMatrix;

    public:
        typedef std::complex<double> cdouble;

        // defaults to attractive force (gravity)
        InvSqKernels(int p, double G = -1);

        double Prefactor(int n, int m) const;
        // surface spherical harmonics via Boost
        cdouble Y(const Vec& v, int n, int m) const;
        // solid spherical harmonics as per Dehnen 2014 (boost)
        cdouble GammaBoost(const Vec& v, int n, int m) const;
        cdouble ThetaBoost(const Vec& v, int n, int m) const;
        cdouble Gamma(const Vec& v, int n, int m) const;
        cdouble Theta(const Vec& v, int n, int m) const;
        void Gamma(const Vec& v, int n);
        void Theta(const Vec& v, int n);
        
        ComplexMatrix GammaCopy(const Vec& v, int n);
        ComplexMatrix ThetaCopy(const Vec& v, int n);

        void AddAccel(Particle& par,
                const ComplexMatrix& F) const override;
        void P2M(Octree* leaf) override;
        void M2M(Octree const* child, Octree* parent) override;
        ComplexMatrix M2X(Octree const* source, const Vec& s) override;
        ComplexMatrix L2X(Octree const* previous, const Vec& sp) override;
};

}
#endif

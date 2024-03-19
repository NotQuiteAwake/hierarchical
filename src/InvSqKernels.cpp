#include <cassert>
#include <cmath>
#include "boost/math/special_functions/spherical_harmonic.hpp"
#include "boost/math/special_functions/factorials.hpp"
#include "InvSqKernels.hpp"

namespace sim {

static constexpr double PI() { return std::atan(1) * 4; }

InvSqKernels::InvSqKernels(int p): Kernels(p) {};

double InvSqKernels::Prefactor(int n, int m) const {
    using boost::math::factorial;
    return sqrt(factorial<double>(n - m) * factorial<double>(n + m));
}

// boost_Y legendre includes sign term
// but still differs by a factor, and a sign...?
// further note the definitions of theta and phi.
InvSqKernels::cdouble InvSqKernels::Y(const Vec& v, int n, int m) const {
    using boost::math::spherical_harmonic;
    double theta = v.GetTheta();
    double phi = v.GetPhi();
    cdouble boost_Y = spherical_harmonic(n, m, theta, phi);
    return boost_Y * sqrt(4 * PI() / (2 * n + 1)) * pow(-1, m);
}

InvSqKernels::cdouble InvSqKernels::Gamma(const Vec& v, int n, int m) const {
    double r = v.GetNorm();
    return 1.0 / Prefactor(n, m) * pow(r, n) * Y(v, n, m);
}

InvSqKernels::cdouble InvSqKernels::Theta(const Vec& v, int n, int m) const {
    double r = v.GetNorm();
    assert(r);
    return Prefactor(n, m) * pow(r, -n - 1) * Y(v, n, m);

}

void InvSqKernels::AddAccel(Particle& par,
        const ComplexMatrix& psi) const {
    // note psi[1][0] should just be real.
    par.accel -= Vec({psi[1][1].real(), psi[1][1].imag(), psi[1][0].real()})
        * par.GetCharge() / par.GetMass();
}

void InvSqKernels::P2M(Octree* leaf) const {
    assert(leaf);
    assert(leaf->IsLeaf());

    for (int n = 0; n <= mP; n++) {
        for (int m = -n; m <= n; m++) {
            for (int soul : leaf->mSouls) {
                const Particle& par = leaf->GetParticle(soul);
                const double q = par.GetCharge();
                const Vec& pos = par.pos;
                const Vec& com = leaf->com;
                leaf->M[n][m] += q * Gamma(pos - com, n, m);
            }
        }
    }
    
}

void InvSqKernels::M2M(Octree const* child, Octree* parent) const {
    assert(parent && child);
    assert(parent == child->GetParent());

    const Vec& zp = parent->com;
    const Vec& z = child->com;
    for (int n = 0; n <= mP; n++) {
        for (int m = -n; m <= n; m++) {
            // I used to loop over all children here when M2M was meant to
            // propagate one child -> one parent only. Royal pain.
            for (int k = 0; k <= n; k++) {
                for (int l = -k; l <= k; l++) {
                    parent->M[n][m] += Gamma(z - zp, k, l) * 
                        child->M[n - k][m - l];
                }
            }
        }
    }
}

ComplexMatrix InvSqKernels::M2X(Octree const* source, const Vec& s) const {
    assert(source);
    ComplexMatrix F(mP, mP);

    const Vec& z = source->com;
    for (int n = 0; n <= mP; n++) {
        for (int m = -n; m <= n; m++) {
            for (int k = 0; k <= mP - n; k++) {
                for (int l = -k; l <= k; l++) {
                    F[n][m] += std::conj(source->M[k][l]) 
                        * Theta(s - z, n + k, m + l);
                }
            }
        }
    }

    return F;
}

ComplexMatrix InvSqKernels::L2X(Octree const* previous, const Vec& sp) const {
    assert(previous);
    ComplexMatrix F(mP, mP);

    const Vec& s = previous->com;
    for (int n = 0; n <= mP; n++) {
        for (int m = -n; m <= n; m++) {
            for (int k = 0; k < mP - n; k++) {
                for (int l = -k; l <= k; l++) {
                    F[n][m] += std::conj(Gamma(s - sp, k, l))
                        * previous->F[n + k][m + l];
                }
            }
        }
    }
    return F;
}

}

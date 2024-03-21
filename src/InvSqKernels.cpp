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
    return boost_Y * sqrt(4 * PI() / (2 * n + 1)) * (m % 2 ? -1.0 : 1.0);
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

// Spherical harmonics via recurrence reln, ref Dehnen 2014 A.4, A.5
ComplexMatrix InvSqKernels::Gamma(const Vec& v, int n) const {
    typedef InvSqKernels::cdouble cdouble;
    using namespace std::complex_literals;

    double r = v.GetNorm();
    // gamma could represent going from the only mass to com - r = 0. so don't
    // assert(r);

    ComplexMatrix mat(n, n);

    // initialise on a boundary "lower triangle".
    double x = v[0], y = v[1], z = v[2];
    mat[0][0] = cdouble(1);
    mat[1][0] = cdouble(z);
    mat[1][1] = (x + y * 1i) / 2.0;
    
    // apply recurrence relations
    // NOTE n is not loop variable! replace n with i!
    for (int i = 2; i <= n; i++) {
        for (int j = 0; j <= n; j++) {
            if (i == j) {
                mat[i][j] = (x + y * 1i) / (2.0 * i) * mat[i - 1][j - 1];
            } else {
                mat[i][j] = 1.0 / (i * i - j * j) *
                   ((2 * i - 1) * z * mat[i - 1][j] - r * r * mat[i - 2][j]);
            }
        }
    }

    // enforce symmetry to get other half (from observing calculated values)
    for (int i = 1; i <= n; i++) {
        for (int j = -n; j < 0; j++) {
            if (j % 2) { mat[i][j] = -std::conj(mat[i][-j]); } // odd j
            else { mat[i][j] = std::conj(mat[i][-j]); } // even j
        }
    }
    return mat;
}

ComplexMatrix InvSqKernels::Theta(const Vec& v, int n) const {
    typedef InvSqKernels::cdouble cdouble;
    using namespace std::complex_literals;
    double r = v.GetNorm();
    assert(r);

    ComplexMatrix mat(n, n);

    double x = v[0], y = v[1], z = v[2];

    mat[0][0] = cdouble(1.0 / r);
    mat[1][0] = cdouble(z / pow(r, 3));
    mat[1][1] = (2.0 * 1 - 1) * (x + 1i * y) / (r * r) * mat[0][0];
    for (int i = 2; i <= n; i++) {
        for (int j = 0; j <= i; j++) {
            if (i == j) {
                mat[i][j] = (2.0 * i - 1) * (x + 1i * y) / (r * r)
                    * mat[i - 1][j - 1]; }
            else {
                mat[i][j] = ((2 * i - 1) * z * mat[i - 1][j]
                              - 1.0 * ((i - 1) * (i - 1) - j * j) * mat[i - 2][j])
                    * 1.0 / (r * r);
            } 
        } 
    }

    for (int i = 1; i <= n; i++) {
        for (int j = -i; j < 0; j++) {
            if (j % 2) { mat[i][j] = -std::conj(mat[i][-j]); }
            else { mat[i][j] = std::conj(mat[i][-j]); }
        }
    }

    return mat;
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

    // TODO: test effect of using new Gamma here
    // This doesn't appear to be the bottleneck for now though...

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
    const ComplexMatrix gamma = Gamma(z - zp, mP);

    for (int n = 0; n <= mP; n++) {
        for (int m = -n; m <= n; m++) {
            // I used to loop over all children here when M2M was meant to
            // propagate one child -> one parent only. Royal pain.
            for (int k = 0; k <= n; k++) {
                for (int l = -k; l <= k; l++) {
                    // originally Gamma() in 4th order loop.
                    // Could have been very costly
                    parent->M[n][m] += gamma[k][l] * 
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
    const ComplexMatrix theta = Theta(s - z, mP);
    for (int n = 0; n <= mP; n++) {
        for (int m = -n; m <= n; m++) {
            for (int k = 0; k <= mP - n; k++) {
                for (int l = -k; l <= k; l++) {
                    F[n][m] += std::conj(source->M[k][l]) 
                        * theta[n + k][m + l];
                    // highest 1st index is still just n
                }
            }
        }
    }

    return F;
}

ComplexMatrix InvSqKernels::L2X(Octree const* previous, const Vec& sp) const {
    assert(previous);
    const Vec& s = previous->com;

    ComplexMatrix F(mP, mP);
    ComplexMatrix gamma = Gamma(s - sp, mP);

    for (int n = 0; n <= mP; n++) {
        for (int m = -n; m <= n; m++) {
            for (int k = 0; k < mP - n; k++) {
                for (int l = -k; l <= k; l++) {
                    F[n][m] += std::conj(gamma[k][l])
                        * previous->F[n + k][m + l];
                }
            }
        }
    }
    return F;
}

}

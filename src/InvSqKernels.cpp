/**
 * @file
 * @brief Implementation of the inverse square multipole expansion kernels.
 */

#include <cassert>
#include <cmath>
#include "boost/math/special_functions/spherical_harmonic.hpp"
#include "boost/math/special_functions/factorials.hpp"
#include "InvSqKernels.hpp"

namespace sim {

/**
 * @class
 * @brief Implementation of the inverse square multipole expansion kernels.
 *
 * As in the Kernels interface, non-const reference or pointer parameters in
 * method signature are usually subject to direct modification, as they are the
 * the "target" of the operation.
 */

static constexpr double PI() { return std::atan(1) * 4; }

/**
 * @brief Initialise InvSqKernels.
 *
 * @param[in] p Order of multipole expansion.
 * @param[in] G Coupling constant in inverse-square law force.
 */
InvSqKernels::InvSqKernels(int p, double G):
    Kernels(p),
    mG(G),
    mTempMatrix(ComplexMatrix(p, p)) {};

double InvSqKernels::Prefactor(int n, int m) const {
    using boost::math::factorial;
    return sqrt(factorial<double>(n - m) * factorial<double>(n + m));
}

/**
 * @brief Calculate surface spherical harmonics via boost
 *
 * My definition of the harmonics here agrees with Dehnen 2014.
 *
 * Interestingly boost documentation says their Legendre has absorbed a sign
 * term, but turns out I still need an extra sign term to get it right. I
 * further need to fix its prefactor to agree with Dehnen 2014.
 *
 * The angle definitions follow those in boost.
 */
InvSqKernels::cdouble InvSqKernels::Y(const Vec& v, int n, int m) const {
    using boost::math::spherical_harmonic;
    double theta = v.GetTheta();
    double phi = v.GetPhi();
    cdouble boost_Y = spherical_harmonic(n, m, theta, phi);
    return boost_Y * sqrt(4 * PI() / (2 * n + 1)) * (m % 2 ? -1.0 : 1.0);
}

/**
 * @brief Calculate solid harmonics gamma from boost functions
 */
InvSqKernels::cdouble InvSqKernels::GammaBoost(const Vec& v, int n, int m) const {
    double r = v.GetNorm();
    return 1.0 / Prefactor(n, m) * pow(r, n) * Y(v, n, m);
}

/**
 * @brief Calculate solid harmonics theta from boost functions
 */
InvSqKernels::cdouble InvSqKernels::ThetaBoost(const Vec& v, int n, int m) const {
    double r = v.GetNorm();
    assert(r);
    return Prefactor(n, m) * pow(r, -n - 1) * Y(v, n, m);
}

/**
 * @brief Calculate a specific solid harmonics gamma via recurrence.
 *
 * Kept for completeness, but there is no need for this - in all applications we
 * just use the preprocessed full matrix which is quicker.
 */
InvSqKernels::cdouble InvSqKernels::Gamma(const Vec& v, int n, int m) const {
    typedef InvSqKernels::cdouble cdouble;
    using namespace std::complex_literals;

    int abs_m = std::abs(m);
    if (abs_m > n) return 0;
    if (useBoostLimit >= 0 && (n + abs_m) >= useBoostLimit) {
        return GammaBoost(v, n, m);
    }

    // O(n + m) evaluation via recurrence
    cdouble cur = 1.0;
    double x = v[0], y = v[1], z = v[2];
    double r = v.GetNorm();

    // get to theta(abs_m, abs_m) first
    for (int i = 1; i <= abs_m; i++) {
        cur = (x + y * 1i) / (2.0 * i) * cur;
    }


    cdouble pre1 = 0; // mat[i - 1][j]
    cdouble pre2 = 0; // mat[i - 2][j]

    if (n >= 1) {
        int start = abs_m + 1; // we are at (abs_m, abs_m) right now.
        if (m == 0) {
            // the recursion formula does not work for n = 1, m = 0.
            // So we do it first and manually displace start
            // this is fine, since we are inside {if (n >= 1)}
            pre1 = cur;
            cur = cdouble(z);
            start = 2;
        }
        for (int i = start; i <= n; i++) {
            pre2 = pre1;
            pre1 = cur;
            cur = 1.0 / (i * i - abs_m * abs_m) *
                ((2 * i - 1) * z * pre1 - r * r * pre2);
        }
    }

    if (m < 0) {
        if (m % 2) { cur = -std::conj(cur); } // odd j
        else { cur = std::conj(cur); } // even j
    }

    return cur;
}

/**
 * @brief Calulate a single solid harmonics theta via recurrence.
 */
InvSqKernels::cdouble InvSqKernels::Theta(const Vec& v, int n, int m) const {
    typedef InvSqKernels::cdouble cdouble;
    using namespace std::complex_literals;

    int abs_m = std::abs(m);
    if (abs_m > n) return 0;
    if (useBoostLimit >= 0 && (n + abs_m) >= useBoostLimit) {
        return ThetaBoost(v, n, m);
    }

    // O(n + m) evaluation via recurrence
    const double x = v[0], y = v[1], z = v[2];
    const double r = v.GetNorm();
    const double r2 = r * r;

    cdouble cur = 1.0 / r;

    // go along diagonal
    for (int i = 1; i <= abs_m; i++) {
        cur = (2.0 * i - 1) * (x + 1i * y) / r2 * cur;
    }
    
    cdouble pre1 = 0; // mat[i - 1][j]
    cdouble pre2 = 0; // mat[i - 2][j]

    if (n >= 1) {
        int start = abs_m + 1; // we are at (abs_m, abs_m) right now.
        if (m == 0) {
            // the recursion formula does not work for n = 1, m = 0.
            // So we do it first and manually displace start
            // this is fine, since we are inside {if (n >= 1)}
            pre1 = cur;
            cur = cdouble(z / (r2 * r));
            start = 2;
        }
        for (int i = start; i <= n; i++) {
            pre2 = pre1;
            pre1 = cur;
            cur = ((2 * i - 1) * z * pre1
                    - 1.0 * ((i - 1) * (i - 1) - abs_m * abs_m) * pre2) * 1.0 / r2;
        }
    }


    if (m < 0) {
        if (m % 2) { cur = -std::conj(cur); }
        else { cur = std::conj(cur); }
    }

    return cur;
}

/**
 * @brief pre-calculate solid harmonics gamma of order n.
 *
 * The recurrence relations follow from Dehnen 2014 A.4 A.5.
 *
 * The values are stored in an internal temporary matrix to reduce allocations,
 * and is made available to all internal functions.
 */
void InvSqKernels::Gamma(const Vec& v, int n) {
    typedef InvSqKernels::cdouble cdouble;
    using namespace std::complex_literals;

    double r = v.GetNorm();
    // gamma could represent going from the only mass so coc - r = 0. so don't
    // assert(r);
    assert(n <= mP);

    ComplexMatrix& mat = mTempMatrix;

    // initialise on a boundary "lower triangle".
    double x = v[0], y = v[1], z = v[2];
    double r2 = r * r;
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
                   ((2 * i - 1) * z * mat[i - 1][j] - r2 * mat[i - 2][j]);
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
}

/**
 * @brief pre-calculate solid harmonics theta of order n.
 *
 * The recurrence relations follow from Dehnen 2014 A.4 A.5.
 *
 * The values are stored in an internal temporary matrix to reduce allocations,
 * and is made available to all internal functions.
 */
void InvSqKernels::Theta(const Vec& v, int n) {
    typedef InvSqKernels::cdouble cdouble;
    using namespace std::complex_literals;
    double r = v.GetNorm();
    assert(r);
    assert(n <= mP);

    ComplexMatrix& mat = mTempMatrix;

    double x = v[0], y = v[1], z = v[2];
    double r2 = r * r;

    mat[0][0] = cdouble(1.0 / r);
    mat[1][0] = cdouble(z / (r2 * r));
    mat[1][1] = (2.0 * 1 - 1) * (x + 1i * y) / r2 * mat[0][0];
    for (int i = 2; i <= n; i++) {
        for (int j = 0; j <= i; j++) {
            if (i == j) {
                mat[i][j] = (2.0 * i - 1) * (x + 1i * y) / r2 
                    * mat[i - 1][j - 1]; }
            else {
                mat[i][j] = ((2 * i - 1) * z * mat[i - 1][j]
                              - 1.0 * ((i - 1) * (i - 1) - j * j) * mat[i - 2][j])
                    * 1.0 / r2;
            } 
        } 
    }

    for (int i = 1; i <= n; i++) {
        for (int j = -i; j < 0; j++) {
            if (j % 2) { mat[i][j] = -std::conj(mat[i][-j]); }
            else { mat[i][j] = std::conj(mat[i][-j]); }
        }
    }
}

/**
 * @brief Calculate and return a copy of gamma coefficients of order n
 *
 * @return solid harmonics gamma of order n.
 */
ComplexMatrix InvSqKernels::GammaCopy(const Vec& v, int n) {
    Gamma(v, n);
    return mTempMatrix;
}

/**
 * @brief Calculate and return a copy of theta coefficients of order n
 *
 * @return solid harmonics theta of order n.
 */
ComplexMatrix InvSqKernels::ThetaCopy(const Vec& v, int n) {
    Theta(v, n);
    return mTempMatrix;
}

/**
 * @brief Calculate and add accel and pot due to psi to particle.
 * 
 * This follows from expression at top of page 3, Dehnen 2014. 
 *
 * @param[in, out] par Target particle.
 * @param[in] psi \f$F_n^m\f$ at particle position.
 */
void InvSqKernels::AddAccel(Particle& par,
        const ComplexMatrix& psi) const {
    // note psi[1][0] should just be real.
    par.accel += Vec({psi[1][1].real(), psi[1][1].imag(), psi[1][0].real()})
        * mG * par.GetCharge() / par.GetMass();
    par.pot += mG * par.GetCharge() * psi[0][0].real();
}

/**
 * @brief P2M kernel. (Eq 3c, Dehnen 2014)
 *
 * @param[in, out] leaf Leaf node to apply P2M kernel to.
 */
void InvSqKernels::P2M(Octree* leaf) {
    assert(leaf);
    assert(leaf->IsLeaf());

    const Vec& coc = leaf->coc;

    for (int soul : leaf->mSouls) {
        const Particle& par = leaf->GetParticle(soul);
        const double q = par.GetCharge();
        const Vec& pos = par.pos;
        
        Gamma(pos - coc, mP);
        const ComplexMatrix& gamma = mTempMatrix;
        
        for (int n = 0; n <= mP; n++) {
            for (int m = -n; m <= n; m++) {
                leaf->M[n][m] += q * gamma[n][m];
            }
        }
    }
    
}

/**
 * @brief M2M kernel (Eq 3d, Dehnen 2014)
 *
 * @param[in] child Child node with known M coefficients
 * @param[in, out] parent parent node to add shifted M coefficients to
 */
void InvSqKernels::M2M(Octree const* child, Octree* parent) {
    assert(parent && child);
    assert(parent == child->GetParent());

    const Vec& zp = parent->coc;
    const Vec& z = child->coc;
    Gamma(z - zp, mP);
    const ComplexMatrix& gamma = mTempMatrix;

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

/**
 * @brief Implements M2L, M2P kernels (Eq 3b, Dehnen 2014)
 *
 * M2L and M2P both finds F coefficients at a target location, so are the same
 * thing.
 *
 * @param[in] source Node with known M coefficients
 * @param[in] Target position
 * @return F coefficients at s position
 */
ComplexMatrix InvSqKernels::M2X(Octree const* source, const Vec& s) {
    assert(source);
    ComplexMatrix F(mP, mP);

    const Vec& z = source->coc;
    Theta(s - z, mP);
    const ComplexMatrix& theta = mTempMatrix;
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

/**
 * @brief Implements L2L, L2P kernels (Eq 3f, Dehnen 2014)
 *
 * @param[in] previous Node with known F coefficients
 * @param[in] sp Target location
 * @return F coefficients at sp
 */
ComplexMatrix InvSqKernels::L2X(Octree const* previous, const Vec& sp) {
    assert(previous);
    const Vec& s = previous->coc;

    ComplexMatrix F(mP, mP);
    Gamma(s - sp, mP);
    const ComplexMatrix& gamma = mTempMatrix;

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

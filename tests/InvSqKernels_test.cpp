#include "doctest/doctest.h"
#include "InvSqKernels.hpp"

namespace sim {

double Factorial(int n) {
    assert(n >= 0);
    if (n == 0) { return 1; }
    return Factorial(n - 1) * n; 
}

double Prefactor(int n, int m) {
    return Factorial(n - m) * Factorial(n + m);
}

double U33(const Vec& vec) {
    const double x = vec[0];
    const double y = vec[1];
    return 15 * x * (pow(x, 2) - 3 * pow(y, 2)) / Prefactor(3, 3);
}

double U3_2(const Vec& vec) {
    return 30 * vec[0] * vec[1] * vec[2] / Prefactor(3, -2);
}

double U2_1(const Vec& vec) {
    return (3 * vec[1] * vec[2]) / Prefactor(2, -1);
}

double T_factor(const Vec& vec, int n, int m) {
    return Prefactor(n, m) / std::pow(vec.GetNorm(), 2 * n + 1);
}

double T33(const Vec& vec) {
    return U33(vec) * T_factor(vec, 3, 3);
}

double T3_2(const Vec& vec) {
    return U3_2(vec) * T_factor(vec, 3, -2);
}

double T2_1(const Vec& vec) {
    return U2_1(vec) * T_factor(vec, 2, -1);
}

double GammaU(const Vec& vec, int n, int m, InvSqKernels& invsq) {
    using namespace std::complex_literals;
    typedef std::complex<double> cdouble;

    cdouble gamma = invsq.Gamma(vec, n, abs(m));

    if (m >= 0) {
        return gamma.real();
    }

    return gamma.imag();
}

double ThetaT(const Vec& vec, int n, int m, InvSqKernels& invsq) {
    using namespace std::complex_literals;
    typedef std::complex<double> cdouble;

    cdouble theta = invsq.Theta(vec, n, abs(m));
    
    if (m >= 0) {
        return theta.real();
    }
    
    return theta.imag();
}

void CHECK_CDOUBLE(std::complex<double> z1, std::complex<double> z2) {
    CHECK(z1.real() == doctest::Approx(z2.real()));
    CHECK(z1.imag() == doctest::Approx(z2.imag()));
}

TEST_CASE("test InvSqKernels") {
    typedef std::complex<double> cdouble;

    int n_max = 10;
    InvSqKernels invsq(n_max); // number doesn't matter in tests below.

    const Vec vec({1, 2, 3});
    
    SUBCASE("Solid harmonics gamma") {
        // U(3, 3) U as per Dehnen 2014 P19
        double Gamma_U33 = GammaU(vec, 3, 3, invsq);
        CHECK(U33(vec) == doctest::Approx(Gamma_U33));
        // U(3, -2)
        double Gamma_U3_2 = GammaU(vec, 3, -2, invsq);
        CHECK(U3_2(vec) == doctest::Approx(Gamma_U3_2));
        // U(2, -1)
        double Gamma_U2_1 = GammaU(vec, 2, -1, invsq);
        CHECK(U2_1(vec) == doctest::Approx(Gamma_U2_1));
    }

    SUBCASE("Solid harmonics theta") {
        // T(3, 3)
        double Theta_T33 = ThetaT(vec, 3, 3, invsq);
        CHECK(T33(vec) == doctest::Approx(Theta_T33));
        // T(3, -2)
        double Theta_T3_2 = ThetaT(vec, 3, -2, invsq);
        CHECK(T3_2(vec) == doctest::Approx(Theta_T3_2));
        // T(2, -1)
        double Theta_T2_1 = ThetaT(vec, 2, -1, invsq);
        CHECK(T2_1(vec) == doctest::Approx(Theta_T2_1));
    }

    SUBCASE("preprocessing & recurrence methods") {
        ComplexMatrix mtheta = invsq.ThetaCopy(vec, n_max); // m for matrix
        ComplexMatrix mgamma = invsq.GammaCopy(vec, n_max);
        for (int n = 0; n <= n_max; n++) {
            for (int m = -n; m <= n; m++) {
                cdouble rtheta = invsq.Theta(vec, n, m); // r for recurrence
                cdouble rgamma = invsq.Gamma(vec, n, m);
                CHECK_CDOUBLE(mtheta[n][m], rtheta);
                CHECK_CDOUBLE(mgamma[n][m], rgamma);
            }
        }
    }

}

}

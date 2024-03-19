#include "doctest/doctest.h"
#include "InvSqKernels.hpp"

namespace sim {

double U33(const Vec& vec) {
    const double& x = vec[0];
    const double& y = vec[1];
    return 15 * x * (pow(x, 2) - 3 * pow(y, 2));
}

double U2_1(const Vec& vec) {
    return 3 * vec[1] * vec[2];
}

double GammaU(const Vec& vec, int n, int m, const InvSqKernels& invsq) {
    using namespace std::complex_literals;
    typedef std::complex<double> cdouble;

    cdouble gamma = invsq.Gamma(vec, n, abs(m));
    cdouble gamma_conj = std::conj(gamma);
    double sign = (m >= 0) ? 1.0 : -1.0;
    cdouble factor = (m >= 0) ? (1.0 / 2) : (1.0 / (2i));

    return (factor * (gamma + sign * gamma_conj)).real()
        * pow(invsq.Prefactor(n, m), 2);
}

TEST_CASE("test InvSqKernels") {
    const InvSqKernels invsq(3); // number doesn't matter in tests below.
    const Vec vec({1, 2, 3});
    
    SUBCASE("test Prefactor") {
        CHECK(invsq.Prefactor(3, 2) ==
                doctest::Approx(sqrt(1 * 5 * 4 * 3 * 2 * 1)));
    }

    SUBCASE("test spherical harmonics Y and Gamma") {
        // U(3, 3) U as per Dehnen 2014 P19
        double Gamma_U33 = GammaU(vec, 3, 3, invsq);
        CHECK(U33(vec) == doctest::Approx(Gamma_U33));
        // U(2, -1)
        double Gamma_U2_1 = GammaU(vec, 2, -1, invsq);
        CHECK(U2_1(vec) == doctest::Approx(Gamma_U2_1));
    }

}

}

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

TEST_CASE("test InvSqKernels") {
    using namespace std::complex_literals;

    const Vec vec({1, 2, 3});

    InvSqKernels invsq(3);
    
    SUBCASE("test Prefactor") {
        CHECK(invsq.Prefactor(3, 2) ==
                doctest::Approx(sqrt(1 * 5 * 4 * 3 * 2 * 1)));
    }

    SUBCASE("test spherical harmonics Y and Gamma") {
        // U(3, 3) U as per Dehnen 2014 P19
        double Gamma_U33 = (1.0 / (2) * (invsq.Gamma(vec, 3, 3)
                + std::conj(invsq.Gamma(vec, 3, 3)))).real()
            * pow(invsq.Prefactor(3, 3), 2);
        CHECK(U33(vec) == doctest::Approx(Gamma_U33));
        // U(2, -1)
        double Gamma_U2_1 = (1.0 / (2i) * (invsq.Gamma(vec, 2, 1)
                    - std::conj(invsq.Gamma(vec, 2, 1)))).real()
                * pow(invsq.Prefactor(2, -1), 2);
        CHECK(U2_1(vec) == doctest::Approx(Gamma_U2_1));
    }

}

}

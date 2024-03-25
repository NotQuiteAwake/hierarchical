#include <cmath>
#include "doctest/doctest.h"
#include "Vec.hpp"

namespace sim {

static constexpr double PI() { return std::atan(1)*4; }

TEST_CASE("testing vector class") {
    Vec v1 = Vec({1, 2, 3});
    Vec v2 = Vec({1, 1, 1});

    // awful things could happen if it was not set by default (assumed
    // behaviour) as we are dynamically allocating the memory
    SUBCASE("test ctor initial value") {
        Vec v = Vec();
        for (int i = 0; i < Vec::mDim; i++) {
            CHECK(v[i] == 0);
        }
    }

    SUBCASE("vector inner product") {
        CHECK(v1.GetNorm() == doctest::Approx(sqrt(14)));
    }

    SUBCASE("vector cross product") {
        Vec res = CrossProduct(v1, v2);
        doctest::Approx exp[3] = {-1, 2, -1};
        for (int i = 0; i < Vec::mDim; i++) {
            CHECK(res[i] == exp[i]);
        }
    }

    SUBCASE("vector assignment") {
        v1 = v2;
        CHECK(v1[2] == v2[2]);
    }

    SUBCASE("vector addition") {
        Vec res = v1 + v2;
        doctest::Approx exp[3] = {2, 3, 4};
        for (int i = 0; i < Vec::mDim; i++) {
            CHECK(res[i] == exp[i]);
        }
    }

    SUBCASE("spherical polar theta") {
        Vec vec({1, 0, 1});
        CHECK(vec.GetTheta() == doctest::Approx(PI() / 4));
        vec = Vec({1, 0, -1});
        CHECK(vec.GetTheta() == doctest::Approx(PI() * 3 / 4));
    }

    SUBCASE("spherical polar phi") {
        Vec vec({1, 1, 0});
        CHECK(vec.GetPhi() == doctest::Approx(PI() / 4));
        vec = Vec({-1, std::sqrt(1.0 / 3), 1});
        CHECK(vec.GetPhi() == doctest::Approx(5 * PI() / 6));
        vec = Vec({-1, -std::sqrt(1.0 / 3), 1});
        CHECK(vec.GetPhi() == doctest::Approx(7 * PI() / 6));
        vec = Vec({std::sqrt(3), -1, 3});
        CHECK(vec.GetPhi() == doctest::Approx(11 * PI() / 6));
    }

}

}

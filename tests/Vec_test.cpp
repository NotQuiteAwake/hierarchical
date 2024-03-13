#include <cmath>
#include "doctest/doctest.h"
#include "Vec.hpp"

namespace sim {

TEST_CASE("testing vector class") {
    Vec v1 = Vec({1, 2, 3});
    Vec v2 = Vec({1, 1, 1});

    SUBCASE("vector inner product") {
        CHECK(v1.CalculateNorm() == doctest::Approx(sqrt(14)));
    }

    SUBCASE("vector cross product") {
        CHECK(CrossProduct(v1, v2)[0] == doctest::Approx(-1));
    }

    SUBCASE("vector assignment") {
        v1 = v2;
        CHECK(v1[2] == v2[2]);
    }

    SUBCASE("vector addition") {
        Vec v3 = v1 + v2;
        CHECK(v3[2] == 4);
    }

}

}

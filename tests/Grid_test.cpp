#include "doctest/doctest.h"
#include "Octant.hpp"
#include "Vec.hpp"

namespace sim {

TEST_CASE("test Octant") {
    Octant octant({{-1, 1}, {-1, 1}, {-1, 1}});

    SUBCASE("test GetOctantNumber") {
        Vec v({0.5, 0.5, 0.5});
        CHECK(octant.GetOctantNumber(v) == 0b111);
        v = Vec({-0.5, 0.5, -0.5});
        CHECK(octant.GetOctantNumber(v) == 0b010);
    } 

    SUBCASE("test GetOctant") {
        Octant new_oct = octant.GetOctant(0b010);
        CHECK(new_oct[1][0] == doctest::Approx(0));
        CHECK(new_oct[1][1] == doctest::Approx(1));
    }


}

}

#include "doctest/doctest.h"
#include "Octant.hpp"
#include "Vec.hpp"

namespace sim {

TEST_CASE("test Octant") {
    Octant octant({{-1, 1}, {-1, 1}, {-1, 1}});

    SUBCASE("test operator==") {
        Octant other_octant({{-1, 1}, {-1, 1}, {-1, 1}});
        CHECK(octant == other_octant);
    }

    SUBCASE("test GetOctantNumber") {
        Vec v({0.5, 0.5, 0.5});
        CHECK(octant.GetOctantNumber(v) == 0b111);
        v = Vec({-0.5, 0.5, -0.5});
        CHECK(octant.GetOctantNumber(v) == 0b010);
    } 

    SUBCASE("test GetOctant") {
        Octant new_oct = octant.GetOctant(0b010);
        Octant expected_oct = Octant({{-1, 0}, {0, 1}, {-1, 0}});
        CHECK(new_oct == expected_oct);
    }


}

}

#include "doctest/doctest.h"
#include "Grid.hpp"
#include "Octant.hpp"
#include "Particle.hpp"

namespace sim {

void AddParticles(int size, Grid& grid) {
    for (int i = 0; i < size; i++) {
        Particle par(1, 1);
        double x = double(i);
        par.pos = Vec({x, x, x});
        grid.AddParticle(par);
    }
}

TEST_CASE("test Grid") {
    SUBCASE("test Grid with maximum limits") {
        Octant max_lim({{1, 4}, {1, 4}, {1, 4}});
        Grid grid(max_lim);

        AddParticles(5, grid);

        CHECK(grid.GetSize() == 3); // first and final particle not in
        CHECK(grid[0].pos[2] == 1);
    }

    SUBCASE("test Grid with no max limits") {
        Grid grid;

        AddParticles(5, grid);

        CHECK(grid.GetSize() == 5);
        CHECK(grid[4].pos[2] == 4);
    }

}

}

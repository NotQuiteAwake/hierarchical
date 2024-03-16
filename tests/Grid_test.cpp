#include "doctest/doctest.h"
#include "Grid.hpp"
#include "Octant.hpp"
#include "Particle.hpp"

namespace sim {

TEST_CASE("test Grid") {
    Octant octant = Octant();
    for (int i = 0; i < Octant::mDim; i++) {
        octant[i][0] = -1;
        octant[i][1] = 1;
    }
    Grid grid = Grid(octant);

    int size = 5;
    for (int i = 0; i < size; i++) {
        Particle par(1, 1);
        double x = double(i);
        par.pos = Vec({x, x, x});
        grid.AddParticle(par);
    }

    SUBCASE("test operator[]") {
        CHECK(grid[3].pos[2] == 3);
    }

    SUBCASE("test GetSize()") {
        CHECK(grid.GetSize() == size);
    }

}

}

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

void CheckVec(Vec vec, const double (&exp)[Vec::mDim]) {
    for (int i = 0; i < Vec::mDim; i++) {
        CHECK(vec[i] == doctest::Approx(exp[i]));
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
    
    SUBCASE("test Grid getters") {
        Grid grid;

        Particle p1(1, 1);
        p1.pos = Vec({1, 0, 0});
        p1.vel = Vec({0, -1, 0});
        p1.pot = -1;
        
        Particle p2(2, 1);
        p2.pos = Vec({0, -1, 0});
        p2.vel = Vec({-0.5, 0, 0});
        p2.pot = -2;

        grid.AddParticle(p1);
        grid.AddParticle(p2);

        SUBCASE("test Grid GetL") {
            double exp[3] = {0, 0, -2};
            CheckVec(grid.GetL(), exp);
        }

        SUBCASE("test Grid GetP") {
            double exp[3] = {-1, -1, 0};
            CheckVec(grid.GetP(), exp);
        }

        SUBCASE("test Grid GetE") {
            CHECK(grid.GetPE() == doctest::Approx(-1.5));
            CHECK(grid.GetKE() == doctest::Approx(1.0 / 4 + 1.0 / 2));
            CHECK(grid.GetE() == doctest::Approx(-1.5 + 1.0 / 4 + 1.0 / 2));
        }

        SUBCASE("test Grid GetCOM") {
            double exp[3] = {1.0 / 3, -2.0 / 3, 0};
            CheckVec(grid.GetCOM(), exp);
        }

    }
        

}

}

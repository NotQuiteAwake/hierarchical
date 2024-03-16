#include <memory>
#include <random>
#include "doctest/doctest.h"
#include "Grid.hpp"
#include "Brute.hpp"
#include "Force.hpp"
#include "Barnes-Hut.hpp"

namespace sim {

Grid MakeTriangleGrid(const Grid& g1) {
    // in practice we want to set charge as gravitational charge.
    // here set charge to zero to make sure we are just calculating gravity from
    // masses alone.
    Particle p[3] = {Particle({2, 0}), Particle({2, 0}), Particle({2, 0})};
    p[0].pos = Vec({1, sqrt(3), 0});
    p[1].pos = Vec({0, 0, 0});
    p[2].pos = Vec({2, 0, 0});

    Grid grid = g1;
    for (int i = 0; i < 3; i++) {
        grid.AddParticle(p[i]);
    }
    return grid;
}

void TestTriangle(std::shared_ptr<Interaction> interaction) {
    const doctest::Approx TriangleExp[3][3] =
    {   {0, -sqrt(3) / 2, 0},
        {3.0 / 4, sqrt(3) / 4, 0},
        {-3.0 / 4, sqrt(3) / 4, 0}
    };
    const Octant octant({{-1, 5}, {-1, 5}, {-1, 5}});

    Grid grid = Grid(octant);
    grid = MakeTriangleGrid(grid);
    MakeTriangleGrid(grid);
    grid = interaction->Calculate(grid);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            CHECK(grid[i].accel[j] == TriangleExp[i][j]);
        }
    }
}

Grid MakeCluster(const Grid& g1, const Vec& centre,
        const double spread, const int num) {
    std::mt19937 mt(0);
    std::uniform_real_distribution<double> dist(-1, 1);

    Grid grid = g1;
    for (int i = 0; i < num; i++) {
        Particle par(dist(mt) * 2, dist(mt) * 2);
        par.pos = centre + 
            Vec({dist(mt) * spread, dist(mt) * spread, dist(mt) * spread});
        grid.AddParticle(par);
    }

    return grid;
}


// this ASSUMES Brute is fully functional and uses it as baseline.
void TestWellSeparated(std::shared_ptr<Interaction> interaction) {
    std::shared_ptr<Force> grav(new Gravity());
    std::shared_ptr<Interaction> brute(new Brute(grav));

    srand(0);
    
    const Octant octant({{-2, 10}, {-2, 10}, {-2, 10}});
    Grid grid = Grid(octant);

    grid = MakeCluster(grid, Vec({2, 2, 2}), 2, 10);
    grid = MakeCluster(grid, Vec({8, 8, 8}), 1, 5);
    grid = MakeCluster(grid, Vec({7, 2, 4}), 1, 5);
    grid = MakeCluster(grid, Vec({4, 6, 3}), 1, 6);

    Grid grid_brute = brute->Calculate(grid);
    Grid grid_alt = interaction->Calculate(grid);

    CHECK(grid_alt.GetSize() == grid_brute.GetSize());
    for (int i = 0; i < grid_alt.GetSize(); i++) {
        for (int j = 0; j < Vec::mDim; j++) {
            CHECK(grid_alt[i].pos[j] ==
                    doctest::Approx(grid_brute[i].pos[j]).epsilon(0.05));
        }
    }
}

TEST_CASE("test Brute with Gravity") {
    std::shared_ptr<Force> grav(new Gravity());
    std::shared_ptr<Interaction> brute(new Brute(grav));

    SUBCASE("test Triangle") {
        TestTriangle(brute);
    }
}

TEST_CASE("test Barnes-Hut with Gravity") {
    std::shared_ptr<Force> grav(new Gravity());
    std::shared_ptr<Interaction> brute(new Brute(grav));
    std::shared_ptr<Interaction> bh(new BarnesHut(1, 0.5, grav));

    SUBCASE("test trianglar masses") {
        TestTriangle(bh);
    }

    SUBCASE("test well-separated batches") {
        TestWellSeparated(bh);
    }
}

}

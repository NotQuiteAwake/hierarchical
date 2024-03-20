#include <doctest/doctest.h>
#include <memory>
#include <random>
#include "Grid.hpp"
#include "Brute.hpp"
#include "Force.hpp"
#include "InvSqKernels.hpp"
#include "Barnes-Hut.hpp"
#include "FMM.hpp"

namespace sim {

static double eps = 0.01;

Grid MakeTriangleGrid(const Grid& g1) {
    // The InvSqKernels is agnostic to the type of inverse square force we have
    // and instead uses charge. For this reason we must set the "gravitational
    // charge", which of course is just mass.
    Particle p[3] = {Particle({2, 2}), Particle({2, 2}), Particle({2, 2})};
    p[0].pos = Vec({1, sqrt(3), 0});
    p[1].pos = Vec({0, 0, 0});
    p[2].pos = Vec({2, 0, 0});

    Grid grid = g1;
    for (int i = 0; i < 3; i++) {
        grid.AddParticle(p[i]);
    }
    return grid;
}

void TestTriangle(Interaction const* interaction) {
    const double TriangleExp[3][3] =
    {   {0, -sqrt(3) / 2, 0},
        {3.0 / 4, sqrt(3) / 4, 0},
        {-3.0 / 4, sqrt(3) / 4, 0}
    };
    Grid grid;
    grid = MakeTriangleGrid(grid);
    grid = interaction->Calculate(grid);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            CHECK(grid[i].accel[j] ==
                    doctest::Approx(TriangleExp[i][j]).epsilon(eps));
        }
    }
}

Grid MakeCluster(const Grid& g1, const Vec& centre,
        const double spread, const int num) {
    std::mt19937 mt(0);
    std::uniform_real_distribution<double> mdist(0, 1); // mass dist
    std::uniform_real_distribution<double> pdist(-0.5, 0.5); // pos dist

    Grid grid = g1;
    for (int i = 0; i < num; i++) {
        double mass = mdist(mt);
        Particle par(mass, mass);
        par.pos = centre + 
            Vec({pdist(mt) * spread, pdist(mt) * spread, pdist(mt) * spread});
        grid.AddParticle(par);
    }

    return grid;
}

// this ASSUMES Brute is fully functional and uses it as baseline.
void TestWellSeparated(
        Interaction const* interaction,
        int numMass1,
        int numMass2
        ) {
    auto grav = std::make_unique<const Gravity>();
    auto brute = std::make_unique<const Brute>(grav.get());

    srand(0);
    
    Grid grid;

    grid = MakeCluster(grid, Vec({2, 2, 2}), 1, numMass1);
    grid = MakeCluster(grid, Vec({10, 10, 10}), 1, numMass2);
    // grid = MakeCluster(grid, Vec({7, 2, 4}), 1, 5);
    // grid = MakeCluster(grid, Vec({4, 6, 3}), 1, 6);

    Grid grid_brute = brute->Calculate(grid);
    Grid grid_alt = interaction->Calculate(grid);

    for (int i = 0; i < grid_alt.GetSize(); i++) {
        // so diagnostic message is more illustrating...
        for (int j = 0; j < 3; j++) {
            CHECK(grid_alt[i].accel[j] ==
                    doctest::Approx(grid_brute[i].accel[j]).epsilon(eps));
        }
    }
}

TEST_CASE("test inverse square Interactions") {
    int p = 6;
    double theta = 0.5;
    auto grav = std::make_unique<const Gravity>();
    auto invsq = std::make_unique<const InvSqKernels>(p);

    SUBCASE("test Brute with triangular masses") {
        auto brute = std::make_unique<Brute>(grav.get());
        TestTriangle(brute.get());
    } 

    SUBCASE("test Barnes-Hut") {
        auto bh = std::make_unique<BarnesHut>(
                p,
                theta,
                invsq.get(),
                grav.get());

        SUBCASE("test BH with triangular masses") {
            TestTriangle(bh.get());
        }

        SUBCASE("test BH with well-separated masses") {
            TestWellSeparated(bh.get(), 1, 1);
            TestWellSeparated(bh.get(), 1, 2);
            TestWellSeparated(bh.get(), 1, 10);
            TestWellSeparated(bh.get(), 10, 10);
        }
    } 

    SUBCASE("test FMM") {  
        const int fmm_params = 4;
        
        std::unique_ptr<FMM> fmm[fmm_params] = {
            // pairwise calculation between diff cells regime
            std::make_unique<FMM>(p, theta, 1, 200, invsq.get(), grav.get()),
            // pairwise calculation all in the same cell regime
            std::make_unique<FMM>(p, theta, 200, 1, invsq.get(), grav.get()),
            // little pairwise interaction as possible regime (force FMM)
            std::make_unique<FMM>(p, theta, 1, 0, invsq.get(), grav.get()),
            // normal configuration
            std::make_unique<FMM>(p, theta, 4, 4, invsq.get(), grav.get())
        };

        for (int i = 0; i < fmm_params; i++) {
            //std::cout << i << std::endl;
            TestTriangle(fmm[i].get());
        }
        for (int i = 0; i < fmm_params; i++) {
            TestWellSeparated(fmm[i].get(), 1, 1); // pairwise
            TestWellSeparated(fmm[i].get(), 1, 2); // verifies recursion 
            TestWellSeparated(fmm[i].get(), 1, 10); // FMM likely active
            TestWellSeparated(fmm[i].get(), 10, 10); // usual case
        }
    }
}

}

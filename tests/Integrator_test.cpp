#include <memory>
#include "doctest/doctest.h"
#include "Integrator.hpp"

namespace sim {
 
Grid AnalyticConstAccel(const Grid& g1, double time) {
    Grid grid = g1;
    for (int i = 0; i < grid.GetSize(); i++) {
        grid[i].vel += g1[i].accel * time;
        // use average velocity
        grid[i].pos += (g1[i].vel + grid[i].vel) / 2 * time;
    }
    return grid;
}

// test with constant acceleration scenario
void TestIntegrator(std::unique_ptr<const Integrator> integrator) {
    Grid grid = Grid(Octant({{-10, 10}, {-10, 10}, {-10, 10}}));
    double time = 1;
    int n_steps = int(1e5);
    double step_size = time / n_steps;

    for (int i = 0; i < 10; i++) {
        double x = i * 0.5;
        double v = i * 0.7;
        double a = i * 2;
        Particle par(1, 1);
        par.pos = Vec({x, x, x});
        par.vel = Vec({v, v, v});
        par.accel = Vec({a, a, a});
        grid.AddParticle(par);
    }

    Grid grid_iter = grid; // keep grid for accel reference
   
    for (int i = 0; i < n_steps; i++) {
        grid_iter = integrator->Evolve(grid_iter, step_size);
        for (int i = 0; i < grid_iter.GetSize(); i++) {
            grid_iter[i].accel = grid[i].accel; // const accel
        }
    }
    
    Grid grid_exact = AnalyticConstAccel(grid, time); 
    
    for (int i = 0; i < grid.GetSize(); i++) {
        for (int j = 0; j < Vec::mDim; j++) {
            CHECK(doctest::Approx(grid_exact[i].pos[j]).epsilon(0.05) == 
                    grid_iter[i].pos[j]);
            CHECK(doctest::Approx(grid_exact[i].accel[j]).epsilon(0.05) == 
                    grid_iter[i].accel[j]);
        }
    }
}

TEST_CASE("test Integrator") {
    // not really significant; Just a placeholder.
    // We use evolve(grid, step) functions to manually overwrite the behaviour.
    double dummy_step = 0.001;
    std::unique_ptr<Integrator> euler_int(new Euler(dummy_step));
    std::unique_ptr<Integrator> lf_int(new LeapFrog(dummy_step));
    
    SUBCASE("test GetStep") {
        CHECK(euler_int->GetStep() == doctest::Approx(dummy_step));
    }

    SUBCASE("test Euler") {
        TestIntegrator(std::make_unique<Euler>(dummy_step));
    }

    SUBCASE("test Leapfrog") {
        TestIntegrator(std::make_unique<LeapFrog>(dummy_step));
    }
}

}

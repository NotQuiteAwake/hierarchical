#include "doctest/doctest.h"
#include "Octree.hpp"
#include "Grid.hpp"

namespace sim {

/* test case
   1 1
   -1 -1 -1
   1 2
   -0.7 -0.7 -0.7
   1 3
   -0.4 -0.4 -0.4
   1 4
   -0.1 -0.1 -0.1
   1 5
   0.2 0.2 0.2
*/
Grid AddTestParticles(const Grid& g1) {
    Grid grid = g1;
    for (int i = 0; i < 5; i++) {
        Particle par(1, i + 1);
        double x = -1 + i * 0.3;
        par.pos = Vec({x, x, x});
        grid.AddParticle(par);
    }
    // IO::DumpGrid(grid);
    return grid;
}

void CheckTree(Octree const* node) {
    assert(node);
    Vec coc = Vec();
    double charge = 0;
    for (int soul : node->mSouls) {
        const Particle& par = node->GetParticle(soul);
        coc += par.pos * par.GetCharge();
        charge += par.GetCharge();
        CHECK(node->GetOctant().Within(par.pos));
    }
    coc = coc / charge;
    for (int i = 0; i < Vec::mDim; i++) {
        CHECK(coc[i] == doctest::Approx(node->coc[i]));
    }

    if (node->IsLeaf()) {
        CHECK(node->GetMaxParticles() >= node->mSouls.size());
        int child_cnt = 0;
        for (int i = 0; i < node->mBoxes; i++) {
            if (node->GetChild(i)) child_cnt++;
        }
        CHECK(child_cnt == 0);
        return;
    } else {
        int child_cnt = 0;
        for (int i = 0; i < node->mBoxes; i++) {
            Octree const* child_node = node->GetChild(i);
            if (!child_node) continue;
            child_cnt++;
            CheckTree(child_node);
        }
        CHECK(child_cnt > 0);
    }
    
}

TEST_CASE("test Octree") {
    Octant octant({{-1, 1}, {-1, 1}, {-1, 1}});
    Grid grid(octant); 
    int p = 3;

    grid = AddTestParticles(grid);

    SUBCASE("test M, F initialisation") {
        std::unique_ptr<Octree> root = Octree::BuildTree(grid, 1, p);
        CHECK(root->M.GetRows() == p);
        CHECK(root->M.GetCols() == p);
    }

    SUBCASE("test GetParticle") {
        std::unique_ptr<Octree> root = Octree::BuildTree(grid, 1, p);
        CHECK(root->GetParticle(3).pos[1] == doctest::Approx(-0.1));
        CHECK(root->GetParticle(1).pos[1] == doctest::Approx(-0.7));
    }

    SUBCASE("test maxParticles = 1") {
        std::unique_ptr<Octree> root = Octree::BuildTree(grid, 1, p);
        CheckTree(root.get());
    }

    SUBCASE("test mP = 2") {
        std::unique_ptr<Octree> root = Octree::BuildTree(grid, 2, p);
        CheckTree(root.get());
    }
}

}

#include <sstream>
#include "doctest/doctest.h"
#include "IO.hpp"

namespace sim {
namespace IO {

typedef std::complex<double> cdouble;
static const int p = 3;

bool IsEmpty(std::stringstream& ss) {
    char tmp_char;
    return !bool(ss >> tmp_char);
}

// fill Octree with random nonsense.
void FillOctree(Octree* node) {
    using namespace std::complex_literals;
    if (node->IsLeaf()) {
        for (int i = 0; i <= p; i++) {
            node->M[i][i] = 1.0 * i + 2i;
            node->F[i][i - 1] = 1.0 * 2i;
        }
        return;
    }
    for (int i = 0; i < node->mBoxes; i++) {
        Octree* child_node = node->GetChild(i);
        if (!child_node) continue;
        FillOctree(child_node);
        for (int j = 0; j <= p; j++) {
            node->M[j][j] += child_node->M[j][j];
            node->F[j][j - 1] += child_node->F[j][j - 1];
        }
    }
}

void CheckVec(const Vec& v1, const Vec& v2) {
    for (int i = 0; i < Vec::mDim; i++) {
        CHECK(v1[i] == doctest::Approx(v2[i]));
    }
}

void CheckParticles(const Particle& p1, const Particle& p2) {
    CHECK(p1.GetCharge() == doctest::Approx(p2.GetCharge()));
    CHECK(p1.GetMass() == doctest::Approx(p2.GetMass()));
    CheckVec(p1.pos, p2.pos);
    CheckVec(p1.vel, p2.vel);
    CheckVec(p1.accel, p2.accel);
}


void CheckOctant(const Octant& o1, const Octant& o2) {
    CHECK(o1.IsInitialised() == o2.IsInitialised());
    for (int i = 0; i < Octant::mDim; i++) {
        for (int j = 0; j < 2; j++) {
            CHECK(o1[i][j] == doctest::Approx(o2[i][j]));
        }
    }
}

void CheckMatrix(const Matrix<cdouble>& m1, const Matrix<cdouble>& m2) {
    int nrow = m1.GetRows();
    int ncol = m1.GetCols();
    CHECK(nrow == m2.GetRows());
    CHECK(ncol == m2.GetCols());
    for (int i = 0; i <= nrow; i++) {
        for (int j = -ncol; j <= ncol; j++) {
            CHECK(m1[i][j] == m2[i][j]);
        }
    }
}

// a single DFS, tracking two nodes at the same time
void CheckOctree(Octree* n1, Octree* n2) {
    if (!n1) {
        CHECK(n2 == nullptr);
        return;
    }

    CheckOctant(n1->GetOctant(), n2->GetOctant());
    for (int i = 0; i < n1->mSouls.size(); i++) {
        CHECK(n1->mSouls[i] == n2->mSouls[i]);
        // don't check mOctantSouls; n2 does not load this.
    }

    if (n1->IsLeaf()) {
        CHECK(n2->IsLeaf());
        return;
    } 

    for (int i = 0; i < n1->mBoxes; i++) {
        // child node 1
        Octree* cn1 = n1->GetChild(i);
        Octree* cn2 = n2->GetChild(i);
        CheckOctree(cn1, cn2);
    }
    
}

TEST_CASE("test IO") {
    std::stringstream ss;

    Vec vec({1.1, 2.3, 3.5});
    Octant octant({{-10, 20}, {-15, 25}, {-5, 30}});

    Particle p1(1, 1);
    p1.pos = Vec({1, 2, 3});
    p1.vel = Vec({2, 3, 4});
    p1.accel = Vec({4, 5, 6});
    Particle p2(2, 3);
    p2.pos = Vec({3, 1, 4});

    Grid grid(octant);
    for (int i = 0; i < 10; i++) {
        Particle p(i + 1, i+1); // must have positive masses.
        double x = i;
        p.pos = Vec({x, x, x + 1});
        p.vel = Vec({x - 1, x, x});
        p.accel = Vec({x + 1, x + 1, x});
        grid.AddParticle(p);
    }

    // pre-process our Octree into a fully filled one
    auto root = Octree::BuildTree(grid, 3, p);
    FillOctree(root.get());

    SUBCASE("test Octant round-trip") {
        SUBCASE("initialised octant") {
            DumpOctant(octant, ss);
            Octant new_oct = LoadOctant(ss);
            CheckOctant(octant, new_oct);
            CHECK(IsEmpty(ss));
        }

        SUBCASE("uninitialised octant") {
            Octant uninit_oct;
            DumpOctant(uninit_oct, ss);
            Octant new_oct = LoadOctant(ss);
            CheckOctant(uninit_oct, new_oct);
            CHECK(IsEmpty(ss));
        }
    }

    SUBCASE("test Vec round-trip") {
        DumpVec(vec, ss);
        Vec new_vec = LoadVec(ss);
        CheckVec(vec, new_vec);
        CHECK(IsEmpty(ss));
    }

    SUBCASE("Test Particle round-trip") {
        DumpParticle(p1, ss);
        Particle new_par = LoadParticle(ss);
        CheckParticles(p1, new_par);
        CHECK(IsEmpty(ss));
    }

    SUBCASE("test Grid round-trip") {
        DumpGrid(grid, ss);
        Grid new_grid = LoadGrid(ss);
        for (int i = 0; i < grid.GetSize(); i++) {
            const Particle& p = grid[i];
            const Particle& new_p = new_grid[i];
            CheckParticles(p, new_p);
        }
        CHECK(IsEmpty(ss));
    }

    SUBCASE("test Octree round-trip") {
        DumpOctree(root.get(), ss);
        auto new_root = LoadOctree(grid, ss); 
        for (int i = 0; i < grid.GetSize(); i++) {
            CheckParticles(root->GetParticle(i), new_root->GetParticle(i));
        }
        CheckOctree(root.get(), new_root.get());
        CHECK(IsEmpty(ss));
    }
}

}
}

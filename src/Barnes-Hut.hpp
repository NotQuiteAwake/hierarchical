#ifndef BARNESHUTHEADERDEF
#define BARNESHUTHEADERDEF

#include "Grid.hpp"
#include "Octree.hpp"

namespace sim {

class BarnesHut {
    private:
        const int mP = 1; // order of multipole expansion
        const double mTheta;

        // step 1 of B-H: generate the tree
        Octree BuildTree(const Grid& g1) const;
        bool IsFarAway(const Particle& p1, const Octree* n1) const;
        // acceleration for a single particle
        double GetAccel(const Grid& g1, int index) const;

    public:
        BarnesHut(int p, double theta);
       
        int GetP() const;
        double GetTheta() const;

        // performs the whole routine, setting accel in each particle in grid
        Grid Calculate(const Grid& g1) const;
};

}

#endif

#ifndef BARNESHUTHEADERDEF
#define BARNESHUTHEADERDEF

#include <memory>
#include "Grid.hpp"
#include "Octree.hpp"
#include "Force.hpp"
#include "Interaction.hpp"

namespace sim {

class BarnesHut : public Interaction {
    private:
        const int mP = 1; // TODO: order of multipole expansion
        const double mTheta;

        double SumChildMass(Octree* node) const;
        double SumChildM(Octree* node) const;
        Vec GetNodeCOM(Octree* node) const;

        void ProcessTree(Octree* node) const;
        // step 1 of B-H: generate the tree
        std::unique_ptr<Octree> BuildTree(const Grid& g1) const;

        bool IsFarAway(const Particle& par, const Octree* const node) const;
        // soul needed to assert identity
        Vec GetAccel(int soul, const Octree* const node) const;

    public:
        BarnesHut(int p, double theta,
                const std::shared_ptr<const Force> forceLaw);
       
        int GetP() const;
        double GetTheta() const;

        // performs the whole routine, setting accel in each particle in grid
        Grid Calculate(const Grid& g1) const;
};

}

#endif

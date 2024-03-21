#ifndef BARNESHUTHEADERDEF
#define BARNESHUTHEADERDEF

#include "Grid.hpp"
#include "Octree.hpp"
#include "Force.hpp"
#include "Interaction.hpp"
#include "Kernels.hpp"

namespace sim {

class BarnesHut : public Interaction {
    private:
        const int mP = 1; // TODO: order of multipole expansion
        const double mTheta;
        Kernels* mKernels;

        bool IsFarAway(Octree const* node, const Particle& par) const;
        // soul needed to assert identity
        Particle& EvaluateAccel(
                int soul,
                Particle& par,
                Octree const* node
        ) const;

    public:
        BarnesHut(int p,
                  double theta,
                  Kernels* mKernels,
                  Force const* forceLaw);
       
        int GetP() const;
        double GetTheta() const;

        // performs the whole routine, setting accel in each particle in grid
        Grid Calculate(const Grid& g1) const override;
};

}

#endif

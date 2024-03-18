#ifndef FMMHEADERDEF
#define FMMHEADERDEF

#include <memory>
#include "Octree.hpp"
#include "Kernels.hpp"
#include "Interaction.hpp"

namespace sim {

class FMM : public Interaction {
    private:
        const int mP;
        const double mTheta;
        const int mMaxPerCell;
        const int mMaxPairwiseLimit;
        const std::unique_ptr<const Kernels> mKernels;

        bool MAC(const Octree* node1, const Octree* node2) const;

        void Interact(Octree* node1, Octree* node2, Grid& grid) const;
        // L push-down at this step
        void EvaluateAccel(Octree* node, Grid& grid) const;

    public:
        FMM(int p,
            double theta,
            int maxPerCell,
            int maxPairwiseLimit,
            const std::unique_ptr<const Kernels> kernels,
            const std::unique_ptr<const Force> forceLaw);

        int GetP() const;
        double GetTheta() const;
        
        Grid Calculate(const Grid& g1) const override;
};

}

#endif

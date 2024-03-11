#ifndef FMMHEADERDEF
#define FMMHEADERDEF

#include "Octree.hpp"

namespace sim {

class FMM {
    private:
        const int mP;
        const int mMaxPerCell;
        const double mTheta;
        
        bool MAC(const Octree* t1, const Octree* t2) const;

        // build initial tree with no calculations done.
        Octree BuildTree(const Grid& g1) const;

        // destructive, since we don't care about half-baked octrees.
        void P2M(Octree* leafT, const Grid& g1) const;
        void M2M(Octree* nonleafT) const;
        void M2L(Octree* t1, Octree* t2, const Grid& g1) const;
        void L2L(Octree* t1, Octree* t2) const;
        void L2P(Octree* leafT, const Grid& g1) const;
        
        void CalculateM(Octree* t1) const;
        void Interact(Octree* t1, Octree* t2) const;
        // L push-down at this step?
        void EvaluateAccel(Octree *t1, const Grid& g1) const;

    public:
        FMM(int p, double theta);

        int GetP() const;
        double GetTheta() const;
        
        Grid Calculate(const Grid& g1) const; 
};

}

#endif

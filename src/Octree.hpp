#ifndef OCTREEHEADERDEF
#define OCTREEHEADERDEF

#include <memory>
#include <vector>
#include "Particle.hpp"
#include "Grid.hpp"
#include "Octant.hpp"

namespace sim {

class Octree {
    public:
        static const int mBoxes = 8; // TODO: rework naming since these are
                                     // public
        static const int mDim = 3;

    private:
        int mMaxParticles;
        bool mIsLeaf = true;
        // we are not responsible for the parent's lifetime (no new/delete)
        Octree* mParent; 
        std::unique_ptr<Octree> mChildren[mBoxes];
        const Grid& mGrid;

        Octree(const Octree& otherTree);  // disallow this with unique_ptr
        void Split(); // split the tree into 8


    public:
        Octant octant;
        std::vector<int> mOctantSouls[mBoxes]; // number pointer to particles
        std::vector<int> mSouls;

        std::unique_ptr<double[]> M; // multipoles; depends on expansion order p
        std::unique_ptr<double[]> L; // local expansion

        Vec com; // centre of mass
        double mass;

        Octree(Octree* parent, const Grid& grid, const Octant& new_oct,
                int maxParticles);
        
        // just want to know what child is, no need for ownership semantics
        Octree* GetChild(int index);
        Octree const* GetChild(int index) const;
        Octree* GetParent();

        double GetMaxLength() const;
        bool IsLeaf() const;
        void AddParticle(const Particle& par, int soul);
        Particle GetParticle(int soul) const;
};

}

#endif

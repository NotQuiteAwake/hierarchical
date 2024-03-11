#ifndef OCTREEHEADERDEF
#define OCTREEHEADERDEF

#include <memory>
#include <vector>
#include "Particle.hpp"

namespace sim {

class Octree {
    private:
        static const int mBoxes = 8;
        static const int mDim = 3;

        bool isLeaf;

        // we are not responsible for the parent's lifetime
        // and don't need to new/delete
        const Octree* mParent; 
        std::unique_ptr<Octree> mChildren[mBoxes];

        Octree(const Octree& otherTree);  // disallow this with unique_ptr

    public:
        const double limits[mDim][2];

        std::vector<int> mOctantSouls[mBoxes]; // number pointer to particles
        std::vector<int> mSouls;

        // mLim[Axes::x][0] <= x < mLim[Axes::x][1]
        double COM[mDim]; // centre of mass
        std::unique_ptr<double[]> M; // multipoles; depends on expansion order p
        std::unique_ptr<double[]> L; // local expansion

        Octree(const Octree* parent, const double (&lim)[mDim][2]); 

        // just want to know what child is, no need for ownership semantics
        // return value could be a nullptr
        Octree* GetChild(int index);
        Octree* GetParent();

        bool IsLeaf() const;

        static void GetOctant(const double (&lim)[mDim][2],
                double (&octLim)[mDim][2], int octant);
        void GetOctant(double (&octLim)[mDim][2], int octant);

        void Split(); // split the tree into 8
        // maintain some statistics along the way
        void AddParticle(const Particle& p, int soul);
};

}

#endif

#ifndef OCTREEHEADERDEF
#define OCTREEHEADERDEF

#include <memory>
#include <vector>
#include "Particle.hpp"
#include "Grid.hpp"
#include "Octant.hpp"
#include "Matrix.hpp"

namespace sim {

class Octree {
    public:
        static const int mBoxes = 8; // TODO: rework naming since these are
                                     // public
        static const int mDim = 3;

    private:
        const int mMaxParticles;
        const int mP;
        // we are not responsible for the parent's lifetime (no new/delete)
        Octree* mParent = nullptr; 
        std::unique_ptr<Octree> mChildren[mBoxes];
        const Grid& mGrid;

        bool mLeafCount = 0;

        Octree(const Octree& otherTree);  // disallow this with unique_ptr
        void GenChildNode(int octantNumber);
        void Split(); // split the tree into the 8 octants

    public:
        Octant octant;
        std::vector<int> mOctantSouls[mBoxes]; // number pointer to particles
        std::vector<int> mSouls;

        ComplexMatrix M; // multipoles; depends on expansion order p
        ComplexMatrix F; // local expansion

        // the next two lines took me a full day to diagnose. Without these lines
        // sometimes com and mass are correctly zeroed sometimes they are not.
        // thanks to valgrind and -fsanitize=address
        // boy is C++ initialization not a bloody mess.
        Vec com = Vec(); // centre of mass
        double mass = 0;

        Octree(Octree* parent,
               const Grid& grid,
               const Octant& new_oct,
               int maxParticles,
               int p);
        
        // just want to know what child is, no need for ownership semantics
        Octree* GetChild(int index);
        void SetChild(int index, Octree* node);
        Octree const* GetChild(int index) const;
        Octree* GetParent();
        Octree const* GetParent() const;

        double GetMaxLength() const;
        bool IsLeaf() const;
        void AddParticle(int soul);
        Particle GetParticle(int soul) const;
        int GetMaxParticles() const;
        int GetP() const;

        static std::unique_ptr<Octree> BuildTree(const Grid& grid,
                                                 int maxParticles,
                                                 int p);
};

}

#endif

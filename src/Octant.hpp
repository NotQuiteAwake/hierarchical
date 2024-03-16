#ifndef OCTANTHEADERDEF
#define OCTANTHEADERDEF

#include "Vec.hpp"

namespace sim {

class Octant {
    public:
        const static int mDim = 3;
        double limits[mDim][2];
        Octant(const double (&lim)[mDim][2]);
        Octant();

        bool Within(const Vec& vec) const;
        int GetOctantNumber(const Vec& vec) const;
        static Octant GetOctant(const Octant& octant, int octantNumber);
        double GetMaxLength() const;
        Octant GetOctant(int octantNumber) const;
        double* operator[](int index);
        double const* operator[](int index) const;
        bool operator==(const Octant& otherOctant) const; // TODO: not floating
};

}

#endif

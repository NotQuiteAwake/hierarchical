#ifndef OCTANTHEADERDEF
#define OCTANTHEADERDEF

#include "Vec.hpp"

namespace sim {

class Octant {
    public:
        const static int mDim = 3;
        constexpr static double margin = 2;

    private:
        double mLimits[mDim][2];
        bool mInitialised = false;

        void SetLim(int axis, int side, double val);

    public:
        Octant(const double (&lim)[mDim][2]);
        Octant();

        void Relax(const Vec& vec);

        bool IsInitialised() const;
        bool Within(const Vec& vec) const;
        int GetOctantNumber(const Vec& vec) const;

        static Octant GetOctant(const Octant& octant, int octantNumber);
        double GetMaxLength() const;
        Octant GetOctant(int octantNumber) const;

        double const* operator[](int index) const;
        bool operator==(const Octant& otherOctant) const; // TODO: not floating
};

}

#endif

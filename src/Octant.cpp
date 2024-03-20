#include <cassert>
#include <algorithm>
#include <cstring>
#include "Octant.hpp"

namespace sim {

Octant::Octant(const double (&lim)[mDim][2]) {
    for (int i = 0; i < mDim; i++) {
        for (int j = 0; j < 2; j++) {
            mLimits[i][j] = lim[i][j];
        }
        double diff = mLimits[i][1] - mLimits[i][0];
        assert(diff >= 0);
    }
    mInitialised = true;
}

Octant::Octant() { std::memset(mLimits, 0, sizeof mLimits); }

void Octant::Relax(const Vec& vec) {
    if (!mInitialised) {
        // only use of mInitialised
        for (int i = 0; i < mDim; i++) {
            mLimits[i][0] = vec[i] - margin;
            mLimits[i][1] = vec[i] + margin;
        }
        mInitialised = true;
    }
    else {
        for (int i = 0; i < Vec::mDim; i++) {
            mLimits[i][0] = std::min(vec[i] - margin, mLimits[i][0]);
            mLimits[i][1] = std::max(vec[i] + margin, mLimits[i][1]);
        }
    }
}

bool Octant::IsInitialised() const {
    return mInitialised;
}

bool Octant::Within(const Vec& vec) const {
    for (int i = 0; i < mDim; i++) {
        // NOTE interval convention
        if (!(mLimits[i][0] <= vec[i] && vec[i] < mLimits[i][1])) return false;
    }
    return true;
}

int Octant::GetOctantNumber(const Vec& vec) const {
    int octant_number = 0;
    assert (Within(vec));
    for (int i = 0; i < mDim; i++) {
        double mid = (mLimits[i][0] + mLimits[i][1]) / 2;
        octant_number |= (mid <= vec[i]) << i;
    }
    return octant_number;
}

Octant Octant::GetOctant(const Octant& octant, int octantNumber) {
    double new_lim[mDim][2];
    const double (&mLimits)[mDim][2] = octant.mLimits;
    for (int i = 0; i < mDim; i++) {
        // from low to high bits, octants encode xyz respectively.
        double mid = (mLimits[i][0] + mLimits[i][1]) / 2;
        if (octantNumber & (1 << i)) {
            new_lim[i][0] = mid;
            new_lim[i][1] = mLimits[i][1];
        } else {
            new_lim[i][0] = mLimits[i][0];
            new_lim[i][1] = mid;
        }
    }
    return Octant(new_lim);
}

double Octant::GetMaxLength() const {
    double max_length = 0;
    for (int i = 0; i < mDim; i++) {
        max_length = std::max(max_length, mLimits[i][1] - mLimits[i][0]);
    }
    return max_length;
}

Octant Octant::GetOctant(int octantNumber) const {
    return GetOctant(*this, octantNumber);
}

double const* Octant::operator[](int index) const {
    return mLimits[index];
}

bool Octant::operator==(const Octant& otherOctant) const {
    for (int i = 0; i < Octant::mDim; i++) {
        for (int j = 0; j < 2; j++) {
            if (mLimits[i][j] != otherOctant[i][j]) {
                return false;
            }
        }
    }
    return true;
}


}

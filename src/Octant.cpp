#include <cassert>
#include <algorithm>
#include "Octant.hpp"

namespace sim {

Octant::Octant(const double (&lim)[mDim][2]) {
    for (int i = 0; i < mDim; i++) {
        for (int j = 0; j < 2; j++) {
            limits[i][j] = lim[i][j];
        }
    }
}

Octant::Octant() {}

bool Octant::Within(const Vec& vec) const {
    for (int i = 0; i < mDim; i++) {
        // NOTE interval convention
        if (!(limits[i][0] <= vec[i] && vec[i] < limits[i][1])) return false;
    }
    return true;
}

int Octant::GetOctantNumber(const Vec& vec) const {
    int octant_number = 0;
    assert (Within(vec));
    for (int i = 0; i < mDim; i++) {
        double mid = (limits[i][0] + limits[i][1]) / 2;
        octant_number |= (mid <= vec[i]) << i;
    }
    return octant_number;
}

Octant Octant::GetOctant(const Octant& octant, int octantNumber) {
    double new_lim[mDim][2];
    const double (&limits)[mDim][2] = octant.limits;
    for (int i = 0; i < mDim; i++) {
        // from low to high bits, octants encode xyz respectively.
        double mid = (limits[i][0] + limits[i][1]) / 2;
        if (octantNumber & (1 << i)) {
            new_lim[i][0] = mid;
            new_lim[i][1] = limits[i][1];
        } else {
            new_lim[i][0] = limits[i][0];
            new_lim[i][1] = mid;
        }
    }
    return Octant(new_lim);
}

double Octant::GetMaxLength() const {
    double lmax = 0;
    for (int i = 0; i < mDim; i++) {
        lmax = std::max(lmax, limits[i][1] - limits[i][0]);
    }
    return lmax;
}

Octant Octant::GetOctant(int octantNumber) const {
    return GetOctant(*this, octantNumber);
}

double* Octant::operator[](int index) {
    return limits[index];
}

double const* Octant::operator[](int index) const {
    return limits[index];
}


bool Octant::operator==(const Octant& otherOctant) const {
    for (int i = 0; i < Octant::mDim; i++) {
        for (int j = 0; j < 2; j++) {
            if (limits[i][j] != otherOctant[i][j]) {
                return false;
            }
        }
    }
    return true;
}


}

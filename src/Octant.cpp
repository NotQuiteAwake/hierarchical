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

int Octant::GetOctantNumber(Vec v1) const {
    int octant_number = 0;
    for (int i = 0; i < mDim; i++) {
        double comp = v1[i];
        assert(limits[i][0] < comp && comp < limits[i][1]);
        double mid = (limits[i][0] + limits[i][1]) / 2;
        octant_number |= (comp > mid) << i;
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

}

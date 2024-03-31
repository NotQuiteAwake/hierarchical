/**
 * @file
 * @brief Data structure representing a cuboid box in space.
 */

#include <cassert>
#include <algorithm>
#include <cstring>
#include "Octant.hpp"

namespace sim {

/**
 * @class
 * @brief Data structure representing a cuboid box in space.
 *
 * Octant has capabilities to obtain octant numbers, and find which octant a
 * vector belongs to, therefore its name.
 */

/**
 * @brief Instantiate an octant from limits directly
 *
 * Octants intantiated this way are considered "initialised."
 *
 * @param[in] lim limits by dimension; second index corresponds to min and max.
 */
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

/**
 * @brief Instantiate an Octant without initialising it.
 */
Octant::Octant() { std::memset(mLimits, 0, sizeof mLimits); }

/**
 * @brief Ask octant to adjust its boundaries to accommodate new vector.
 *
 * @param[in] vec displacement of a vector
 */
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

/**
 * @brief Test if a vector is within the Octant boundary.
 */
bool Octant::Within(const Vec& vec) const {
    for (int i = 0; i < mDim; i++) {
        // NOTE interval convention
        if (!(mLimits[i][0] <= vec[i] && vec[i] < mLimits[i][1])) return false;
    }
    return true;
}

/**
 * @brief Get which octant a vector lies within the current Octant
 *
 * The octant number is a three-bit mask. The i-th bit is 1 if vec[i] lies above
 * the midpoint of Octant in that dimension; Else 0.
 *
 * @return Number from 0 to 7 encoding an octant.
 */
int Octant::GetOctantNumber(const Vec& vec) const {
    int octant_number = 0;
    assert (Within(vec));
    for (int i = 0; i < mDim; i++) {
        double mid = (mLimits[i][0] + mLimits[i][1]) / 2;
        octant_number |= (mid <= vec[i]) << i;
    }
    return octant_number;
}

/**
 * @brief Get one octant of the provided octant, by its octant number
 *
 * This typically called after GetOctantNumber(). See documentation of that
 * function for how octant number is defined.
 *
 * @param[in] octant Octant in which to find the sub-octants.
 * @param[in] octantNumber A number from 0 to 7 showing which sub-octant it is.
 * @return The sub-octant represented by octantNumber.
 */
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

/**
 * @brief Get one octant of the current Octant, by its octant number
 *
 * This typically called after GetOctantNumber(). See documentation of that
 * function for how octant number is defined.
 *
 * @param[in] octantNumber A number from 0 to 7 showing which sub-octant it is.
 * @return The sub-octant represented by octantNumber.
 */
Octant Octant::GetOctant(int octantNumber) const {
    return GetOctant(*this, octantNumber);
}

/**
 * @brief Get length of the longest side of Octant.
 */
double Octant::GetMaxLength() const {
    double max_length = 0;
    for (int i = 0; i < mDim; i++) {
        max_length = std::max(max_length, mLimits[i][1] - mLimits[i][0]);
    }
    return max_length;
}

/**
 * @brief Get the limits of the current Octant in the specified dimension.
 *
 * Do not allow modification, for this will break internal tracking of whether
 * octant is initialised or not.
 *
 * @param[in] index Dimension to access
 * @return array of length 2, with lower and upper limits.
 */
double const* Octant::operator[](int index) const {
    return mLimits[index];
}

/**
 * @brief Two octants are equal if each of their limits agree.
 *
 * This operation is not safe against floating point errors and is only meant to
 * check whether two octants are identical copies of each other.
 */
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

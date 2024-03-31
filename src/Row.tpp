#ifndef ROWHEADERDEF_TPP
#define ROWHEADERDEF_TPP

/**
 * @file
 * @brief Data structure representing a row in the Matrix class.
 */

#include "Row.hpp"

#include <cmath>

namespace sim {

/**
 *
 * @class
 * @brief Data structure representing a row in the Matrix class.
 *
 * The index runs from -ncol to ncol inclusive for the representation of
 * multipole expansion coefficients. See Matrix class documentation for more
 * details.
 */

/**
 * @brief Initialise Row to have property ncols.
 *
 * @tparam T a numerical type
 * @param[in] ncols Allows for access to elements -ncol to ncol inclusive.
 */
template<typename T> void Row<T>::Resize(int ncols) {
    // must work from -ncols to ncols, inclusive
    mRow.resize(ncols * 2 + 1);
    mColCnt = ncols;
}

/**
 * @brief Resize Row to have property ncols.
 *
 * @tparam T a numerical type
 * @param[in] ncols Allows for access to elements -ncol to ncol inclusive.
 */
template<typename T> Row<T>::Row(int ncols) {
    Resize(ncols);
}

template<typename T> Row<T>::Row(): mColCnt(0) {};

template<typename T> int Row<T>::GetSize() const {
    return mColCnt;
}

/**
 * @brief Access element at colIndex.
 *
 * Note that 0 will be returned for an out-of-bounds access - this is because
 * the mathematics we are trying to represent, solid harmonics of order p goes
 * to zero for m > n, or colIndex > ncols.
 *
 * @return element value if colIndex <= ncols, otherwise 0.
 */
template<typename T> T Row<T>::operator[](int colIndex) const {
    // for this we implicitly require type T to have a conversion operator, and
    // this forbids us from return by const ref (of this local var)
    if (std::abs(colIndex) > mColCnt) return T(0.0);
    // shift all indices to >= 0
    return mRow[colIndex + mColCnt];
}

template<typename T> T& Row<T>::operator[](int colIndex) {
    assert(std::abs(colIndex) <= mColCnt);
    return mRow[colIndex + mColCnt];
}

/**
 * @brief Perform elementwise addition
 */
template<typename T> Row<T> Row<T>::operator+(const Row<T>& otherRow) const {
    assert(mColCnt == otherRow.GetSize());
    Row new_row(mColCnt);
    
    for (int i = -mColCnt; i <= mColCnt; i++) {
        new_row[i] = (*this)[i] + otherRow[i];
    }

    return new_row;
}

/**
 * @brief Elementwise +=
 */
template<typename T> Row<T>& Row<T>::operator+=(const Row<T>& otherRow) {
    assert(mColCnt == otherRow.GetSize());
    for (int i = -mColCnt; i <= mColCnt; i++) {
        (*this)[i] += otherRow[i];
    }
    return *this;
}

}

#endif

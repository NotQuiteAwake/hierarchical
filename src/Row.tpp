#ifndef ROWHEADERDEF_TPP
#define ROWHEADERDEF_TPP

#include "Row.hpp"

#include <cmath>

namespace sim {

template<typename T> void Row<T>::Resize(int ncols) {
    // must work from -ncols to ncols, inclusive
    mRow.resize(ncols * 2 + 1);
    mColCnt = ncols;
}

template<typename T> Row<T>::Row(int ncols) {
    Resize(ncols);
}

template<typename T> Row<T>::Row(): mColCnt(0) {};

template<typename T> int Row<T>::GetSize() const {
    return mColCnt;
}

template<typename T> T Row<T>::operator[](int colIndex) const {
    // IMPORTANT out-of-bounds returns 0. Mathematically these elements go to 0
    // as the spherical harmonics = 0 for such n, m.
    // for this we also implicitly require type T to have a conversion operator.
    // and this forbids us from return by const ref (of this local var)
    if (std::abs(colIndex) > mColCnt) return T(0.0);
    // shift all indices to >= 0
    return mRow[colIndex + mColCnt];
}

template<typename T> T& Row<T>::operator[](int colIndex) {
    return mRow[colIndex + mColCnt];
}

template<typename T> Row<T> Row<T>::operator+(const Row<T>& otherRow) const {
    assert(mColCnt == otherRow.GetSize());
    Row new_row(mColCnt);
    
    for (int i = -mColCnt; i <= mColCnt; i++) {
        new_row[i] = (*this)[i] + otherRow[i];
    }

    return new_row;
}

template<typename T> Row<T>& Row<T>::operator+=(const Row<T>& otherRow) {
    assert(mColCnt == otherRow.GetSize());
    for (int i = -mColCnt; i <= mColCnt; i++) {
        (*this)[i] += otherRow[i];
    }
    return *this;
}

}

#endif

#ifndef MATRIXHEADERDEF_TPP
#define MATRIXHEADERDEF_TPP

#include "Matrix.hpp"

namespace sim {

// pre-allocate
template<typename T> Matrix<T>::Matrix(int nrows, int ncols) {
    Resize(nrows, ncols);
}

template<typename T> void Matrix<T>::Resize(int nrows, int ncols) {
    // must work for 0 to nrows inclusive
    mRows.resize(nrows + 1);
    for (Row<T>& row : mRows) {
        row.Resize(ncols);
    }
    mRowCnt = nrows;
    mColCnt = ncols;
}

template<typename T> const Row<T>& Matrix<T>::operator[](int rowIndex) const {
    return mRows[rowIndex];
}

template<typename T> Row<T>& Matrix<T>::operator[](int rowIndex) {
    return mRows[rowIndex];
}

template<typename T> Matrix<T>& Matrix<T>::operator+=(
        const Matrix<T>& otherMatrix
        ) {
    assert(mRowCnt == otherMatrix.GetRows());
    assert(mColCnt == otherMatrix.GetCols());
    for (int i = 0; i <= mRowCnt; i++) {
        mRows[i] += otherMatrix[i];
    }
    return *this;
}

template<typename T> int Matrix<T>::GetRows() const {
    return mRowCnt;
}

template<typename T> int Matrix<T>::GetCols() const {
    return mColCnt;
}

}

#endif

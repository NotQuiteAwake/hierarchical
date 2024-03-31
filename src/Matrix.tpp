#ifndef MATRIXHEADERDEF_TPP
#define MATRIXHEADERDEF_TPP

/**
 * @file
 * @brief Implementation of the Matrix class for multipole coefficient storage.
 */

#include <cassert>
#include "Matrix.hpp"

namespace sim {

/**
 * @class
 * @brief Implementation of the Matrix class for multipole coefficient storage.
 *
 * We note that nrows and ncols have special meanings: One can access row
 * indices 0 to nrows inclusive, and column indices -ncols ... +ncols inclusive.
 *
 * This is for compatibility with multipole expansion coefficients, which do
 * have this index range.
 */

/**
 * @brief Initialise a Matrix with parameters nrows and ncols.
 *
 * For nrows and ncols see class documentation.
 *
 * @tparam T A numerical type (double, std::complex<double>, ...)
 */
template<typename T> Matrix<T>::Matrix(int nrows, int ncols) {
    Resize(nrows, ncols);
}

/**
 * @brief Resizwe the matrix to work with parameters nrows and ncols
 *
 * For nrows and ncols see class documentation.
 */
template<typename T> void Matrix<T>::Resize(int nrows, int ncols) {
    assert(nrows >= 0);
    assert(ncols >= 0);

    mRows.resize(nrows + 1);
    for (Row<T>& row : mRows) {
        row.Resize(ncols);
    }
    mRowCnt = nrows;
    mColCnt = ncols;
}

// const& as user don't need to deal with Row, rather just its elements
// and the Row operator[] returns by value anyway
template<typename T> const Row<T>& Matrix<T>::operator[](int rowIndex) const {
    assert(0 <= rowIndex && rowIndex <= mRowCnt);
    return mRows[rowIndex];
}

template<typename T> Row<T>& Matrix<T>::operator[](int rowIndex) {
    assert(0 <= rowIndex && rowIndex <= mRowCnt);
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

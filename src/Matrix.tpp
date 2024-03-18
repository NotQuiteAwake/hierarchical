#ifndef MATRIXHEADERDEF_TPP
#define MATRIXHEADERDEF_TPP

#include "Matrix.hpp"

namespace sim {

template<typename T> void Row<T>::Resize(int ncols) {
    // must work from -ncols to ncols, inclusive
    mRow.resize(ncols * 2 + 1);
    mColCnt = ncols;
}

template<typename T> Row<T>::Row() {}

template<typename T> Row<T>::Row(int ncols) {
    Resize(ncols);
}

template<typename T> const T& Row<T>::operator[](int colIndex) const {
    // shift all indices to >= 0
    return mRow[colIndex + mColCnt];
}

template<typename T> T& Row<T>::operator[](int colIndex) {
    return mRow[colIndex + mColCnt];
}

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

template<typename T> int Matrix<T>::GetRows() const {
    return mRowCnt;
}

template<typename T> int Matrix<T>::GetCols() const {
    return mColCnt;
}

}

#endif

#ifndef MATRIXHEADERDEF_HPP
#define MATRIXHEADERDEF_HPP

#include <vector>
#include <complex>
#include "Row.hpp"

namespace sim {

template<typename T> class Matrix {
    private: 
        std::vector<Row<T>> mRows;
        int mRowCnt;
        int mColCnt;

    public:
        Matrix(int nrows, int ncols);
        void Resize(int nrows, int ncols);

        const Row<T>& operator[](int rowIndex) const;
        Row<T>& operator[](int rowIndex);
        Matrix<T>& operator+=(const Matrix<T>& otherMatrix);

        int GetRows() const;
        int GetCols() const;
};

typedef Matrix<std::complex<double>> ComplexMatrix;
typedef Matrix<double> DoubleMatrix;

}

#include "Matrix.tpp"

#endif

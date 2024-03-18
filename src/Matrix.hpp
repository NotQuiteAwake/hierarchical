#ifndef MATRIXHEADERDEF_HPP
#define MATRIXHEADERDEF_HPP

#include <vector>
#include <complex>

namespace sim {

template<typename T> class Matrix;

template<typename T> class Row {
    private:
        std::vector<T> mRow;
        int mColCnt;

        void Resize(int ncols);

    public:
        // alternatively employ Matrix* mat
        Row(int ncols);
        Row();

        const T& operator[](int colIndex) const;
        T& operator[](int colIndex);

    friend class Matrix<T>;
};


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

        int GetRows() const;
        int GetCols() const;
};

typedef Matrix<std::complex<double>> ComplexMatrix;
typedef Matrix<double> DoubleMatrix;

}

#include "Matrix.tpp"

#endif

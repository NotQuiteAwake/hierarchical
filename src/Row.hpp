#ifndef ROWHEADERDEF_HPP
#define ROWHEADERDEF_HPP

#include <vector>

namespace sim {

template<typename T> class Matrix; // forward declaration

template<typename T> class Row {
    private:
        std::vector<T> mRow;
        int mColCnt;

        void Resize(int ncols);

    public:
        // alternatively employ Matrix* mat
        Row(int ncols);
        Row();

        // for all purposes this would appear like a normal vector
        // but its actual range is modified.
        int GetSize() const;

        T operator[](int colIndex) const;
        T& operator[](int colIndex);
        Row<T> operator+(const Row<T>& otherRow) const;
        Row<T> operator-(const Row<T>& otherRow) const;
        Row<T>& operator+=(const Row<T>& otherRow);

    friend class Matrix<T>;
};

}

#include "Row.tpp"

#endif

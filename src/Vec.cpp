#include "Vec.hpp"
#include <cmath>

namespace sim {

Vec::Vec() {
    mCoords.reset(new double[3]);
}

Vec::Vec(const Vec& otherVec) {
    mCoords.reset(new double[3]);
    for (int i = 0; i < mDim; i++) {
        mCoords[i] = otherVec[i];
    } 
}

Vec::Vec(const double (&coords)[mDim]) {
    mCoords.reset(new double[3]);
    for (int i = 0; i < mDim; i++) {
        mCoords[i] = coords[i];
    }
}

double& Vec::operator[](int index) {
    return mCoords[index];
}

double Vec::operator[](int index) const {
    return mCoords[index];
}

Vec Vec::operator+() const {
    return Vec(*this);
}

Vec Vec::operator-() const {
    Vec new_vec = Vec();
    for (int i = 0; i < mDim; i++) {
        new_vec[i] = -mCoords[i];
    }
    return new_vec;
}

Vec Vec::operator+(const Vec& otherVec) const {
    Vec new_vec = Vec();
    for (int i = 0; i < Vec::mDim; i++) {
        new_vec[i] = mCoords[i] + otherVec[i];
    }
    return new_vec;
}

Vec Vec::operator-(const Vec& otherVec) const {
    return *this + (-otherVec);
}

Vec Vec::operator*(double factor) const {
    Vec new_vec = Vec();
    for (int i = 0; i < Vec::mDim; i++) {
        new_vec[i] = mCoords[i] * factor;
    }
    return new_vec;
}

Vec Vec::operator/(double factor) const {
    Vec new_vec = Vec();
    for (int i = 0; i < Vec::mDim; i++) {
        new_vec[i] = mCoords[i] / factor;
    }
    return new_vec;
}

Vec& Vec::operator=(const Vec& otherVec) {
    for (int i = 0; i < Vec::mDim; i++) {
        mCoords[i] = otherVec[i];
    } 
    return *this;
}

// TODO: test this
Vec& Vec::operator+=(const Vec& otherVec) {
    for (int i = 0; i < Vec::mDim; i++) {
        mCoords[i] += otherVec[i];
    }
    return *this;
}

double Vec::CalculateNorm() const {
    double norm = 0;
    for (int i = 0; i < Vec::mDim; i++) {
        norm += mCoords[i] * mCoords[i];
    }
    norm = sqrt(norm);
    return norm;
}

double DotProduct(const Vec& v1, const Vec& v2) {
    double product = 0;
    for (int i = 0; i < Vec::mDim; i++) {
        product += v1[i] * v2[i];  
    }
    return product;
}

Vec CrossProduct(const Vec& v1, const Vec& v2) {
    Vec product = Vec();
    for (int i = 0; i < Vec::mDim; i++) {
        int j = (i + 1) % Vec::mDim;
        int k = (i + 2) % Vec::mDim;
        product[i] = v1[j] * v2[k] - v2[j] * v1[k];
    }
    return product;
}

std::ostream& operator<<(std::ostream& output, const Vec& v1) {
    output << "(" << v1[0];
    for (int i = 1; i < Vec::mDim; i++) {
        output << ", " << v1[i];
    }
    output << ")";
    return output;
}

}

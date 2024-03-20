#include <cmath>
#include <iostream>
#include <cstring>
#include "Vec.hpp"

namespace sim {

static constexpr double PI() { return std::atan(1)*4; }

Vec::Vec() {
    std::memset(mCoords, 0, sizeof mCoords);
}

Vec::Vec(const Vec& otherVec) {
    for (int i = 0; i < mDim; i++) {
        mCoords[i] = otherVec[i];
    } 
}

Vec::Vec(const double (&coords)[mDim]) {
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

Vec& Vec::operator+=(const Vec& otherVec) {
    for (int i = 0; i < Vec::mDim; i++) {
        mCoords[i] += otherVec[i];
    }
    return *this;
}

Vec& Vec::operator-=(const Vec& otherVec) {
    for (int i = 0; i < Vec::mDim; i++) {
        mCoords[i] -= otherVec[i];
    }
    return *this;
}

// exact identity. Don't care about floating point error.
bool Vec::operator==(const Vec& otherVec) const {
    for (int i = 0; i < Vec::mDim; i++) {
        if (mCoords[i] != otherVec[i]) return false;
    }
    return true;
}

double Vec::GetNorm() const {
    double norm = 0;
    for (int i = 0; i < Vec::mDim; i++) {
        norm += mCoords[i] * mCoords[i];
    }
    norm = sqrt(norm);
    return norm;
}

// this is not useful for all cases where coords change,
// so don't maintain a value on change with operator[]
double Vec::GetTheta() const {
    double norm = GetNorm();
    if (!norm) return 0;
    return PI() / 2 - std::asin(mCoords[2] / GetNorm());
}

// conform to boost's spherical_harmonics convention
double Vec::GetPhi() const {
    const double& x = mCoords[0];
    const double& y = mCoords[1];
    if (!x) {
        if (y == 0) return 0;
        if (y > 0) return PI() / 2;
        else return PI() * 3 / 2;
    } else {
        double phi = std::atan(y / x);
        if (y > 0 && x < 0) {
            // direction is wrong in 2nd quadrant
            phi += PI();
        } else if (y < 0 && x < 0) {
            // direction is wrong in 3rd quadrant
            phi += PI();
        } else if (y < 0 && x > 0) {
            // bring 4th quadrant to within 0, 2pi
            phi += 2 * PI();
        }
        return phi;
    }
}

Vec Vec::FromSpherical(double r, double theta, double phi) {
    double x = r * std::sin(theta) * std::cos(phi);
    double y = r * std::sin(theta) * std::sin(phi);
    double z = r * std::cos(theta);
    return Vec({x, y, z});
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

#ifndef VECHEADERDEF
#define VECHEADERDEF

#include <memory>
#include <iostream>

namespace sim {

class Vec {
    private:
        static const int mDim = 3; // TODO: should this be public?
        std::unique_ptr<double[]> mCoords;

    public:
        Vec();
        Vec(const Vec& otherVec);
        Vec(const double (&coords)[mDim]);

        double& operator[](int index);
        double operator[](int index) const;

        Vec operator+() const;
        Vec operator-() const;

        Vec operator+(const Vec& otherVec) const;
        Vec operator-(const Vec& otherVec) const;

        Vec operator*(double factor) const;
        Vec operator/(double factor) const;

        Vec& operator=(const Vec& otherVec);
        Vec& operator+=(const Vec& otherVec);

        double CalculateNorm() const;
        
        friend double DotProduct(const Vec& v1, const Vec& v2);
        friend Vec CrossProduct(const Vec& v1, const Vec& v2);

        friend std::ostream& operator<<(std::ostream& output,
                                const Vec& otherVec);
};

}

#endif

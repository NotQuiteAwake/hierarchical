#ifndef VECHEADERDEF
#define VECHEADERDEF

#include <memory>

namespace sim {

class Vec {
    private:
        static const int mDim = 3;
        std::unique_ptr<double[]> mCoords;

    public:
        Vec();
        Vec(const Vec& otherVec);
        Vec(const double (&coords)[mDim]);

        Vec operator+() const;
        Vec operator-() const;

        Vec operator+(const Vec& otherVec) const;
        Vec operator-(const Vec& otherVec) const;

        Vec operator*(double factor) const;

        Vec& operator=(const Vec& otherVec);

        double& operator[](int index);
        double operator[](int index) const;

        double CalculateNorm() const;
        
        friend double DotProduct(const Vec& v1, const Vec& v2);
        friend Vec CrossProduct(const Vec& v1, const Vec& v2);
};

}

#endif

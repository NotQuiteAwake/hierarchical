#ifndef VECHEADERDEF
#define VECHEADERDEF

#include <iostream>

namespace sim {

class Vec {
    public:
        static const int mDim = 3;

    private:
        double mCoords[mDim];

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
        Vec& operator-=(const Vec& otherVec);

        bool operator==(const Vec& otherVec) const;

        double GetNorm() const;
        // spherical polar
        double GetTheta() const;
        double GetPhi() const;
        
        static Vec FromSpherical(double r, double theta, double phi);

        friend double DotProduct(const Vec& v1, const Vec& v2);
        friend Vec CrossProduct(const Vec& v1, const Vec& v2);

        friend std::ostream& operator<<(std::ostream& output,
                                const Vec& otherVec);
};

}

#endif

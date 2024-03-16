#ifndef PARTICLEHEADERDEF
#define PARTICLEHEADERDEF

#include "Vec.hpp"

namespace sim {

class Particle {
    private:
        double mMass;
        double mCharge;

    public:
        Vec pos;
        Vec vel;
        Vec accel;

        Particle(double mass, double charge);

        double GetMass() const;
        double GetCharge() const;
};

}

#endif

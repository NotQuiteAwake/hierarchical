#ifndef PARTICLEHEADERDEF
#define PARTICLEHEADERDEF

#include "Vec.hpp"

namespace sim {

class Particle {
    private:
        double mMass;
        double mCharge;

    public:
        Vec pos = Vec();
        Vec vel = Vec();
        Vec accel = Vec();

        Particle(double mass, double charge);

        double GetMass() const;
        double GetCharge() const;
};

}

#endif

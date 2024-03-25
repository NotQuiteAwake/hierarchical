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
        double pot = 0;

        Particle(double mass, double charge);

        double GetMass() const;
        double GetCharge() const;
        double GetPE() const;
        double GetKE() const;
        Vec GetP() const;
};

}

#endif

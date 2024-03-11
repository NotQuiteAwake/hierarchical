#ifndef PARTICLEHEADERDEF
#define PARTICLEHEADERDEF

#include "Vec.hpp"

namespace sim {

class Particle {
    private:
        double mMass;
        double mCharge;
        Vec mPos;
        Vec mVel;
        Vec mAccel;

    public:
        double GetMass() const;
        double GetCharge() const;

        Vec GetPos() const;
        Vec GetVel() const;
        Vec GetAccel() const;

        void SetPos(const Vec& pos);
        void SetVel(const Vec& vel);
        void SetAccel(const Vec& accel);
};

}

#endif

#ifndef GRIDHEADERDEF
#define GRIDHEADERDEF

#include <memory>
#include "Particle.hpp"

namespace sim {

class Grid {
    private:
        const int mSize;
        const double mLength;
        std::unique_ptr<Particle> mParticles;

    public:
        Grid(int size, double length);
        Grid(const Grid& otherGrid);
        Grid& operator=(const Grid& otherGrid);

        int GetSize() const;
        double GetLength() const;

        Particle& operator[](int index);
        Particle operator[](int index) const;
};

}

#endif


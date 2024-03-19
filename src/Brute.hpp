#ifndef BRUTEHEADERDEF
#define BRUTEHEADERDEF

#include "Grid.hpp"
#include "Interaction.hpp"

namespace sim {

class Brute : public Interaction {
    public:
        Brute(Force const* forceLaw);
        Grid Calculate(const Grid& g1) const override;
};

}

#endif

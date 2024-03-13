#ifndef BRUTEHEADERDEF
#define BRUTEHEADERDEF

#include "Grid.hpp"
#include "Interaction.hpp"

namespace sim {

class Brute : public Interaction {
    public:
        Grid Calculate(const Grid& g1) const override;
};

}

#endif

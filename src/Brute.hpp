#ifndef BRUTEHEADERDEF
#define BRUTEHEADERDEF

#include "Grid.hpp"
#include "Interaction.hpp"

namespace sim {

class Brute : public Interaction {
    public:
        Brute(std::shared_ptr<Force> forceLaw);
        Grid Calculate(const Grid& g1) const override;
};

}

#endif

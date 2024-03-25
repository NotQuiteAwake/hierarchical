#ifndef EXPERIMENTSHEADERDEF
#define EXPERIMENTSHEADERDEF

#include <string>
#include <cmath>
#include "Interaction.hpp"

namespace sim {
namespace exp {

namespace {

static const std::string dataDir = "./data/";
static constexpr double PI() { return std::atan(1) * 4; }

void Err(const std::string& errMessage);
bool CheckFile(const std::string& fileName, bool expect);
int TimeUniformRun(Interaction const* interaction, int n, int seed);

}

void TimeComplexity();
void ExpansionOrderComplexity();
void ThetaComplexity();

namespace {
void SimulationHelper(const std::string dump_dir,
                      const Grid& grid,
                      const double scale
        );
}

void ColdStartSim();
void ThinDiskSim();
void GalaxySim();
void TwoGalaxiesSim();

}
}

#endif

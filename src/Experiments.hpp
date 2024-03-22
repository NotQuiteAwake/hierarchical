#include <string>
#include "Interaction.hpp"

namespace sim {
namespace exp {

namespace {
static const std::string dataDir = "./data/";

void Log(const std::string& logMessage);
void Err(const std::string& errMessage);

bool CheckFile(const std::string& fileName, bool expect);

int TimeUniformRun(Interaction const* interaction, int n, int seed);
}


void TimeComplexity();

}
}


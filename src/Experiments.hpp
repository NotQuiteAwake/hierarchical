#ifndef EXPERIMENTSHEADERDEF
#define EXPERIMENTSHEADERDEF

#include <string>
#include <cmath>
#include <memory>
#include <fstream>
#include "Interaction.hpp"
#include "InvSqKernels.hpp"
#include "Integrator.hpp"

namespace sim {
namespace exp {

namespace {

constexpr double PI() { return std::atan(1) * 4; }

void Err(const std::string& errMessage);
bool CheckFile(const std::string& fileName, bool expect);
bool StreamSetup(std::ofstream& stream);


const struct GenParams {
    int n = 1000;
    int p = 3;
    double theta = 0.5;
    int repeats = 10;
    int maxPerCell = 5;
    int maxPairwiseLimit = 3;
    double G = -1;
    double step = 0.01;
    double evo_time = 10;
    double scale = 10;

    static const int intTypes = 3;
    std::string intNames[intTypes]{"brute", "bh", "fmm"};
    enum intVals {brute, bh, fmm}; // a bijection...

} genDefParams;

struct FileParams {
    inline static const std::string dataDir = "./data/";
    std::string runDir = "";
    std::string comp = "complexity.out";
    std::string dump = "dump/";

    FileParams(const std::string runDir);
    std::string GetCompName() const;
    std::string GetDumpName() const;
    std::string GetRunDirName() const;
};

bool CompStreamSetup(std::ofstream& stream,
                     FileParams fileParams);

class Benchmarker {
    private:
        inline static auto dummyStream = std::ofstream();
        int TimeUniformRun(Interaction const* interaction,
                int n,
                int seed,
                std::ofstream& stream) const;
    public:
        std::unique_ptr<LeapFrog> integrator;
        std::unique_ptr<InvSqForce> force;
        std::unique_ptr<InvSqKernels> kernels;
        std::unique_ptr<Interaction> interactions[genDefParams.intTypes];
        std::ofstream& stream;
        GenParams genParams;
        FileParams fileParams;
        Benchmarker(GenParams genParams,
                    FileParams fileParams,
                    std::ofstream& stream = dummyStream
                );
        void Run(double param, bool skipBrute) const;
        void Evo(const Grid& grid) const;
};

}

void NComplexity();
void PComplexity();
void ThetaComplexity();

void ColdStartSim();
void ThinDiskSim();
void GalaxySim();
void TwoGalaxiesSim();

}
}

#endif

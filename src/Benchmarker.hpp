#ifndef BENCHMARKERHEADERDEF
#define BENCHMARKERHEADERDEF

#include <string>
#include <memory>

#include "Interaction.hpp"
#include "Integrator.hpp"
#include "InvSqKernels.hpp"

namespace sim {

namespace exp {

const struct GenParams {
    int n = 1000;               // number of particles
    int p = 3;                  // order of multipole expansion
    double theta = 0.5;         // opening angle
    int repeats = 10;           // number of repeats
    int maxPerCell = 5;         // max number of particles per cell (FMM)
    int maxPairwiseLimit = 3;   // limit below which brute-force allowed (FMM)
    double G = -1;              // force coupling constant
    double step = 0.01;         // step size (integrator)
    double evo_time = 10;       // evolution time
    double scale = 10;          // typical scale in simulation

    static const int intTypes = 3; // types of interactions,
    std::string intNames[intTypes]{"brute", "bh", "fmm"}; // their names,
    enum intVals {brute, bh, fmm}; // and how the name maps to numbers...

} genDefParams;

struct FileParams {
    inline static const std::string dataDir = "./data/";
    std::string runDir = "";                // folder to run this sim in
    std::string comp = "complexity.out";    // complexity log file name
    std::string dump = "dump/";             // full grid dump relative path

    FileParams(const std::string runDir);
    std::string GetCompName() const;
    std::string GetDumpName() const;
    std::string GetRunDirName() const;
};

class Benchmarker {
    private:
        int TimeUniformRun(Interaction const* interaction,
                int n,
                int seed,
                std::ofstream& stream) const;
    public:
        std::unique_ptr<LeapFrog> integrator;
        std::unique_ptr<InvSqForce> force;
        std::unique_ptr<InvSqKernels> kernels;
        std::unique_ptr<Interaction> interactions[genDefParams.intTypes];
        std::ostream& stream;
        GenParams genParams;
        FileParams fileParams;
        Benchmarker(GenParams genParams,
                    FileParams fileParams,
                    std::ostream& stream = std::cout
                );
        void Run(double param, bool skipBrute) const;
        void Evo(const Grid& grid) const;
};

}

}

#endif

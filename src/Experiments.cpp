#include <string>
#include <cmath>
#include <memory>
#include <chrono>
#include <fstream>

#include "Experiments.hpp"

#include "Interaction.hpp"
#include "Brute.hpp"
#include "Barnes-Hut.hpp"
#include "FMM.hpp"
#include "InvSqKernels.hpp"
#include "Distribution.hpp"
#include "IO.hpp"

namespace sim {
namespace exp {

namespace {
void Log(const std::string& logMessage) {
    std::clog << "(I) " << logMessage << std::endl;
}

void Err(const std::string& errMessage) {
    std::cerr << "(E) " << errMessage << std::endl;
}

bool CheckFile(const std::string& fileName, bool expect) {
    bool exists = IO::FileExists(fileName);

    if (expect && !exists) {
        std::string msg = "The directory/file " + fileName;
        msg += " does not exist. The called function may refuse to run.";
        Err(msg);
    } else if (!expect && exists) {
        std::string msg = "The directory/file " + fileName;
        msg += " already exists. The called function may refuse to run.";
        Err(msg);
    }
    return exists; 
}

int TimeUniformRun(Interaction const* interaction, int n, int seed) {
    namespace time = std::chrono;

    // prepare grid
    dist::SetSeed(seed);
    Vec centre({0, 0, 0});
    double spread[Vec::mDim]{10, 10, 10};
    std::vector<Particle> particles = dist::MakeUniformMass(10, 5, n);
    particles = dist::SetUniformPos(centre, spread, particles);
    Grid grid = particles;

    // actual benchmarking
    auto start = time::high_resolution_clock::now();
    auto new_grid = interaction->Calculate(grid);
    auto end = time::high_resolution_clock::now();
    auto duration =
        time::duration_cast<std::chrono::microseconds>(end - start);

    return duration.count();
}

}

void TimeComplexity() {
    Log("Begin TimeComplexity.");
    Log("Make sure the governor is set to performance.");

    // need to find correct directory
    if (!CheckFile(dataDir, true)) { return; }

    const std::string dir_name = dataDir + "complexity/";
    const std::string comp_file_name = dir_name + "complexity.out";
    IO::MakeDir(dir_name);

    // do not overwrite
    if (CheckFile(comp_file_name, false)) { return; }

    // file IO
    std::ofstream stream;
    stream.open(comp_file_name);
    assert(stream.good());

    // parameters
    const int p = 3;
    const double theta = 0.5;
    const int maxPerCell = 5;
    const int maxPairwiseLimit = 3;
    const double G = -1;
    const auto invsq_force = std::make_unique<InvSqForce>(G);
    auto invsq_ker = std::make_unique<InvSqKernels>(p, G);

    // interaction specific settings
    const int int_types = 3;
    std::unique_ptr<Interaction> interactions[int_types] = {
        std::make_unique<Brute>(invsq_force.get()),
        std::make_unique<BarnesHut>(p,
                                    theta,
                                    invsq_ker.get(),
                                    invsq_force.get()),
        std::make_unique<FMM>(p,
                              theta,
                              maxPerCell,
                              maxPairwiseLimit,
                              invsq_ker.get(),
                              invsq_force.get())
    };
    const std::string int_names[int_types]{"brute", "bh", "fmm"};
    const int min_n = 10;
    const int int_max_n[int_types]{int(1e5), int(1e5), int(1e5)};

    // iteration settings
    const int runs_per_fac_10 = 4; 
    const double factor = std::pow(10, 1.0 / runs_per_fac_10);

    // repeats for each iteration
    const int repeats = 10;

    // log to file
    stream << int_types << std::endl;

    for(int i = 0; i < int_types; i++) {
        // prepare for iterations
        int max_n = int_max_n[i];
        double dn = min_n; // "double" n
        int runs = std::log(double(max_n) / min_n) / std::log(factor);

        Log("Time measurements for " + int_names[i]);

        stream << int_names[i] << " " << runs << std::endl;

        for (int j = 0; j <= runs; j++) {
            int n = int(dn);
            Log("Run " + std::to_string(j + 1) + " of "
                + std::to_string(runs + 1) + ", n = " + std::to_string(n));

            // timers for repeats
            long long time_sum = 0;
            
            stream << n << " " << repeats << std::endl;

            for (int k = 0; k < repeats; k++) {
                Log("Repeat " + std::to_string(k + 1) + " of "
                        + std::to_string(repeats)
                        );
                
                int repeat_time = TimeUniformRun(interactions[i].get(), n, k);
                time_sum += repeat_time;
                Log("Timed: " + std::to_string(repeat_time) + " us");

                stream << repeat_time << " ";
            }

            Log("Repeats finished. Average time: " + 
                    std::to_string(double(time_sum) / repeats));

            // measure next n
            dn *= factor;

            stream << std::endl;
        }
        
        // measure next interaction type
    }

}

}
}

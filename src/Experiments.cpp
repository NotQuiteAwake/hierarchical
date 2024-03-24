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

int TimeUniformRun(Interaction const* interaction,
                   int n,
                   int seed,
                   const std::string& fileName // dump destination
                   ) {
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

    // save results for further analysis
    std::ofstream stream;
    stream.open(fileName, std::ios::app);
    assert(stream.good());
    IO::SetHexfloat(stream);

    IO::DumpGrid(new_grid, stream);

    return duration.count();
}

}

void TimeComplexity() {
    std::clog << "Begin TimeComplexity." << std::endl;
    std::clog << "Make sure the governor is set to performance." << std::endl;

    // need to find correct directory
    if (!CheckFile(dataDir, true)) { return; }

    const std::string comp_dir_name = dataDir + "complexity/";
    const std::string comp_file_name = comp_dir_name + "complexity.out";
    const std::string dump_dir_name = comp_dir_name + "dump/";
    IO::MakeDir(comp_dir_name);

    // do not overwrite
    if (CheckFile(dump_dir_name, false)) { return; }
    if (CheckFile(comp_file_name, false)) { return; }

    // file IO
    IO::MakeDir(dump_dir_name);
    std::ofstream stream;
    stream.open(comp_file_name);
    assert(stream.good());
    IO::SetHexfloat(stream);

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
    const int runs_per_fac_10 = 8; 
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

        std::clog << "Time measurements for " + int_names[i] << std::endl;

        stream << int_names[i] << " " << runs << std::endl;

        for (int j = 0; j <= runs; j++) {
            int n = int(dn);
            std::clog << "Run " << j + 1 << " of " << runs + 1 << ", n = " << n
                << std::endl;

            // timers for repeats
            long long time_sum = 0;
            
            stream << n << " " << repeats << std::endl;

            for (int k = 0; k < repeats; k++) {
                std::clog << "Repeat " << k + 1 << " of " << repeats
                    << std::endl;
                
                std::stringstream dump_file_name_ss;
                dump_file_name_ss << dump_dir_name << int_names[i] << "_"
                    << n << ".dump";
                int repeat_time = TimeUniformRun(interactions[i].get(),
                                                 n,
                                                 k,
                                                 dump_file_name_ss.str());
                time_sum += repeat_time;
                std::clog << "Timed: " << repeat_time << " us" << std::endl;

                stream << repeat_time << " ";
            }

            std::clog << "Repeats finished. Average time: " << 
                (double(time_sum) / repeats) << std::endl;

            // measure next n
            dn *= factor;

            stream << std::endl;
        }
        
        // measure next interaction type
    }

}

void ExpansionOrderComplexity() {
    std::clog << "Begin ExpansionOrderComplexity." << std::endl;
    std::clog << "Make sure the governor is set to performance." << std::endl;

    // need to find correct directory
    if (!CheckFile(dataDir, true)) { return; }

    const std::string comp_dir_name = dataDir + "p-complexity/";
    const std::string comp_file_name = comp_dir_name + "complexity.out";
    const std::string dump_dir_name = comp_dir_name + "dump/";
    IO::MakeDir(comp_dir_name);

    // do not overwrite
    if (CheckFile(dump_dir_name, false)) { return; }
    if (CheckFile(comp_file_name, false)) { return; }

    // file IO
    IO::MakeDir(dump_dir_name);
    std::ofstream stream;
    stream.open(comp_file_name);
    assert(stream.good());
    IO::SetHexfloat(stream);

    // parameters
    const int n = 1000;
    const int p_min = 1;
    const int p_max = 10;
    const double theta = 0.5;
    const int maxPerCell = 5;
    const int maxPairwiseLimit = 3;
    const double G = -1;

    // repeats for each iteration
    const int repeats = 10;

    // interaction specific settings
    const int int_types = 3;
    const std::string int_names[int_types]{"brute", "bh", "fmm"};

    // log to file
    stream << (p_max - p_min + 1) << " " << n << std::endl;

    for (int p = p_min; p <= p_max; p++) {
        assert(p > 0);

        // Instantiate Interaction for the particular p
        // must do it for each p since by design we don't expect p to change
        // "half-way" in a calculation so it's meant to be immutable
        const auto invsq_force = std::make_unique<InvSqForce>(G);
        auto invsq_ker = std::make_unique<InvSqKernels>(p, G);

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

        // for p = p_min, calculate all three types: 3
        // for p != p_min, don't run brute-force again: 2
        stream << p << " " << (int_types - 1 + (p == p_min)) << std::endl;

        for(int i = 0; i < int_types; i++) {
            std::clog << int_names[i] << " p = " << p
                << ", n = " << n << std::endl;

            if (int_names[i] == "brute" && p != p_min) {
                // simply no point in repeating the same calculation
                continue;
            }

            // timers for repeats
            long long time_sum = 0;

            stream << int_names[i] << " " << repeats << std::endl;

            for (int k = 0; k < repeats; k++) {
                std::clog << "Repeat " << k + 1 << " of " << repeats
                    << std::endl;

                std::stringstream dump_file_name_ss;
                dump_file_name_ss << dump_dir_name << int_names[i] << "_"
                    << n << "_" << p << ".dump";
                int repeat_time = TimeUniformRun(interactions[i].get(),
                        n,
                        k,
                        dump_file_name_ss.str());
                time_sum += repeat_time;
                std::clog << "Timed: " << repeat_time << " us" << std::endl;

                stream << repeat_time << " ";
            }

            std::clog << "Repeats finished. Average time: " << 
                (double(time_sum) / repeats) << std::endl;

            stream << std::endl;
        }

        // measure next interaction type
    }
}

}
}

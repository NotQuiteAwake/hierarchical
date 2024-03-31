/**
 * @file
 * @brief Implementation of Benchmarker for testing and benchmarking.
 */

#include <chrono>
#include <fstream>
#include <iomanip>

#include "Benchmarker.hpp"
#include "Distribution.hpp"
#include "IO.hpp"
#include "Brute.hpp"
#include "Barnes-Hut.hpp"
#include "FMM.hpp"

namespace sim {

namespace exp {


/**
 * @class FileParams
 * @brief Parameters related to file storage in an experiment.
 */
FileParams::FileParams(std::string runDir): runDir(runDir) {};

/**
 * @brief Get relative path to the directory where experiment will run.
 */
std::string FileParams::GetRunDirName() const {
    return dataDir + '/' + runDir;
}

/**
 * @brief Get relative path to file where complexity timing should go
 */
std::string FileParams::GetCompName() const {
    return GetRunDirName() + '/' + comp;
}

/**
 * @brief Get relative path to directory where grids should be dumped
 *
 */
std::string FileParams::GetDumpName() const {
    return GetRunDirName() + '/' + dump;
}

/**
 * @class GenParams
 * @brief General parameters to be used in an experiment.
 */

/**
 * @class Benchmarker
 * @brief Class to help benchmark and test different algorithms.
 */

/**
 * @brief Time a single run with a uniform distritbuion
 *
 * @param[in] interaction Force calculation method to use
 * @param[in] n Number of particles in distribution
 * @param[in] seed Seed for random distribution generation
 * @param[in, out] stream File stream to dump processed grid to
 * @return Timing result in microseconds
 */
int Benchmarker::TimeUniformRun(Interaction const* interaction,
                   int n,
                   int seed,
                   std::ofstream& stream
                   ) const {
    namespace time = std::chrono;

    // prepare grid
    dist::SetSeed(seed);
    Vec centre({0, 0, 0});
    double spread[Vec::mDim]{genParams.scale, genParams.scale, genParams.scale};
    std::vector<Particle> particles = dist::MakeUniformMass(10, 5, n);
    particles = dist::SetUniformPos(centre, spread, particles);
    Grid grid = particles;

    // actual benchmarking
    auto start = time::high_resolution_clock::now();
    auto new_grid = interaction->Calculate(grid);
    auto end = time::high_resolution_clock::now();
    auto duration =
        time::duration_cast<std::chrono::microseconds>(end - start);

    IO::DumpGrid(new_grid, stream);

    return duration.count();
}

/**
 * @brief Initialised Benchmarker
 *
 * @param[in] genParams General parameters for simulations
 * @param[in] fileParams File-related parameters
 * @param[in, out] stream Stream to dump finished grids to
 */
Benchmarker::Benchmarker(GenParams genParams,
                         FileParams fileParams,
                         std::ostream& stream):
        stream(stream),
        genParams(genParams),
        fileParams(fileParams) {
    
    force.reset(new InvSqForce(genParams.G));
    kernels.reset(new InvSqKernels(genParams.p, genParams.G));
    integrator.reset(new LeapFrog(genParams.step));

    interactions[genParams.brute] =
        std::make_unique<Brute>(force.get());
    interactions[genParams.bh] =
        std::make_unique<BarnesHut>(genParams.p,
                genParams.theta,
                kernels.get(),
                force.get());
    interactions[genParams.fmm] =
        std::make_unique<FMM>(genParams.p,
                genParams.theta,
                genParams.maxPerCell,
                genParams.maxPairwiseLimit,
                kernels.get(),
                force.get());
}


/**
 * @brief Run the repeats by configuration in genParams and fileParams.
 *
 * The Benchmarker itself is agnostic to which parameter is actually being
 * varied, but this needs to logged to file so we provide the value (doesn't
 * matter which param it actually is) to the method for it to output.
 *
 * @param[in] param Value of parameter being varied in this run
 * @param[in] skipBrute Whether to skip run with brute-force
 */
void Benchmarker::Run(double param, bool skipBrute) const {
    stream << param << " " << (genParams.intTypes - skipBrute) << std::endl;
    for(int i = 0; i < genParams.intTypes; i++) {
        const std::string& int_name = genParams.intNames[i];

        if (skipBrute && i == genParams.brute) { continue; }

        // timers for repeats
        long long time_sum = 0;

        std::clog << int_name << " param " << param << std::endl;
        stream << int_name << " " << genParams.repeats << std::endl;

        std::stringstream dump_file_name_ss;
        dump_file_name_ss << fileParams.GetDumpName() << int_name << "_"
            << genParams.n << "_" << param << ".dump";
        std::ofstream dump_stream(dump_file_name_ss.str());
        assert(dump_stream.good());
        dump_stream << genParams.repeats << std::endl;
        IO::SetHexfloatOut(dump_stream);

        for (int j = 0; j < genParams.repeats; j++) {
            std::clog << "Repeat " << j + 1 << " of " << genParams.repeats
                << std::endl;

            int repeat_time = TimeUniformRun(interactions[i].get(),
                    genParams.n,
                    j,
                    dump_stream);
            time_sum += repeat_time;

            std::clog << "Timed: " << repeat_time << " us" << std::endl;
            stream << repeat_time << " ";
        }
        std::clog << "Repeats finished. Average time: " << 
            (double(time_sum) / genParams.repeats) << std::endl;
        stream << std::endl;
    }
}

/**
 * @brief Perform time evolution by parameters in genParams and fileParams
 *
 * @param[in] grid Initial condition for simulation.
 */
void Benchmarker::Evo(const Grid& grid) const {
    if (!IO::CheckFile(FileParams::dataDir, true)) { return; }
    IO::MakeDir(fileParams.GetRunDirName());
    if (IO::CheckFile(fileParams.GetDumpName(), false)) { return; }
    IO::MakeDir(fileParams.GetDumpName());

    const int steps_cnt = genParams.evo_time / genParams.step;
   
    namespace time = std::chrono;
    std::clog << std::fixed << std::setprecision(2);

    for (int i = 0; i < genParams.intTypes; i++) {
        const std::string& int_name = genParams.intNames[i];

        Grid sim_grid = grid;
        Interaction* const interaction = interactions[i].get();
        
        const std::string file_name = fileParams.GetDumpName() 
            + genParams.intNames[i] + ".dump";
        if (IO::CheckFile(file_name, false)) { return; }
        std::ofstream stream(file_name);
        assert(stream.good());
        IO::SetHexfloatOut(stream);

        std::clog << int_name << std::endl;
        stream << int_name << std::endl;
        stream << steps_cnt << " " << genParams.step <<
            " " << genParams.scale << std::endl;

        sim_grid = interaction->Calculate(grid);
        for (int j = 0; j < steps_cnt; j++) {
            const double t = j * genParams.step;
            if (j > 0 && j % 10 == 0) std::clog << std::endl;
            std::clog << t << " ";
            IO::DumpGrid(sim_grid, stream);

            auto start = time::high_resolution_clock::now();
            sim_grid = integrator->Evolve(sim_grid);
            sim_grid = interaction->Calculate(sim_grid);
            auto end = time::high_resolution_clock::now();
            auto duration =
                time::duration_cast<std::chrono::microseconds>(end - start);

            stream << duration.count() << " " << sim_grid.GetPE() << " "
                   << sim_grid.GetKE() << std::endl;
            IO::DumpVec(sim_grid.GetL(), stream);
            IO::DumpVec(sim_grid.GetP(), stream);
        }

        std::clog << std::endl;    
    }
}



}

}

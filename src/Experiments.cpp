#include <string>
#include <cmath>
#include <memory>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "Experiments.hpp"

#include "Interaction.hpp"
#include "Integrator.hpp"
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

FileParams::FileParams(std::string runDir): runDir(runDir) {};

std::string FileParams::GetRunDirName() const {
    return dataDir + '/' + runDir;
}

std::string FileParams::GetCompName() const {
    return GetRunDirName() + '/' + comp;
}

std::string FileParams::GetDumpName() const {
    return GetRunDirName() + '/' + dump;
}

bool CompStreamSetup(std::ofstream& stream,
                     FileParams fileParams
        ) {
    const std::string comp_dir_name = fileParams.GetRunDirName();
    const std::string comp_file_name = fileParams.GetCompName();
    const std::string dump_dir_name = fileParams.GetDumpName();

     // need to find correct directory
    if (!CheckFile(FileParams::dataDir, true)) { return false; }

    IO::MakeDir(comp_dir_name);

    // do not overwrite
    if (CheckFile(dump_dir_name, false)) { return false; }
    if (CheckFile(comp_file_name, false)) { return false; }

    // file IO
    IO::MakeDir(dump_dir_name);
    stream.open(comp_file_name);
    assert(stream.good());
    IO::SetHexfloat(stream);

    return true;
}

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

Benchmarker::Benchmarker(GenParams genParams,
                         FileParams fileParams,
                         std::ofstream& stream):
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
        IO::SetHexfloat(dump_stream);

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

void Benchmarker::Evo(const Grid& grid) const {
    if (!CheckFile(FileParams::dataDir, true)) { return; }
    IO::MakeDir(fileParams.GetRunDirName());
    if (CheckFile(fileParams.GetDumpName(), false)) { return; }
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
        if (CheckFile(file_name, false)) { return; }
        std::ofstream stream(file_name);
        assert(stream.good());
        IO::SetHexfloat(stream);

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

} // namespace

void NComplexity() {
    std::clog << "Begin NComplexity." << std::endl;
    std::clog << "Make sure the governor is set to performance." << std::endl;

    FileParams fileParams("n-complexity/");
    std::ofstream stream;
    if (!CompStreamSetup(stream, fileParams)) { return; }

    const int min_n = 10;
    const int max_n = int(1e5 + 0.5);
    const int runs_per_fac_10 = 8; 
    const double factor = std::pow(10, 1.0 / runs_per_fac_10);
    const int num_n = std::log(double(max_n) / min_n) / std::log(factor) + 1;

    double n_double = min_n;

    stream << num_n << std::endl;

    for (int i = 0; i < num_n; i++) {
        int n = int(n_double + 0.5);
        
        GenParams params;
        params.n = n;
        Benchmarker bm(params, fileParams, stream);
        
        // don't skipBrute, we need that in measurements as well
        bm.Run(n, false);

        n_double *= factor;
    }
}

void PComplexity() {
    std::clog << "Begin PComplexity." << std::endl;
    std::clog << "Make sure the governor is set to performance." << std::endl;
    
    FileParams fileParams("p-complexity/");
    std::ofstream stream;
    if (!CompStreamSetup(stream, fileParams)) { return; }

    const int p_min = 1;
    const int p_max = 10;
    const int num_p = p_max - p_min + 1;
     
    // TODO: REMOVED N from output
    stream << num_p << std::endl;

    for (int i = 0; i < num_p; i++) {
        const int p = p_min + i;
        assert(p > 0);
        
        GenParams params;
        params.p = p;
        Benchmarker bm(params, fileParams, stream);

        bm.Run(p, i != 0);
    }
}

void ThetaComplexity() {
    std::clog << "Begin ThetaComplexity." << std::endl;
    std::clog << "Make sure the governor is set to performance." << std::endl;

    FileParams fileParams("theta-complexity/");
    std::ofstream stream;
    if (!CompStreamSetup(stream, fileParams)) { return; }

    const double theta_min = 0.001;
    const double theta_max = 1;
    const int runs_per_fac_10 = 8; 
    const double factor = std::pow(10, 1.0 / runs_per_fac_10);
    const int num_theta = std::log(theta_max / theta_min)
        / std::log(factor) + 1;

    double theta = theta_min;

    // log to file
    stream << num_theta << std::endl;

    for (int i = 0; i < num_theta; i++) {
        assert(theta > 0);
    
        GenParams params;
        params.theta = theta;
        Benchmarker bm(params, fileParams, stream);

        bm.Run(theta, i != 0);

        theta *= factor;
    }
}

void ColdStartSim() {
    GenParams genParams;
    genParams.scale = 20;
    FileParams fileParams("cold/");
    Benchmarker bm(genParams, fileParams);

    const double seed = 1;
    dist::SetSeed(seed);

    Grid grid(genParams.n);
    const double mean_mass = 10;
    const double sigma_mass = 1;
    const Vec centre({0, 0, 0});
    const double R = genParams.scale / 2;
    auto par_list = dist::MakeNormalMass(mean_mass, sigma_mass, genParams.n);
    par_list = dist::SetSphericalPos(centre, 0, R, par_list);
   
    grid.AddParticles(par_list);
    bm.Evo(grid);
}

void ThinDiskSim() {
    GenParams genParams;
    genParams.scale = 20;
    genParams.n = 500;
    FileParams fileParams("disk/");
    Benchmarker bm(genParams, fileParams);

    const double seed = 10;
    dist::SetSeed(seed);

    Grid grid(genParams.n);

    const double mean_mass = 10;
    const double sigma_mass = 1;

    const Vec centre({0, 0, 0});
    const Vec axis({0, 0, 1});
    const double R = genParams.scale;
    const double z_spread = 0.5;

    auto par_list = dist::MakeNormalMass(mean_mass, sigma_mass, genParams.n);
    par_list = dist::SetDiskPos(centre, axis, 0, R, z_spread, par_list);

    {
        // find PE. REQUIRES SimulationHelper to have G = -1.
        Grid grid(par_list);
        auto force = std::make_unique<InvSqForce>();
        auto brute = std::make_unique<Brute>(force.get());
        grid = brute->Calculate(grid);
        const double PE = grid.GetPE();

        // expectation from virial thm
        const double KE = -1.0 / 2 * PE;
        const double I = 1.0 / 2 * (mean_mass * genParams.n) * (R * R); // MOI
        const Vec omega = axis * std::sqrt(2 * KE / I);
        par_list = dist::SetUniformRotVel(centre, omega, par_list);
    }

    grid.AddParticles(par_list);
    bm.Evo(grid);
}

void GalaxySim() {
    GenParams genParams;
    genParams.scale = 50;
    FileParams fileParams("galaxy/");
    Benchmarker bm(genParams, fileParams);

    const double seed = 10;
    dist::SetSeed(seed);

    Grid grid(genParams.n);

    const Vec centre({0, 0, 0});
    const Vec axis({0, 0, 1});

    const double mean_mass = 40;
    const double sigma_mass = 7;
    const double smbh_mass = 5 * mean_mass * genParams.n;
    Particle smbh(smbh_mass, smbh_mass);
    smbh.pos = centre;

    const double scale = genParams.scale;
    const double r0 = scale / 10;
    const double r1 = scale;
    const double z_spread = 5;

    auto par_list = dist::MakeNormalMass(mean_mass, sigma_mass, genParams.n);
    par_list = dist::SetDiskPos(centre, axis, r0, r1, z_spread, par_list);

    const double G = 1;
    const double alpha = G * smbh_mass;
    par_list = dist::SetCircVel(centre, axis, alpha, par_list);

    par_list.push_back(smbh);

    grid.AddParticles(par_list);
    bm.Evo(grid);
}

void TwoGalaxiesSim() {
    GenParams genParams;
    genParams.scale = 75;
    FileParams fileParams("two-galaxies/");
    Benchmarker bm(genParams, fileParams);

    const double seed = 2333;
    dist::SetSeed(seed);

    const int n = genParams.n;
    Grid grid(2 * n);

    const double mean_mass = 40;
    const double sigma_mass = 7;
    const double smbh_mass = 5 * mean_mass * n;

    const double scale = genParams.scale;
    const double r1 = scale / 2;
    const double r0 = r1 / 7;
    const double z_spread = r1 / 10;
    const Vec centres[2] = {Vec({-scale / 2.5, 0, 0}), Vec({scale / 2.5, 0, 0})};
    const double dist_centres = scale / 2;
    const Vec axes[2] = {Vec({1, -0.5, 1}), Vec({-0.5, -0.5, 1})};

    for (int i = 0; i < 2; i++) {
        const Vec& centre = centres[i];
        const Vec& axis = axes[i];

        Particle smbh(smbh_mass, smbh_mass);
        smbh.pos = centre;

        auto par_list = dist::MakeNormalMass(mean_mass, sigma_mass, n);
        par_list = dist::SetDiskPos(centre, axis, r0, r1, z_spread, par_list);

        const double G = 1;
        const double alpha = G * smbh_mass;
        par_list = dist::SetCircVel(centre, axis, alpha, par_list);

        par_list.push_back(smbh);
        // circular motion of the smbh
        const Vec boost = Vec({0, 1, 0}) * (i ? 1 : -1)
                          * std::sqrt(G * smbh_mass / (dist_centres * 2));
        par_list = dist::AddUniformVel(boost, par_list);
        
        grid.AddParticles(par_list);
    }

    bm.Evo(grid);
}

}
}

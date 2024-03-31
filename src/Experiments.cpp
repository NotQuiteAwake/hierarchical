/**
 * @file
 * @brief Collection of functions for running different experiments on the code.
 */

#include <cmath>
#include <memory>
#include <fstream>
#include <iostream>

#include "Experiments.hpp"

#include "Brute.hpp"
#include "Distribution.hpp"
#include "Benchmarker.hpp"
#include "IO.hpp"

namespace sim {
namespace exp {

namespace {
/**
 * @brief Setup file stream for complexity analyses
 *
 * Performs a number of checks to prevent overwriting existing experimental
 * data, and create directories if they don't already exist.
 *
 * @param[in, out] stream File stream to be initialised
 * @param[in] fileParams Parameters for file storage
 * @return Status of whether the setup attempt was successful.
 */
bool CompStreamSetup(std::ofstream& stream,
                     FileParams fileParams
        ) {
    const std::string comp_dir_name = fileParams.GetRunDirName();
    const std::string comp_file_name = fileParams.GetCompName();
    const std::string dump_dir_name = fileParams.GetDumpName();

     // need to find correct directory
    if (!IO::CheckFile(FileParams::dataDir, true)) { return false; }

    IO::MakeDir(comp_dir_name);

    // do not overwrite
    if (IO::CheckFile(dump_dir_name, false)) { return false; }
    if (IO::CheckFile(comp_file_name, false)) { return false; }

    // file IO
    IO::MakeDir(dump_dir_name);
    stream.open(comp_file_name);
    assert(stream.good());
    IO::SetHexfloatOut(stream);

    return true;
}

}

/**
 * @brief Measure complexity of three "Interaction"s against number of particles
 *
 */
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

/**
 * @brief Measure complexity of algorithms against order of multipole expansion
 *
 */
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

/**
 * @brief Measure complexity of algorithms against opening angle
 *
 */
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

/**
 * @brief Generate and simulate with a cold start initial condition.
 *
 */
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

/** 
 * @brief Generate and simulate with a thin disk initial condition.
 *
 * In this condition, particles of similar masses are uniformly added into a
 * disk, and set to rotate at the same angular velocity, approximated from a use
 * of the Virial theorem.
 *
 * This is not presented or used in the end as this condition is unstable; The
 * system simply disintegrates, so there's little to see. An earlier attempt to
 * simply model the "thin disk" as a "cylinder" and use angular velocity from
 * that also failed tragically (of course that would...)
 *
 * It is this that prompted me to use the SMBHs at the centre.
 *
 */
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

/**
 * @brief Generate and simulate with a single galaxy initial condition
 *
 * Thin disk with SMBH at its centre, so velocities are assumed to be just from
 * circular orbits around the SMBH.
 */
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

/**
 * @brief Generate and simulate with initial conditions with two galaxies
 *
 * Each galaxy has a SMBH at the centre, and they are orbiting each other.
 */
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

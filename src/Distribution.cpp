/** 
 * @file 
 * @brief Generator of Particle distributions for use in experiments.
 *
 * */

#include <cassert>
#include <random>

#include "Distribution.hpp" 

namespace sim {

namespace dist {

namespace {

static constexpr double PI() { return std::atan(1) * 4; }

double massCutoff = 0.05;
int seed = 0;
std::mt19937 engine;

/**
 * @brief Generate normal distribution
 *
 * @param[in] mean mean value for the distribution
 * @param[in] stdev standard deviation from mean
 * @param[in] n number of numbers to generate
 * @return list of random numbers
 */
std::vector<double> GenNormal(double mean, double stdev, int n) {
    assert(stdev >= 0);
    std::vector<double> results;
    results.reserve(n);
    std::normal_distribution<double> gauss_dist(mean, stdev);
    for (int i = 0; i < n; i++) results.push_back(gauss_dist(engine));
    return results;
}

/**
 * @brief generate uniform distribution from left and right boundaries
 *
 * @param[in] left left limit of distribution
 * @param[in] right right limit of distribution
 * @param[in] n number of numbers to generate
 * @return list of random numbers
 */
std::vector<double> GenUniformLR(double left, double right, int n) {
    assert(left <= right);
    std::vector<double> results;
    results.reserve(n);
    std::uniform_real_distribution<double> uniform_dist(left, right);
    for (int i = 0; i < n; i++) results.push_back(uniform_dist(engine));
    return results;
}

/**
 * @brief Generate uniform distribution from mean and spread
 *
 * @param[in] mean centre of uniform distribution
 * @param[in] spread mean +- spread gives distribution boudaries
 * @param[in] n number of numbers to generate
 * @return list of random numbers
 */
std::vector<double> GenUniformMS(double mean, double spread, int n) {
    assert(spread >= 0);
    double left = mean - spread;
    double right = mean + spread;
    return GenUniformLR(left, right, n);
}

/**
 * @brief Set masses in list to massCutoff if it went below. (destructive)
 *
 * This is so that all masses will be positive, even if some continuous
 * distributions yield masses below zero (e.g. normal distribution).
 * 
 * @param[in, out] masses List of masses to be checked and modified.
 */
std::vector<double>& TrimMaases(std::vector<double>& masses) {
    for (int i = 0; i < masses.size(); i++) {
        if (masses[i] < massCutoff) { masses[i] = massCutoff; }
    } 
    return masses;
}

}

/**
 * @brief Set random number seed
 */
void SetSeed(int newSeed) { 
    seed = newSeed;
    engine.seed(seed);
}

int GetSeed() { return seed; }

void SetMassCutoff(double newCutoff) { massCutoff = newCutoff; }

double GetMassCutoff() { return massCutoff; }

/**
 * @brief Make a list of particles whose mass obeys uniform statistics
 *
 * @param[in] mean Centre of uniform distribution
 * @param[in] spread mean +- spread gives distribution boudaries
 * @param[in] n Number of particles to generate masses for
 * @return A list of particles, in which only mass (= charge) are set
 */
std::vector<Particle> MakeUniformMass(double mean, double spread, int n) {
    std::vector<double> masses = GenUniformMS(mean, spread, n);
    masses = TrimMaases(masses);
    std::vector<Particle> results;
    for (double mass : masses) { results.push_back(Particle(mass, mass)); }
    return results;
}

/**
 * @brief Make a list of particles whose mass obeys normal distribution
 *
 * @param[in] mean Centre of normal distribution
 * @param[in] sigma Standard deviation
 * @param[in] n Number of particles to generate masses for
 * @return A list of particles, in which only mass (= charge) are set
 */
std::vector<Particle> MakeNormalMass(double mean, double sigma, int n) {
    std::vector<double> masses = GenNormal(mean, sigma, n);
    masses = TrimMaases(masses);
    std::vector<Particle> results;
    for (double mass : masses) { results.push_back(Particle(mass, mass)); }
    return results;
}

/**
 * @brief Set positions of particles in a list by a uniform distribution
 *
 * @param[in] posCentre Displacement of the centre of the distribution
 * @param[in] spread Spread by axis
 * @param[in, out] list List of particles whose positions will be populated
 * @return Modified list of particles
 */
std::vector<Particle>& SetUniformPos(
        const Vec& posCentre,
        double (&spread)[Vec::mDim], // spread in all three axes
        std::vector<Particle>& list
        ) {
    int n = list.size();
    for (int i = 0; i < Vec::mDim; i++) {
        std::vector<double> x_i = GenUniformMS(posCentre[i], spread[i], n);
        for (int j = 0; j < n; j++) {
            list[j].pos[i] = x_i[j];
        } 
    } 
    return list;
}

/**
 * @brief Set positions of particles in a list by a normal distribution
 * 
 * @param[in] posCentre Displacement of the centre of the distribution
 * @param[in] sigma Standard deviation by axis
 * @param[in, out] list List of particles whose positions will be populated
 * @return Modified list of particles
 */
std::vector<Particle>& SetNormalPos(
        const Vec& posCentre,
        double (&sigma)[Vec::mDim],
        std::vector<Particle>& list 
        ) {
    int n = list.size();
    for (int i = 0; i < Vec::mDim; i++) {
        std::vector<double> x_i = GenNormal(posCentre[i], sigma[i], n);
        for (int j = 0; j < n; j++) {
            list[j].pos[i] = x_i[j];
        } 
    } 
    return list;
}

/**
 * @brief Set position of particles by uniform distribution in a sphere
 *
 * @param[in] posCenter Centre of the sphere
 * @param[in] r0 Minimum allowed radius for position
 * @param[in] r1 Maximum allowed radius for position
 * @param[in, out] list List of particles whose positions will be populated
 * @return Modified list of particles
 */
std::vector<Particle>& SetSphericalPos(
        const Vec& posCentre,
        double r0,
        double r1,
        std::vector<Particle>& list
        ) {
    int n = list.size();
    std::vector<double> r = GenUniformLR(r0, r1, n);
    std::vector<double> theta = GenUniformLR(0, PI(), n);
    std::vector<double> phi = GenUniformLR(0, 2 * PI(), n);
    
    for (int i = 0; i < n; i++) {
        list[i].pos = Vec::FromSpherical(r[i], theta[i], phi[i]);
        // velocity is 0 in cold start, as is also initialised to be.
    }
    return list; 
}

/**
 * @brief Set particles velocities by a uniform angular velocity about axis
 *
 * @param[in] posCenter Centre of rotation
 * @param[in] omega Angular momentum vector
 * @param[in, out] list List of particles whose positions will be populated
 * @return Modified list of particles
 */
std::vector<Particle>& SetUniformRotVel(
        const Vec& posCentre,
        const Vec& omega,
        std::vector<Particle>& list
        ) {
    for (Particle& par : list) {
        const Vec r = par.pos - posCentre;
        par.vel = CrossProduct(omega, r);
    }    
    return list;
}

/**
 * @brief Set particle positions by random distribution on a disk.
 *
 * @param[in] posCenter Centre of the disk
 * @param[in] axis Axis perpendicular to plane of disk
 * @param[in] r0 Minimum allowed radius
 * @param[in] r1 Maximum allowed radius
 * @param[in] zSigma Standard deviation for height fluctuation along axis
 * @param[in, out] list List of particles whose positions will be populated
 * @return Modified list of particles
 */
std::vector<Particle>& SetDiskPos(
        const Vec& posCentre,
        const Vec& axis,
        double r0,
        double r1,
        double zSigma,
        std::vector<Particle>& list
        ) {
    list = SetSphericalPos(Vec(), r0, r1, list);

    int n = list.size();
    std::vector<double> z = GenNormal(0, zSigma, n);

    for (int i = 0; i < n; i++) {
        Particle& par = list[i];
        Vec& pos = par.pos;
        double par_norm = pos.GetNorm();

        // map into plane of disk
        pos = CrossProduct(pos, axis) / axis.GetNorm();
        // but maintain the uniform distribution in r
        pos = pos / pos.GetNorm() * par_norm;
        // finally shift from origin and offset along axis
        pos += posCentre + axis * z[i];
    }
    return list;
}

/**
 * @brief Set particles to velocities required for circular orbits
 * 
 * This assumes an inverse square law, in which angular velocity at a distance r
 * from the central massive body is \f$\omega = \sqrt{\alpha / r^3}\f$.
 *
 * @param[in] centre Centre of circular orbits
 * @param[in] axis Rotational axis
 * @param[in] alpha A strengh parameter. For gravity alpha = GM
 * @param[in, out] list List of particles whose positions will be populated
 * @return Modified list of particles
 */
std::vector<Particle>& SetCircVel(
        const Vec& centre,
        const Vec& axis,
        const double alpha,
        std::vector<Particle>& list
        ) {
    for (Particle& par : list) {
        const Vec r = par.pos - centre;
        const double dist = r.GetNorm();
        const double omega_val = sqrt(alpha / (dist * dist * dist));
        const Vec omega = axis / axis.GetNorm() * omega_val;
        par.vel = CrossProduct(omega, r);
    }
    return list;
}

/**
 * @brief Boost all particles by velocity
 *
 * @param boost Velocity to boost particles by
 * @param[in, out] list List of particles whose positions will be populated
 * @return Modified list of particles
 */
std::vector<Particle>& AddUniformVel(
        const Vec& boost,
        std::vector<Particle>& list
        ) {
    for (Particle& par : list) {
        par.vel += boost;
    }
    return list;
}

}
}

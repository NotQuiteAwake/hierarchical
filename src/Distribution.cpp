#include "Distribution.hpp" 

namespace sim {

namespace dist {

namespace {

static constexpr double PI() { return std::atan(1) * 4; }

std::vector<double> GenNormal(double mean, double stdev, int n) {
    std::vector<double> results;
    results.reserve(n);
    std::normal_distribution<double> gauss_dist(mean, stdev);
    for (int i = 0; i < n; i++) results.push_back(gauss_dist(engine));
    return results;
}

std::vector<double> GenUniform(double left, double right, int n) {
    std::vector<double> results;
    results.reserve(n);
    std::uniform_real_distribution<double> uniform_dist(left, right);
    for (int i = 0; i < n; i++) results.push_back(uniform_dist(engine));
    return results;
}

std::vector<double>& TrimMaases(std::vector<double>& masses) {
    for (int i = 0; i < masses.size(); i++) {
        if (masses[i] < massCutoff) { masses[i] = massCutoff; }
    } 
    return masses;
}

}

void SetSeed(int newSeed) { 
    seed = newSeed;
    engine.seed(seed);
}

void SetMassCutoff(double newCutoff) { massCutoff = newCutoff; }

std::vector<Particle> MakeUniformMass(double mean, double spread, int n) {
    std::vector<double> masses = GenUniform(mean, spread, n);
    masses = TrimMaases(masses);
    std::vector<Particle> results;
    for (double mass : masses) { results.push_back(Particle(mass, mass)); }
    return results;
}

std::vector<Particle> MakeNormalMass(double mean, double sigma, int n) {
    std::vector<double> masses = GenNormal(mean, sigma, n);
    masses = TrimMaases(masses);
    std::vector<Particle> results;
    for (double mass : masses) { results.push_back(Particle(mass, mass)); }
    return results;
}

std::vector<Particle>& SetUniformPos(
        const Vec& posCentre,
        double (&spread)[Vec::mDim], // spread in all three axes
        std::vector<Particle>& list
        ) {
    int n = list.size();
    for (int i = 0; i < Vec::mDim; i++) {
        std::vector<double> x_i = GenUniform(posCentre[i], spread[i], n);
        for (int j = 0; j < n; j++) {
            list[j].pos[i] = x_i[j];
        } 
    } 
    return list;
}

std::vector<Particle>& SetNormalPos(
        const Vec& posCentre,
        double (&sigma)[Vec::mDim],
        std::vector<Particle>& list 
        ) {
    int n = list.size();
    for (int i = 0; i < Vec::mDim; i++) {
        std::vector<double> x_i = GenUniform(posCentre[i], sigma[i], n);
        for (int j = 0; j < n; j++) {
            list[j].pos[i] = x_i[j];
        } 
    } 
    return list;
}

std::vector<Particle>& SetSphericalPos(
        const Vec& posCentre,
        double r0,
        double r1,
        std::vector<Particle>& list
        ) {
    int n = list.size();
    std::vector<double> r = GenUniform(r0, r1, n);
    std::vector<double> theta = GenUniform(0, PI(), n);
    std::vector<double> phi = GenUniform(0, 2 * PI(), n);
    
    for (int i = 0; i < n; i++) {
        list[i].pos = Vec::FromSpherical(r[i], theta[i], phi[i]);
        // velocity is 0 in cold start, as is also initialised to be.
    }
    return list; 
}

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

std::vector<Particle>& SetDiskPos(
        const Vec& posCentre,
        const Vec& axis,
        double r0,
        double r1,
        double z_spread,
        std::vector<Particle>& list
        ) {
    list = SetSphericalPos(Vec(), r0, r1, list);

    int n = list.size();
    std::vector<double> z = GenNormal(0, z_spread, n);

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

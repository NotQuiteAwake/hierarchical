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

std::vector<Particle>& SetColdStart(
        const Vec& posCentre,
        double R,
        std::vector<Particle>& list
        ) {
    int n = list.size();
    std::vector<double> r = GenUniform(0, R, n);
    std::vector<double> theta = GenUniform(0, PI(), n);
    std::vector<double> phi = GenUniform(0, 2 * PI(), n);
    
    for (int i = 0; i < n; i++) {
        list[i].pos = Vec::FromSpherical(r[i], theta[i], phi[i]);
        // velocity is 0 in cold start, as is also initialised to be.
    }
    return list; 
}
    
}
}

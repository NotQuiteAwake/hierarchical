#ifndef DISTRIBUTIONHEADERDEF
#define DISTRIBUTIONHEADERDEF

#include <vector>
#include <random>
#include "Particle.hpp"
#include "Octant.hpp"

namespace sim {

namespace dist {

namespace {
double massCutoff = 0.05;
int seed = 0;
std::mt19937 engine;

std::vector<double> GenNormal(double mean, double stdev, int n);
std::vector<double> GenUniform(double left, double right, int n);
std::vector<double>& TrimMasses(std::vector<double>& masses); 

}

void SetSeed(int newSeed);
void SetMassCutoff(double newCutoff);

// TODO: works for gravity only right now.
std::vector<Particle> MakeUniformMass(double mean, double spread, int n);
std::vector<Particle> MakeNormalMass(double mean, double sigma, int n);

std::vector<Particle>& SetUniformPos(
        const Vec& posCentre,
        double (&spread)[Vec::mDim],
        std::vector<Particle>& list
        );
std::vector<Particle>& SetNormalPos(
        const Vec& posCentre,
        double (&sigma)[Vec::mDim],
        std::vector<Particle>& list
        );
std::vector<Particle>& SetSphericalPos(
        const Vec& posCentre,
        double R,   // maximum radius
        std::vector<Particle>& list
        );
std::vector<Particle>& SetUniformRotVel(
        const Vec& posCentre,
        const Vec& omega,
        std::vector<Particle>& list
        );
std::vector<Particle>& SetPlummerSphere(
        const Vec& posCentre,
        double a, // Plummer parameter a
        std::vector<Particle>& list
        );
std::vector<Particle>& SetDiskPos(
        const Vec& posCentre,
        const Vec& axis, // axis of disk
        double r0, // r_0
        double z_spread, // spread along axis
        std::vector<Particle>& list
        );
std::vector<Particle>& SetUniformVel(
        const Vec& velCentre,
        double spread,
        std::vector<Particle>& list
        );
std::vector<Particle>& SetNormalVel(
        const Vec& velCentre,
        double spread,
        std::vector<Particle>& list
        );
}

}

#endif


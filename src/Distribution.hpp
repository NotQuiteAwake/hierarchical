#ifndef DISTRIBUTIONHEADERDEF
#define DISTRIBUTIONHEADERDEF

#include <vector>
#include "Particle.hpp"

namespace sim {

namespace dist {

void SetSeed(int newSeed);
int GetSeed();
void SetMassCutoff(double newCutoff);
double GetMassCutoff();

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
        double r0,   // minimum radius
        double r1,   // maximum radius
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
        double r0, // min radius
        double r1, // max radius
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
std::vector<Particle>& SetCircVel(
        const Vec& centre,
        const Vec& axis,
        const double alpha,
        std::vector<Particle>& list
        );
std::vector<Particle>& AddUniformVel(
        const Vec& boost,
        std::vector<Particle>& list
        );
}

}
#endif


#ifndef IOHEADERDEF
#define IOHEADERDEF

#include <iostream>
#include <string>
#include <memory>
#include "Octree.hpp"
#include "Grid.hpp"

namespace sim {

namespace IO {

template<typename T> static std::vector<T> LoadVector(std::istream& stream);

template<typename T> static void DumpVector(const std::vector<T> vector,
        std::ostream& stream);

Octant LoadOctant(std::istream& stream);
Vec LoadVec(std::istream& stream);
Particle LoadParticle(std::istream& stream);


void DumpOctant(const Octant& octant, std::ostream& stream);
void DumpVec(const Vec& vec, std::ostream& stream);
void DumpParticle(const Particle& par, std::ostream& stream);
void DumpParticle(const Particle& par);

void DumpGrid(const Grid& grid, std::ostream& stream);
void DumpGrid(const Grid& grid); // to stdout
void DumpGrid(const Grid& grid, const std::string& file); // to file

Grid LoadGrid(std::istream& stream);
Grid LoadGrid(); // from stdin
Grid LoadGrid(const std::string& file); // from file

static void DumpOctreeNode(Octree const* node, std::ostream& stream);
static void DumpOctreeHelper(Octree const* node, std::ostream& stream);
void DumpOctree(Octree const* octree, std::ostream& stream);
void DumpOctree(Octree const* octree); // to stdout
void DumpOctree(Octree const* octree, const std::string& file); // to file

static void LoadOctreeNode(Octree* node, std::istream& stream);
static void LoadOctreeHelper(Octree* node, const Grid& grid,
        std::istream& stream);
std::unique_ptr<Octree> LoadOctree(const Grid& grid, std::istream& stream);
std::unique_ptr<Octree> LoadOctree(const Grid& grid); // from stdin
std::unique_ptr<Octree> LoadOctree(const Grid& grid,
        const std::string& fileName);

}

}

#endif

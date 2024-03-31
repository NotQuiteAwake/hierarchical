#ifndef IOHEADERDEF
#define IOHEADERDEF

#include <iostream>
#include <string>
#include <memory>
#include "Octree.hpp"
#include "Grid.hpp"

namespace sim {

namespace IO {

void Err(const std::string& errMessage, std::ostream& stream = std::cerr);
bool FileExists(const std::string& fileName);
bool CheckFile(const std::string& fileName, bool expect);
void MakeDir(const std::string& dirName);

void SetHexfloatOut(std::ostream& stream);

Octant LoadOctant(std::istream& stream = std::cin);
Vec LoadVec(std::istream& stream = std::cin);
Particle LoadParticle(std::istream& stream = std::cin);

void DumpOctant(const Octant& octant, std::ostream& stream = std::cout);
void DumpVec(const Vec& vec, std::ostream& stream = std::cout);
void DumpParticle(const Particle& par, std::ostream& stream = std::cout);

void DumpGrid(const Grid& grid, std::ostream& stream = std::cout);
void DumpGrid(const Grid& grid, const std::string& file); // to file

Grid LoadGrid(std::istream& stream = std::cin);
Grid LoadGrid(const std::string& file); // from file

void DumpOctree(Octree const* octree, std::ostream& stream = std::cout);
void DumpOctree(Octree const* octree, const std::string& file); // to file

std::unique_ptr<Octree> LoadOctree(
        const Grid& grid,
        std::istream& stream = std::cin
        );
std::unique_ptr<Octree> LoadOctree(const Grid& grid,
        const std::string& fileName);

}

}

#endif

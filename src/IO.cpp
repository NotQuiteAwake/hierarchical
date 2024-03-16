// insert remarks about why dumping classes as json files in python is the best
// human invention ever :)

#include <string>
#include <cassert>
#include <fstream>
#include "IO.hpp"
#include "Octant.hpp"

namespace sim {

namespace IO {

template<typename T> static std::vector<T> LoadVector(std::istream& stream) {
    int cnt; stream >> cnt;
    std::vector<T> vector;
    for (int i = 0; i < cnt; i++) {
       T value; stream >> value;
       vector.push_back(value);
    }
    return vector;
}

template<typename T> static void DumpVector(const std::vector<T> vector,
        std::ostream& stream) {
    stream << vector.size();
    for (int i = 0; i < vector.size(); i++) {
        stream << " " << vector[i];
    }
    stream << std::endl;
}


Octant LoadOctant(std::istream& stream) { // TODO: does istream need to &
    Octant octant = Octant();
    for (int i = 0; i < Octant::mDim; i++) {
        stream >> octant[i][0] >> octant[i][1];
    }
    return octant;
}

Vec LoadVec(std::istream& stream) {
    Vec vec = Vec();
    for (int i = 0; i < Vec::mDim; i++) {
        stream >> vec[i];
    }
    return vec;
}

Particle LoadParticle(std::istream& stream) {
    double mass, charge;
    stream >> mass >> charge;
    Particle par = Particle(mass, charge);
    par.pos = LoadVec(stream);
    par.vel = LoadVec(stream);
    par.accel = LoadVec(stream);
    return par;
}

void DumpOctant(const Octant& octant, std::ostream& stream) {
    for (int i = 0; i < Octant::mDim; i++) {
        stream << octant[i][0] << " " << octant[i][1] << std::endl;
    }
}

void DumpVec(const Vec& vec, std::ostream& stream) {
    stream << vec[0];
    for (int i = 1; i < Vec::mDim; i++) {
        stream << " " << vec[i];
    }
    stream << std::endl;
}

void DumpParticle(const Particle& par, std::ostream& stream) {
    stream << par.GetMass() << " " << par.GetCharge() << std::endl;
    DumpVec(par.pos, stream);
    DumpVec(par.vel, stream);
    DumpVec(par.accel, stream);
}

void DumpParticle(const Particle& par) {
    DumpParticle(par, std::cout);
}

void DumpGrid(const Grid& grid, std::ostream& stream) {
    stream << "BEGIN GRID" << std::endl;

    // print limits to grid
    DumpOctant(grid.GetLimits(), stream);
    
    // print all particles in grid
    stream << grid.GetSize() << std::endl;
    for (int i = 0; i < grid.GetSize(); i++) {
        DumpParticle(grid[i], stream);
    }

    stream << "END GRID" << std::endl;
    
}

void DumpGrid(const Grid& grid) {
    DumpGrid(grid, std::cout);
}

void DumpGrid(const Grid& grid, const std::string& fileName) {
    std::ofstream stream;
    stream.open(fileName, std::ios::out); // TODO: ios::noreplace?
    assert(stream.good());
    DumpGrid(grid, stream);
    stream.close();
    
}

Grid LoadGrid(std::istream& stream) {
    std::string line;
    stream >> line;
    assert(line == "BEGIN GRID");
        
    Grid grid = Grid(LoadOctant(stream));
    int size;
    stream >> size;
    for (int i = 0; i < size; i++) {
        grid.AddParticle(LoadParticle(stream));
    }
    
    stream >> line;
    assert(line == "END GRID");
    
    return grid;
}

Grid LoadGrid() {
    return LoadGrid(std::cin);
}

Grid LoadGrid(const Grid& grid, const std::string& fileName) {
    std::ifstream stream;
    stream.open(fileName, std::ios::in);
    assert(stream.good()); // TODO: CHECK THIS.
    return LoadGrid(stream);
}

static void DumpOctreeNode(Octree const* node, std::ostream& stream) {
    DumpOctant(node->octant, stream);
    DumpVec(node->com, stream);
    stream << node->mass << std::endl;

    DumpVector(node->mSouls, stream);
    DumpVector(node->M, stream);
    DumpVector(node->L, stream);
}

static void DumpOctreeHelper(Octree const* node, std::ostream& stream) {
    assert(node);
    DumpOctreeNode(node, stream);

    // the below always guarantee that on reconstruction we can read a number
    if (node->IsLeaf()) {
        stream << -1 << std::endl;
        return;
    }
    for (int i = 0; i < node->mBoxes; i++) {
        Octree const* child_node = node->GetChild(i);
        if (!child_node) continue;
        stream << i << std::endl;
        DumpOctreeHelper(child_node, stream);
    } 
    stream << -1 << std::endl;
}

void DumpOctree(Octree const* octree, std::ostream& stream) {
    stream << "BEGIN OCTREE" << std::endl;
    stream << octree->GetMaxParticles() << std::endl;
    DumpOctreeHelper(octree, stream);
    stream << "END OCTREE" << std::endl;
}

void DumpOctree(Octree const* octree) {
    DumpOctree(octree, std::cout);
}

void DumpOctree(Octree const* octree, const std::string& fileName) {
    std::ofstream stream;
    stream.open(fileName, std::ios::out);
    assert(stream.good());
    DumpOctree(octree, stream);
}

static void LoadOctreeNode(Octree* node, std::istream& stream) {
    node->octant = LoadOctant(stream);
    node->com = LoadVec(stream);
    stream >> node->mass;
    
    node->mSouls = LoadVector<int>(stream);
    // mOctantSouls is not recovered as a compromise for readability.
    // it is irrelevant for any complete tree.
    node->M = LoadVector<double>(stream);
    node->L = LoadVector<double>(stream);
}

static void LoadOctreeHelper(Octree* node, const Grid& grid, std::istream& stream) {
    LoadOctreeNode(node, stream);
    int child_octant_num; stream >> child_octant_num;
    while (child_octant_num != -1) {
        node->SetChild(child_octant_num,
                new Octree(node, grid, Octant(), node->GetMaxParticles()));
        LoadOctreeNode(node->GetChild(child_octant_num), stream);
        stream >> child_octant_num;
    } 
}

std::unique_ptr<Octree> LoadOctree(const Grid& grid, std::istream& stream) {
    std::string line;
    stream >> line;
    assert(line == "BEGIN OCTREE");

    int max_particles;
    stream >> max_particles;

    std::unique_ptr<Octree> root; 
    root.reset(new Octree(nullptr, grid, Octant(), max_particles));
    LoadOctreeHelper(root.get(), grid, stream);

    stream >> line;
    assert(line == "END OCTREE");
    
    return root;
}

std::unique_ptr<Octree> LoadOctree(const Grid& grid) {
    return LoadOctree(grid, std::cin);
}

std::unique_ptr<Octree> LoadOctree(const Grid& grid,
        const std::string& fileName) {
    std::ifstream stream;
    stream.open(fileName);
    assert(stream.good());
    return LoadOctree(grid, stream);

}

}

}

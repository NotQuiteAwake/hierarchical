/**
 * @file
 * @brief Library for manipulating data into / from text.
 *
 * The function names and parameters should be self explanatory. When an
 * std::istream is not supplied, std::cin is used by default. Similarly if
 * std::ostream is absent, std::cout is used. Alternatively a path to a file can
 * be supplied with the fileName parameter.
 *
 * Insert remarks about why dumping classes as json files in python is the best
 * human invention ever :)
 */

#include <string>
#include <cassert>
#include <fstream>
#include <filesystem>
#include "IO.hpp"
#include "Octant.hpp"

namespace sim {

namespace IO {

namespace {
template<typename T> std::vector<T> LoadVector(std::istream& stream) {
    int cnt; stream >> cnt;
    std::vector<T> vector;
    for (int i = 0; i < cnt; i++) {
       T value; stream >> value;
       vector.push_back(value);
    }
    return vector;
}

template<typename T> Matrix<T> LoadMatrix(std::istream& stream) {
    int nrows, ncols;
    stream >> nrows >> ncols;
    Matrix<T> matrix(nrows, ncols);
    for (int i = 0; i <= nrows; i++) {
        // note our convention
        for (int j = -ncols; j <= ncols; j++) {
            stream >> matrix[i][j];
        }
    }
    return matrix;
}

template<typename T> void DumpVector(
        const std::vector<T> vector,
        std::ostream& stream
        ) {
    stream << vector.size();
    for (int i = 0; i < vector.size(); i++) {
        stream << " " << vector[i];
    }
    stream << std::endl;
}

template<typename T> void DumpMatrix(
        const Matrix<T> matrix,
        std::ostream& stream
        ) {
    int nrows = matrix.GetRows();
    int ncols = matrix.GetCols();
    stream << nrows << " " << ncols << std::endl;
    for (int i = 0; i <= nrows; i++) {
        for (int j = -ncols; j <= ncols; j++) {
            stream << matrix[i][j] << " ";
        }
        stream << std::endl;
    }
}
}

void Err(const std::string& errMessage, std::ostream& stream) {
    stream << "(E) " << errMessage << std::endl;
}

bool FileExists(const std::string& fileName) {
    // fileName automatically converted to std::filesystem::path
    return std::filesystem::exists(fileName);
}

/**
 * @brief Check if file or folder exists.
 *
 * If the existence state is different from the expect value, complain in log.
 * The function that receives the unexpected state will likely just abort.
 *
 * @param[in] fileName Path to file or folder.
 * @param[in] expect Whether we expect the file or folder to be there or not
 * @return Whether the folder or file exists or not.
 */
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

void MakeDir(const std::string& dirName) {
    std::filesystem::create_directory(dirName);
}

/**
 * @brief Set the ostream to output floating points in hexfloat format.
 *
 * Note that for istream's, the use of std::hexfloat seems to be less
 * straightforward. There appears to be a bug in g++:
 *
 * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=81122
 *
 * However my attemps in clang++ has been equally unsucessful.
 *
 * @param[in, out] stream Output stream to be acted on.
 */
void SetHexfloatOut(std::ostream& stream) {
    stream << std::hexfloat;
}

Octant LoadOctant(std::istream& stream) {
    bool is_initialised;
    stream >> is_initialised;
    if (!is_initialised) return Octant();

    double lim[Octant::mDim][2];
    for (int i = 0; i < Octant::mDim; i++) {
        stream >> lim[i][0] >> lim[i][1];
    }
    return Octant(lim);
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
    stream << octant.IsInitialised() << std::endl;
    if (!octant.IsInitialised()) return;
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

void DumpGrid(const Grid& grid, std::ostream& stream) {
    stream << "BEGIN GRID" << std::endl;

    // print limits to grid.
    DumpOctant(grid.GetLimits(), stream);
    DumpOctant(grid.GetOctant(), stream);
    
    // print all particles in grid
    stream << grid.GetSize() << std::endl;
    for (int i = 0; i < grid.GetSize(); i++) {
        DumpParticle(grid[i], stream);
    }

    stream << "END GRID" << std::endl;
    
}

void DumpGrid(const Grid& grid, const std::string& fileName) {
    std::ofstream stream;
    stream.open(fileName, std::ios::out); // TODO: ios::noreplace?
    assert(stream.good());
    DumpGrid(grid, stream);
    stream.close();
    
}

Grid LoadGrid(std::istream& stream) {
    // TODO: maintain maxLim property
    std::string flag, token;
    stream >> flag >> token;
    assert(flag == "BEGIN" && token == "GRID");
        
    Grid grid = Grid(LoadOctant(stream));
    grid.SetOctant(LoadOctant(stream));
    int size;
    stream >> size;
    for (int i = 0; i < size; i++) {
        grid.AddParticle(LoadParticle(stream));
    }
    
    stream >> flag >> token;
    assert(flag == "END" && token == "GRID");
    
    return grid;
}

Grid LoadGrid(const Grid& grid, const std::string& fileName) {
    std::ifstream stream;
    stream.open(fileName, std::ios::in);
    assert(stream.good()); // TODO: CHECK THIS.
    return LoadGrid(stream);
}

namespace {
void DumpOctreeNode(Octree const* node, std::ostream& stream) {
    DumpOctant(node->GetOctant(), stream);
    DumpVec(node->coc, stream);
    stream << node->charge << std::endl;

    DumpVector(node->mSouls, stream);
    DumpMatrix(node->M, stream);
    DumpMatrix(node->F, stream);
}

/**
 * @brief Dump the Octree recursively in a pre-order DFS
 *
 */
void DumpOctreeHelper(Octree const* node, std::ostream& stream) {
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
}

void DumpOctree(Octree const* root, std::ostream& stream) {
    stream << "BEGIN OCTREE" << std::endl;
    stream << root->GetMaxParticles() << " ";
    stream << root->GetP() << std::endl;
    DumpOctreeHelper(root, stream);
    stream << "END OCTREE" << std::endl;
}

void DumpOctree(Octree const* octree, const std::string& fileName) {
    std::ofstream stream;
    stream.open(fileName, std::ios::out);
    assert(stream.good());
    DumpOctree(octree, stream);
}

namespace {
void LoadOctreeNode(Octree* node, std::istream& stream) {
    typedef std::complex<double> cdouble;
    node->SetOctant(LoadOctant(stream));
    node->coc = LoadVec(stream);
    stream >> node->charge;
    
    node->mSouls = LoadVector<int>(stream);
    // mOctantSouls is not recovered as a compromise for readability.
    // it is irrelevant for any complete tree.
    node->M = LoadMatrix<cdouble>(stream);
    node->F = LoadMatrix<cdouble>(stream);
}

/**
 * @brief Load Octree recursively in a pre-order fashion
 *
 */
void LoadOctreeHelper(Octree* node, const Grid& grid, std::istream& stream) {
    LoadOctreeNode(node, stream);
    int child_octant_num; stream >> child_octant_num;
    while (child_octant_num != -1) {
        node->SetChild(child_octant_num,
                new Octree(node,
                    grid,
                    Octant(),
                    node->GetMaxParticles(),
                    node->GetP())
                );
        LoadOctreeHelper(node->GetChild(child_octant_num), grid, stream);
        stream >> child_octant_num;
    } 
}
}

std::unique_ptr<Octree> LoadOctree(const Grid& grid, std::istream& stream) {
    std::string flag, token;
    stream >> flag >> token;
    assert(flag == "BEGIN" && token == "OCTREE");

    int max_particles, p;
    stream >> max_particles >> p;

    std::unique_ptr<Octree> root; 
    root.reset(new Octree(nullptr, grid, Octant(), max_particles, p));
    LoadOctreeHelper(root.get(), grid, stream);

    stream >> flag >> token;
    assert(flag == "END" && token == "OCTREE");
    
    return root;
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

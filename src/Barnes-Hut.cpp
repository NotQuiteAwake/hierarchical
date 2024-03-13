#include "Barnes-Hut.hpp"

namespace sim {

double BarnesHut::SumChildM(Octree* node) const {
    double M = 0;
    for (int i = 0; i < node->mBoxes; i++) {
        Octree* child_node = node->GetChild(i);
        if (child_node) {
            M += child_node->M[0];
        }
    }
    return M;
}

double BarnesHut::SumChildMass(Octree* node) const {
    double mass = 0;
    for (int i = 0; i < node->mBoxes; i++) {
        Octree* child_node = node->GetChild(i);
        if (child_node) {
            mass += child_node->mass;
        }
    }
    return mass;
}

Vec BarnesHut::GetNodeCOM(Octree* node) const {
    Vec new_com = Vec();
    for (int i = 0; i < node->mBoxes; i++) {
        const Octree* child_node = node->GetChild(i);
        if (child_node) {
            new_com += child_node->com * child_node->M[0];
        }
    }
    new_com = new_com / SumChildM(node);
    return new_com;
}

void BarnesHut::ProcessTree(Octree* node) const {
    if (!node) { return; }
    node->M.reset(new double[mP]);

    if (node -> IsLeaf()) {
        int soul = node->mSouls[0];
        const Particle& par = node->GetParticle(soul);
        node->M[0] = par.GetCharge(); // only grav charge is mass.
        node->com = par.GetPos();
        node->mass = par.GetMass();
    } else {
        node->M[0] = SumChildM(node);
        node->com = GetNodeCOM(node);
        node->mass = SumChildMass(node);
    }
}

std::unique_ptr<Octree> BarnesHut::BuildTree(const Grid& g1) const {
    std::unique_ptr<Octree> octree(new Octree(nullptr, g1, g1.GetLimits(), 1));
    for (int i = 0; i < g1.GetSize(); i++) {
        const Particle& par = g1[i];
        octree -> AddParticle(par, i);
    }

    return octree;
}

bool BarnesHut::IsFarAway(const Particle& par, const Octree* const node) const {
    const Vec dr = par.GetPos() - node->com;
    const double dist = dr.CalculateNorm();
    double l = node->GetMaxLength();
    return l < dist * mTheta; // for better double precision
}

Vec BarnesHut::GetAccel(int soul, const Octree* const node) const {
    const Particle& par = node->GetParticle(soul);

    if (!node) { return Vec(); }
    if (node->IsLeaf()) {
        if (node->mSouls[0] == soul) { return Vec(); }
        else {
            const Particle& par2 = node->GetParticle(node->mSouls[0]);
            return GetForce(par, par2);
        }
    }

    if (IsFarAway(par, node)) {
        Particle virt_par = Particle(node->mass, node->M[0]);
        virt_par.SetPos(node->com);
        return GetForce(par, virt_par) / par.GetMass();
    }

    Vec accel = Vec();
    for (int i = 0; i < node->mBoxes; i++) {
        const Octree* child_node = node->GetChild(i);
        accel += GetAccel(soul, child_node);
    }
    return accel;
}

BarnesHut::BarnesHut(int p, double theta,
        const std::shared_ptr<const Force> forceLaw) :
    Interaction(forceLaw),
    mP(p),
    mTheta(theta) {

}

int BarnesHut::GetP() const {
    return mP;
}

double BarnesHut::GetTheta() const {
    return mTheta;
}

Grid BarnesHut::Calculate(const Grid& g1) const {
    Grid g2 = g1;
    std::unique_ptr<Octree> octree = BuildTree(g1);
    ProcessTree(octree.get());

    for (int i = 0; i < g1.GetSize(); i++) {
        const Particle& par = g1[i];
        g2[i].SetAccel(GetAccel(i, octree.get()));
    }
}

}

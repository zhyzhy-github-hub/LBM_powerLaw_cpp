// -------------By ZhangYu HIT
// -------------1225839327@qq.com
#ifndef _NSEQULBM_H_
#define _NSEQULBM_H_

#include <iostream>
#include <cmath>
#include "lattice.h"
#include "mesh.h"

using std::cout;
using std::endl;

class nseD2Q9 : public Lattice::D2Q9Lattice
{
public:
    nseD2Q9(MeshBasic *mesh_, const double dx_, const double dt_, const double nu_);
    ~nseD2Q9();
    void initScalarAndVector();
    void initScalarAndVector(const double tau_);
    void initFeq();

    void streaming();
    void collision();
    void boundaryCondition(const Vector2D u);
    void outputTec(double t) const;

    void collisionNonNewton();
    stressTensor computeStressTensor();
    void computeStressTensor(stressTensor &D, const distributionFunction &fStr,
                             const scalar &rho_, const Vector2D &u_, const scalar tau_);

    scalar getGamma(const stressTensor D) const { return 1 * sqrt(1 * D.DProductD()); };
    scalar getTau(const double gamma, const double m_, const double n_, const double rho_) const;
    void getTauAll(const double m_, const double n_);

    void evolutionCollisionAndStream(const Vector2D u);
    void evolutionCollisionAndStreamNonNewton(const Vector2D u, const double m_, const double n_);

    size_t getScalarCoordinate(const int i, const int j) const { return mesh->scalarIndex(i, j); }
    scalar computeAllOfrho() const;
    Vector2D computeAllOfMomentum() const;
    Vector2D computeAllOfVelocity() const;

    Vector2D getVelocity(const int i, const int j) const { return velocity[getScalarCoordinate(i, j)]; }

private:
    double dx, dt, c;
    double nu, tau;

    MeshBasic *mesh;
    // static const MeshBasic *mesh;
    scalar *rho;
    scalar *tauNonNewton;
    Vector2D *velocity;
    Vector2D *force;
    distributionFunction *funCollision;
    distributionFunction *funStreaming;

private:
    double feq(int k, scalar rho, Vector2D u) const;
    void collisionKernel(distributionFunction &fCol, const distributionFunction &fStr, const scalar &rho_, const Vector2D &u_);
    void collisionKernel(distributionFunction &fCol, const distributionFunction &fStr, const scalar &rho_, const Vector2D &u_, const scalar &tau);
    double getFneq(const distributionFunction &f, const int k, const scalar T, const Vector2D u) const { return (f.f[k] - feq(k, T, u)); }
    scalar getFneq(const int i, const int j, const int k) const;
    double nonEquilibriumExtraPolation(int k, scalar rho, Vector2D &u, scalar fneq) const { return (feq(k, rho, u) - fneq); }
};

#endif
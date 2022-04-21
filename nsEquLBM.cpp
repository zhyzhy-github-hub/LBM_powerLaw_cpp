// -------------By ZhangYu HIT
// -------------1225839327@qq.com
#include <fstream>
#include <sstream>
#include "nsEquLBM.h"

nseD2Q9::nseD2Q9(MeshBasic *mesh_, const double dx_, const double dt_, const double nu_)
    : dx(dx_), dt(dt_), nu(nu_)
{
    MeshBasic *meshTemp = mesh_;
    mesh = meshTemp;

    c = dx / dt;
    tau = 0.5 + nu / (Cs2 * c * c * dt);

    const int NX = mesh->getNX();
    const int NY = mesh->getNY();
    const int NN = NX * NY;

    rho = new scalar[NN]();
    tauNonNewton = new scalar[NN]();
    velocity = new Vector2D[NN]();
    funCollision = new distributionFunction[NN]();
    funStreaming = new distributionFunction[NN]();
    force = new Vector2D[NN]();

    cout << "c = " << c << endl;
    cout << "tau = " << tau << endl;
}

nseD2Q9::~nseD2Q9()
{
    delete[] rho;
    delete[] tauNonNewton;
    delete[] velocity;
    delete[] funCollision;
    delete[] funStreaming;
    delete[] force;

    cout << "-----------------------------NSE---------------------------" << endl;
    cout << "Free memory allocated by D2Q9." << endl;
}
void nseD2Q9::initScalarAndVector(const double tau_)
{
    const int NX = mesh->getNX();
    const int NY = mesh->getNY();
    for (size_t j = 0; j < NY; ++j)
    {
        for (size_t i = 0; i < NX; ++i)
        {
            size_t index = getScalarCoordinate(i, j);

            rho[index] = 1.0;
            tauNonNewton[index] = tau;
            velocity[index] = {0.0, 0.0};
            force[index] = {0.0, 0.0};
        }
    }
}

void nseD2Q9::initScalarAndVector()
{
    const int NX = mesh->getNX();
    const int NY = mesh->getNY();
    for (size_t j = 0; j < NY; ++j)
    {
        for (size_t i = 0; i < NX; ++i)
        {
            size_t index = getScalarCoordinate(i, j);

            rho[index] = 1.0;
            tauNonNewton[index] = tau;
            // if (j > 0.45 * NY && j < 0.55 * NY && i > 0.45 * NY && i < 0.55 * NY)
            // {
            //     rho[index] = 1.1;
            // }
            velocity[index] = {0.0, 0.0};
            force[index] = {0.0, 0.0};
        }
    }
}

void nseD2Q9::initFeq()
{
    const int NX = mesh->getNX();
    const int NY = mesh->getNY();
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = getScalarCoordinate(i, j);
            scalar rho_ = rho[index];
            Vector2D u_ = velocity[index];
            for (int k = 0; k < Q; ++k)
            {
                double feqQ = feq(k, rho_, u_);
                funStreaming[index].f[k] = feqQ;
                funCollision[index].f[k] = feqQ;
            }
        }
    }
}

double nseD2Q9::feq(int k, scalar rho, Vector2D u) const
{
    double feq;
    double CU_Cs2inv;
    const double Cs2_inv = 1.0 / (Cs2 * c * c);
    const double ux = u.x;
    const double uy = u.y;
    const double Omusq = 1.0 - (ux * ux + uy * uy) * 0.5 * Cs2_inv;
    switch (k)
    {
    case 0:
        feq = w0 * rho * Omusq;
        break;
    case 1:
        CU_Cs2inv = ux * Cs2_inv * c;
        feq = wc * rho * (CU_Cs2inv + 0.5 * (CU_Cs2inv) * (CU_Cs2inv) + Omusq);
        break;
    case 2:
        CU_Cs2inv = uy * Cs2_inv * c;
        feq = wc * rho * (CU_Cs2inv + 0.5 * (CU_Cs2inv) * (CU_Cs2inv) + Omusq);
        break;
    case 3:
        CU_Cs2inv = -ux * Cs2_inv * c;
        feq = wc * rho * (CU_Cs2inv + 0.5 * (CU_Cs2inv) * (CU_Cs2inv) + Omusq);
        break;
    case 4:
        CU_Cs2inv = -uy * Cs2_inv * c;
        feq = wc * rho * (CU_Cs2inv + 0.5 * (CU_Cs2inv) * (CU_Cs2inv) + Omusq);
        break;
    case 5:
        CU_Cs2inv = (ux + uy) * Cs2_inv * c;
        feq = wd * rho * (CU_Cs2inv + 0.5 * (CU_Cs2inv) * (CU_Cs2inv) + Omusq);
        break;
    case 6:
        CU_Cs2inv = (-ux + uy) * Cs2_inv * c;
        feq = wd * rho * (CU_Cs2inv + 0.5 * (CU_Cs2inv) * (CU_Cs2inv) + Omusq);
        break;
    case 7:
        CU_Cs2inv = (-ux - uy) * Cs2_inv * c;
        feq = wd * rho * (CU_Cs2inv + 0.5 * (CU_Cs2inv) * (CU_Cs2inv) + Omusq);
        break;
    case 8:
        CU_Cs2inv = (ux - uy) * Cs2_inv * c;
        feq = wd * rho * (CU_Cs2inv + 0.5 * (CU_Cs2inv) * (CU_Cs2inv) + Omusq);
        break;
    }
    return feq;
}

void nseD2Q9::collisionKernel(distributionFunction &fCol, const distributionFunction &fStr, const scalar &rho_, const Vector2D &u_, const scalar &tau_)
{
    // const double tau_ = 0.5 + nu_ / (Cs2 * c * c * dt);
    const double tauInv = 1.0 / tau_;
    const double oneMtauInv = 1.0 - tauInv;
    for (int k = 0; k < Q; ++k)
    {
        fCol.f[k] = oneMtauInv * fStr.f[k] + tauInv * feq(k, rho_, u_);
    }
}

void nseD2Q9::collisionKernel(distributionFunction &fCol, const distributionFunction &fStr, const scalar &rho_, const Vector2D &u_)
{
    const double tauInv = 1.0 / tau;
    const double oneMtauInv = 1.0 - 1.0 / tau;
    for (int k = 0; k < Q; ++k)
    {
        fCol.f[k] = oneMtauInv * fStr.f[k] + tauInv * feq(k, rho_, u_);
    }
}
void nseD2Q9::collisionNonNewton()
{
    const int NX = mesh->getNX();
    const int NY = mesh->getNY();
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = getScalarCoordinate(i, j);
            collisionKernel(funCollision[index], funStreaming[index], rho[index], velocity[index], tauNonNewton[index]);
        }
    }
}
void nseD2Q9::collision()
{
    const int NX = mesh->getNX();
    const int NY = mesh->getNY();
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = getScalarCoordinate(i, j);
            collisionKernel(funCollision[index], funStreaming[index], rho[index], velocity[index]);
        }
    }
}
void nseD2Q9::computeStressTensor(stressTensor &D, const distributionFunction &fStr,
                                  const scalar &rho_, const Vector2D &u_, const scalar tau_)
{
    const double fneq1 = getFneq(fStr, 1, rho_, u_);
    const double fneq2 = getFneq(fStr, 2, rho_, u_);
    const double fneq3 = getFneq(fStr, 3, rho_, u_);
    const double fneq4 = getFneq(fStr, 4, rho_, u_);
    const double fneq5 = getFneq(fStr, 5, rho_, u_);
    const double fneq6 = getFneq(fStr, 6, rho_, u_);
    const double fneq7 = getFneq(fStr, 7, rho_, u_);
    const double fneq8 = getFneq(fStr, 8, rho_, u_);
    stressTensor D_temp{
        fneq1 + fneq3 + fneq5 + fneq6 + fneq7 + fneq8,
        fneq2 + fneq4 + fneq5 + fneq6 + fneq7 + fneq8,
        fneq5 - fneq6 + fneq7 - fneq8};

    D_temp *= (-1.0 / (2 * tau_ * rho_ * Cs2 * c * c * dt));

    D = D_temp;
}
scalar nseD2Q9::getFneq(const int i, const int j, const int k) const
{
    const size_t index = getScalarCoordinate(i, j);
    double feq_ = feq(k, rho[index], velocity[index]);
    return (funStreaming[index].f[k] - feq_);
}

scalar nseD2Q9::getTau(const double gamma, const double m_, const double n_, const double rho_) const
{
    double nM1 = n_ - 1;
    double mu_ = m_ * pow(gamma, nM1);

    double nu_ = mu_ / rho_;

    double tau = 0.5 + nu_ / (Cs2 * c * c * dt);
    return tau;

    // if (nu_ > this->nu * 20)
    //    nu_ = 20 * nu;
    // if (nu_ < this->nu * 0.05)
    //    nu_ = this->nu * 0.05;
    // if (nu_ < 0.001)
    //     nu_ = 0.001;
    // if (nu_ > 1.0)
    //     nu_ = 1.0;
}

void nseD2Q9::getTauAll(const double m_, const double n_)
{
    const int NX = mesh->getNX();
    const int NY = mesh->getNY();
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            const size_t index = getScalarCoordinate(i, j);

            stressTensor DTemp;
            computeStressTensor(DTemp, funStreaming[index], rho[index], velocity[index], tauNonNewton[index]);
            double gammaTemp = getGamma(DTemp);
            tauNonNewton[index] = getTau(gammaTemp, m_, n_, rho[index]);
        }
    }
}

void nseD2Q9::streaming()
{
    const int NX = mesh->getNX();
    const int NY = mesh->getNY();
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            const size_t index = getScalarCoordinate(i, j);
            unsigned int ip1 = (i + 1) % (NX);
            unsigned int jp1 = (j + 1) % (NY);
            unsigned int im1 = (NX + i - 1) % (NX);
            unsigned int jm1 = (NY + j - 1) % (NY);

            double ft0 = funCollision[getScalarCoordinate(i, j)].f[0];
            double ft1 = funCollision[getScalarCoordinate(im1, j)].f[1];
            double ft2 = funCollision[getScalarCoordinate(i, jm1)].f[2];
            double ft3 = funCollision[getScalarCoordinate(ip1, j)].f[3];
            double ft4 = funCollision[getScalarCoordinate(i, jp1)].f[4];
            double ft5 = funCollision[getScalarCoordinate(im1, jm1)].f[5];
            double ft6 = funCollision[getScalarCoordinate(ip1, jm1)].f[6];
            double ft7 = funCollision[getScalarCoordinate(ip1, jp1)].f[7];
            double ft8 = funCollision[getScalarCoordinate(im1, jp1)].f[8];

            funStreaming[index].f[0] = ft0;
            funStreaming[index].f[1] = ft1;
            funStreaming[index].f[2] = ft2;
            funStreaming[index].f[3] = ft3;
            funStreaming[index].f[4] = ft4;
            funStreaming[index].f[5] = ft5;
            funStreaming[index].f[6] = ft6;
            funStreaming[index].f[7] = ft7;
            funStreaming[index].f[8] = ft8;

            //-------------------------------- computing macro valiables
            double rho_ = ft0                     //
                          + ft1 + ft2 + ft3 + ft4 //
                          + ft5 + ft6 + ft7 + ft8;

            double rhoinv = 1.0 / rho_;
            double ux_ = rhoinv * (ft1 + ft5 + ft8 - (ft3 + ft6 + ft7)) * c;
            double uy_ = rhoinv * (ft2 + ft5 + ft6 - (ft4 + ft7 + ft8)) * c;

            rho[index] = rho_;

            double force_x_ = force[index].x;
            double force_y_ = force[index].y;

            velocity[index].x = ux_ + force_x_ * 0.5 * rhoinv;
            velocity[index].y = uy_ + force_y_ * 0.5 * rhoinv;
            // ux[index] = ux_ + force_x_ * dt * 0.5 * rhoinv;
            // uy[index] = uy_ + force_y_ * dt * 0.5 * rhoinv;
        }
    }
}
void nseD2Q9::evolutionCollisionAndStreamNonNewton(const Vector2D u, const double m_, const double n_)
{
    streaming();
    boundaryCondition(u);
    collisionNonNewton();
    getTauAll(m_, n_);
}

void nseD2Q9::evolutionCollisionAndStream(const Vector2D u)
{
    streaming();
    boundaryCondition(u);
    collision();
}

Vector2D nseD2Q9::computeAllOfVelocity() const
{
    const int NX = mesh->getNX();
    const int NY = mesh->getNY();
    Vector2D velocityAll{0, 0};
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = getScalarCoordinate(i, j);
            velocityAll.x += velocity[index].x;
            velocityAll.y += velocity[index].y;
        }
    }
    return velocityAll;
}
Vector2D nseD2Q9::computeAllOfMomentum() const
{
    const int NX = mesh->getNX();
    const int NY = mesh->getNY();
    Vector2D momentumAll{0, 0};
    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = getScalarCoordinate(i, j);
            momentumAll.x += rho[index] * velocity[index].x;
            momentumAll.y += rho[index] * velocity[index].y;
        }
    }
    return momentumAll;
}
scalar nseD2Q9::computeAllOfrho() const
{
    const int NX = mesh->getNX();
    const int NY = mesh->getNY();
    scalar all = 0;

    for (int j = 0; j < NY; ++j)
    {
        for (int i = 0; i < NX; ++i)
        {
            size_t index = getScalarCoordinate(i, j);
            all += rho[index];
        }
    }
    return all;
}

void nseD2Q9::boundaryCondition(const Vector2D u)
{
    const int NX = mesh->getNX();
    const int NY = mesh->getNY();

    for (int j = 1; j < NY - 1; ++j)
    {
        size_t left = getScalarCoordinate(0, j);
        size_t leftP1 = getScalarCoordinate(1, j);
        // velocity[left] = {0,0.1};
        velocity[left] = {0, 0};
        rho[left] = rho[leftP1];
        // rho[left] = 1.0;
        for (int k = 0; k < Q; ++k)
        {
            funStreaming[left].f[k] = feq(k, rho[left], velocity[left]) + getFneq(funStreaming[leftP1], k, rho[leftP1], velocity[leftP1]);
        }

        size_t right = getScalarCoordinate(NX - 1, j);
        size_t rightM1 = getScalarCoordinate(NX - 2, j);
        // velocity[right] = {velocity[rightM1].x, 0};
        velocity[right] = {0, 0};
        rho[right] = rho[rightM1];
        // rho[right] = 1.0;
        for (int k = 0; k < Q; ++k)
        {
            funStreaming[right].f[k] = feq(k, rho[right], velocity[right]) + getFneq(funStreaming[rightM1], k, rho[rightM1], velocity[rightM1]);
        }
    }
    for (int i = 0; i < NX - 0; ++i)
    {
        size_t topM0 = getScalarCoordinate(i, NY - 1);
        size_t topM1 = getScalarCoordinate(i, NY - 2);

        Vector2D uTop{u};
        velocity[topM0] = uTop;
        rho[topM0] = rho[topM1];

        for (int k = 0; k < Q; ++k)
        {
            funStreaming[topM0].f[k] = feq(k, rho[topM0], uTop) + getFneq(funStreaming[topM1], k, rho[topM1], velocity[topM1]);
        }

        size_t bottomP0 = getScalarCoordinate(i, 0);
        size_t bottomP1 = getScalarCoordinate(i, 1);
        Vector2D uBottom{0, 0};
        velocity[bottomP0] = uBottom;
        rho[bottomP0] = rho[bottomP1];

        for (int k = 0; k < Q; ++k)
        {
            funStreaming[bottomP0].f[k] = feq(k, rho[bottomP0], uBottom) + getFneq(funStreaming[bottomP1], k, rho[bottomP1], velocity[bottomP1]);
        }
    }
}

void nseD2Q9::outputTec(double t) const
{
    const int NX = mesh->getNX();
    const int NY = mesh->getNY();
    std::ostringstream name;
    name << "NSE_t_" << std::setw(8) << std::setfill('0') << t << "_.dat";
    std::ofstream out(filePath + name.str().c_str());
    out << "Title= \"Lid Driven Flow\"\n";
    out << "VARIABLES = \"X\", \"Y\", \"rho\", \"ux\", \"uy\", \"forcex\", \"forcey\", \"tau\" \n";
    out << "ZONE T= \"BOX\",I=" << NX << ",J=" << NY << ",F=	POINT" << std::endl;

    for (size_t j = 0; j < NY; ++j)
    {
        for (size_t i = 0; i < NX; ++i)
        {
            const size_t index = getScalarCoordinate(i, j);
            const Coordinate coord = mesh->getCoordinate(i, j);
            out << std::fixed << std::setprecision(10)
                << " " << static_cast<double>(i) / (NX - 1) << " " << static_cast<double>(j) / (NX - 1)
                // << " " << coord.x << " " << coord.y
                << " " << rho[index]
                << " " << velocity[index].x
                << " " << velocity[index].y
                << " " << force[index].x
                << " " << force[index].y
                << " " << tauNonNewton[index]
                // << " " << tau
                << std::endl;
        }
    }
    out.close();
}
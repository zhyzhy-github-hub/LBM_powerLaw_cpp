// -------------By ZhangYu HIT
// -------------1225839327@qq.com
#include <fstream>
#include <sstream>
#include "mesh.h"
#include "nsEquLBM.h"

int main(int argc, char **argv)
{

    // AirPhysicalParameters::printAirPhysicalParameters();
    // LbParameters::printLBParameters();

    // const int NY = LbParameters::NY + 1;
    // const int NX = (NY - 1) * 5 + 1;

    const int NY = 128 * 2 + 1;
    const int NX = (NY - 1 * 1) + 1;


    double nNu = 0.5;
    cout << "Please input nNu = :" << endl;
    std::cin >> nNu;
    // double nNu = 1.0 - 0.70;

    const scalar uTop = 0.1;
    double Re = 400;
    Vector2D uIn{uTop, 0};
    double nu = pow(uTop, 2.0 - nNu) * pow(NY - 1, nNu) / Re;
    // double nu = pow(uTop, 2.0 - nNu) * pow(NY - 1, nNu) / Re;
    double m = nu;

    cout << "Re = " << Re
         << ", uTop = " << uTop
         << ", nNu = " << nNu
         << ", nu = " << nu
         << endl;

    MeshBasic mesh0(NX, NY);
    mesh0.initMesh();
    mesh0.outputMeshTec();

    nseD2Q9 fluid0(&mesh0, 1.0, 1.0, nu);

    // fluid0.initScalarAndVector(1.0 / 1.5);
    fluid0.initScalarAndVector();
    fluid0.initFeq();
    fluid0.outputTec(0);

    for (int n1 = 0; n1 < 1000000; ++n1)
    {
        int n = n1 + 1;
        fluid0.evolutionCollisionAndStreamNonNewton(uIn, nu, nNu);

        if ((n) % 5000 == 0)
        {

            Vector2D mAll = fluid0.computeAllOfMomentum();
            Vector2D velocityAll = fluid0.computeAllOfVelocity();
            cout << " n = " << n
                 << " , all of rho= " << fluid0.computeAllOfrho()
                 << " , all of momentum = (" << mAll.x << ", " << mAll.y << ") "
                 << " , all of velocity = (" << velocityAll.x << ", " << velocityAll.y << ") "
                 << endl;

            fluid0.outputTec(n);
        }
    }
}

//void output_Tec_result(const double t,
//                       nseD2Q9 &nse, adeD2Q5& ade)
//{
//
//    std::ostringstream name;
//    name << "Natural_convection_" << t << ".dat";
//    std::ofstream out(name.str().c_str());
//    out << "Title= \"NC cylinder Flow\"\n";
//    out << "VARIABLES = \"X\", \"Y\", \"rho\", \"ux\", \"uy\", \"C\"  \n";
//    out << "ZONE T= \"BOX\",I=" << NX << ",J=" << NY << ",F=	POINT" << std::endl;
//
//    for (size_t j = 0; j < NY; ++j)
//    {
//        for (size_t i = 0; i < NX; ++i)
//        {
//            const size_t index = mesh0.get_scalar_index(i, j);
//            out << std::fixed << std::setprecision(16)
//                << " " << mesh0.xNorm[index] << " " << mesh0.yNorm[index]
//                << " " << nse0.get_rho(i, j)
//                << " " << nse0.get_ux(i, j)
//                << " " << nse0.get_uy(i, j)
//                << " " << cde.get_C(i, j)
//                << std::endl;
//        }
//    }
//    out.close();
//}
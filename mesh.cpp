// -------------By ZhangYu HIT
// -------------1225839327@qq.com
#include <iostream>
#include "mesh.h"

using std::cout;
using std::endl;

MeshBasic::MeshBasic(int NX_, int NY_) : NX(NX_), NY(NY_)
{
    mesh = new Coordinate[NX * NY]();
    solidFlag = new LBI[NX * NY];
}

MeshBasic::~MeshBasic()
{
    delete[] mesh;
    delete[] solidFlag;
    cout << "The basic mesh is free!!!";
}

size_t MeshBasic::scalarIndex(int i, int j) const { return (NX * j + i); }

void MeshBasic::setCoordinate()
{
    for (size_t y = 0; y < NY; ++y)
    {
        for (size_t x = 0; x < NX; ++x)
        {
            size_t index = scalarIndex(x, y);
            // mesh[index].x = x * LbParameters::Cl;
            // mesh[index].y = y * LbParameters::Cl;
            mesh[index].x = x * 1;
            mesh[index].y = y * 1;
        }
    }
}

void MeshBasic::setFlag()
{
    for (size_t y = 0; y < NY; ++y)
    {
        for (size_t x = 0; x < NX; ++x)
        {
            size_t index = scalarIndex(x, y);
            solidFlag[index] = FLUID;

            if (y == 0 && x != 0 && x != NX - 1)
            {
                // --------------------------------------lower plate
                solidFlag[index] = BOTTOM;
            }
            else if (x == 0 && y == 0)
            {
                // ***********************************left bottom
                solidFlag[index] = LEFT_BOTTOM;
            }
            else if (x == NX - 1 && y == 0)
            {
                // ******************************right bottom
                solidFlag[index] = RIGHT_BOTTOM;
            }
            else if (x == 0 && y == NY - 1)
            {
                //********************************** left top
                solidFlag[index] = LEFT_TOP;
            }
            else if (x == NX - 1 && y == NY - 1)
            {
                // *********************************right top
                solidFlag[index] = RIGHT_TOP;
            }
            else if (x == 0 && y != 0 && y != NY - 1)
            {
                // -----------------------------------left plate
                solidFlag[index] = LEFT;
            }
            else if (x == NX - 1 && y != 0 && y != NY - 1)
            {
                // ---------------------------------right plate
                solidFlag[index] = RIGHT;
            }
            else if (y == NY - 1 && x != 0 && x != NX - 1)
            {
                // ---------------------------------upper plate
                solidFlag[index] = TOP;
            }
        }
    }
}

void MeshBasic::initMesh()
{
    setCoordinate();
    setFlag();
}

void MeshBasic::outputMeshTec() const
{
    std::ofstream out(filePath + "Mesh_MSG.dat");

    out << "Title= \"Lid Driven Flow\"\n";
    out << "VARIABLES = \"X\", \"Y\", \"Solid_Flag\" \n";
    out << "ZONE T= \"BOX\",I=" << NX << ",J=" << NY << ",F=	POINT" << std::endl;

    for (size_t j = 0; j < NY; ++j)
    {
        for (size_t i = 0; i < NX; ++i)
        {
            const size_t index = scalarIndex(i, j);
            out << std::fixed << std::setprecision(6)
                << " " << mesh[index].x << " " << mesh[index].y
                << std::setw(2)
                << " " << solidFlag[index]
                << std::endl;
        }
    }
}
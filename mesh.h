// -------------By ZhangYu HIT
// -------------1225839327@qq.com
#ifndef _MESH_H_
#define _MESH_H_

// #include <memory>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "lattice.h"
// #include "lattice.hh"

// using PtrFloat = std::unique_ptr<LBF[]>;
// using PtrInt = std::unique_ptr<LBI[]>;

struct Coordinate
{
    LBF x;
    LBF y;
};

class MeshBasic
{
public:
    MeshBasic(int NX_, int NY_);
    ~MeshBasic();

    void setCoordinate();
    void setFlag();
    void initMesh();
    void outputMeshTec() const;

    int getNX() const { return NX; };
    int getNY() const { return NY; };
    size_t scalarIndex(int i, int j) const;
    int getFlagOfComp(const size_t i, const size_t j) const { return solidFlag[scalarIndex(i, j)]; };
    Coordinate getCoordinate(const size_t i, const size_t j) const { return mesh[scalarIndex(i, j)]; };

private:
    enum Comp
    {
        OBSTACLE = -2,
        SILID = -1,
        FLUID,
        RIGHT,
        TOP,
        LEFT,
        BOTTOM,
        RIGHT_TOP,
        LEFT_TOP,
        LEFT_BOTTOM,
        RIGHT_BOTTOM
    };

    int NX, NY;
    Coordinate *mesh;
    LBI *solidFlag;
};

class adeD2Q5mesh : public ::MeshBasic
{
    private:
    
};

#endif
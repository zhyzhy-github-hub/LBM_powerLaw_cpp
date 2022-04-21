// -------------By ZhangYu HIT
// -------------1225839327@qq.com
#ifndef _LATTICE_H_
#define _LATTICE_H_

#include <vector>
#include <string>
#include <cmath>

typedef double scalar;
typedef int LBI;
typedef float LBF;

using std::vector;

struct Vector2D
{
    scalar x;
    scalar y;
};

const std::string filePath = "./";
// const std::string filePath = "./result/";

namespace Lattice
{
    class D2Q9Lattice
    {
    protected:
        enum D2Q9
        {
            D = 2,
            Q = 9
        };
        const double w0{4.0 / 9};
        const double wc{1.0 / 9};
        const double wd{1.0 / 36};
        const vector<int> Cx{0, 1, 0, -1, 0, 1, -1, -1, 1};
        const vector<int> Cy{0, 0, 1, 0, -1, 1, 1, -1, -1};
        const vector<double> w{w0, wc, wc, wc, wc, wd, wd, wd, wd};
        const double Cs2{1.0 / 3};
        const double Cs2inv{3.0};

        struct distributionFunction
        {
            scalar f[Q];
            scalar computeAll()
            {
                scalar all = 0;
                for (int i = 0; i < Q; ++i)
                {
                    all += f[i];
                }
                return all;
            }
        } disFun;

        struct stressTensor
        {
            scalar Dxx;
            scalar Dyy;
            scalar Dxy;

            template <typename T>
            stressTensor operator*=(const T b)
            {
                this->Dxx *= b;
                this->Dyy *= b;
                this->Dxy *= b;
                return *this;
            }
            template <typename T>
            stressTensor operator*(const T b)
            {
                stressTensor result{
                    this->Dxx * b,
                    this->Dyy * b,
                    this->Dxy * b};
                return result;
            }
            double DProductD() const { return (pow(this->Dxx, 2) + pow(this->Dyy, 2) + 2 * pow(this->Dxy, 2) + 1e-20); }
        } strTen;
    };

    class D2Q5Lattice
    {
    protected:
        enum D2Q5
        {
            D = 2,
            Q = 5
        };
        const double w0{1.0 / 3};
        const double wc{1.0 / 6};
        const vector<int> Cx{0, 1, 0, -1, 0};
        const vector<int> Cy{0, 0, 1, 0, -1};
        const vector<double> w{w0, wc, wc, wc, wc};
        const double Cs2{1.0 / 3};
        const double Cs2inv{3.0};
        // static const double w0, wc;
        // static const vector<int> Cx;
        // static const vector<int> Cy;
        // static const vector<double> w;
        // static const double cs2;
        // static const double cs2inv;

        struct distributionFunction
        {
            scalar f[Q];
            scalar computeAll()
            {
                scalar all = 0;
                for (int i = 0; i < Q; ++i)
                {
                    all += f[i];
                }
                return all;
            }
        } disFun;
    };

    template <typename T>
    inline T p2(const T x)
    {
        return (x * x);
    }

    template <typename T>
    T oneOrderExtra(const T y1, const T y2, const T q, const int dirPosition = 1)
    {
        double y0;
        y0 = y1 + q * 1;
        return y0;
    }

    // the dir is only 0 or 1
    template <typename T>
    T upwind(const T y1, const T y2, const T q, const int dirPosition = 1)
    {
        double y0;
        switch (dirPosition)
        {
        case 1 /* dirPosition = 0, it is positive direction */:
            y0 = (-2 * q + 4 * y1 - y2) / 3.0;
            break;
        case 0 /* dirPosition = 0, it is negitive direction */:
            y0 = (2 * q + 4 * y1 - y2) / 3.0;
            break;

        default:
            break;
        }
        return y0;
    }
}

#endif
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <fstream>

using namespace std;

bool FuzzyEqual(double a, double b)
{
    double maxOneAB = std::max( { 1.0, std::fabs(a) , std::fabs(b) } );
    return std::fabs(a - b) <= std::numeric_limits<double>::epsilon()*maxOneAB;
}

struct point
{
    point(double x, double y, double z)
        : x(x)
        , y(y)
        , z(z)
    {}

    const double x;
    const double y;
    const double z;
};

static bool operator==(const point & p1, const point & p2)
{
    return FuzzyEqual(p1.x, p2.x) && FuzzyEqual(p1.y, p2.y) && FuzzyEqual(p1.z, p2.z);
}

struct vec
{
    vec(double x, double y, double z)
        : x(x)
        , y(y)
        , z(z)
    {}

    vec(const point & from, const point & to)
        : x(to.x - from.x)
        , y(to.y - from.y)
        , z(to.z - from.z)
    {}

    double len() const
    {
        return std::sqrt(x*x + y*y + z*z);
    }

    const double x;
    const double y;
    const double z;
};

struct section
{
    section(const point & p1, const point &p2)
        : p1(p1)
        , p2(p2)
    {}

    const point p1;
    const point p2;
};

struct triangle
{
    triangle(const point & p1, const point & p2, const point &p3)
        : p1(p1)
        , p2(p2)
        , p3(p3)
    {}

    const point p1;
    const point p2;
    const point p3;
};


double Dot(const vec & v1, const vec & v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

vec Cross(const vec & v1, const vec & v2)
{
    return vec(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

double Triple(const vec & v1, const vec & v2, const vec & v3)
{
    return Dot(v1, Cross(v2, v3));
}

bool CheckInterval(double x, double a, double b)
{
    return (a < x && b > x) || FuzzyEqual(a, x) || FuzzyEqual(b, x);
}

bool Intersect(const point & p1, const point & p2)
{
    return p1 == p2;
}

bool Intersect(const section & s, const point & p)
{
    double dx = (s.p2.x - p.x) / (s.p2.x - s.p1.x);
    double dy = (s.p2.y - p.y) / (s.p2.y - s.p1.y);
    double dz = (s.p2.z - p.z) / (s.p2.z - s.p1.z);
    return CheckInterval(dx, 0.0, 1.0) && FuzzyEqual(dx, dy) && FuzzyEqual(dx, dz);
//    return CheckInterval(dx, 0.0, 1.0) && CheckInterval(dy, 0.0, 1.0) && CheckInterval(dz, 0.0, 1.0)&& FuzzyEqual(dx, dy) && FuzzyEqual(dx, dz) && FuzzyEqual(dy, dz);
}

bool Intersect(const point & p, const section & s)
{
    return Intersect(s, p);
}

bool Intersect(const section & s1, const section & s2)
{
    point p11 = s1.p1;
    point p12 = s1.p2;
    point p21 = s2.p1;
    point p22 = s2.p2;

    // find a and b from [0, 1]
    // a*p11 + (1-a)*p12 = b*p21 + (1-b)*p22

    double bXY =   ((p22.y - p12.y) * (p11.x - p12.x) - (p22.x - p12.x) * (p11.y - p12.y)) / ((p21.x - p22.x) * (p11.y - p12.y) - (p21.y - p22.y) * (p11.x - p12.x));
    double bXZ =   ((p22.z - p12.z) * (p11.x - p12.x) - (p22.x - p12.x) * (p11.z - p12.z)) / ((p21.x - p22.x) * (p11.z - p12.z) - (p21.z - p22.z) * (p11.x - p12.x));
    if (!FuzzyEqual(bXY, bXZ) || !CheckInterval(bXY, 0.0, 1.0))
        return false;

    double a = (bXY * (p21.x - p22.x) - (p22.x - p12.x)) / (p11.x - p12.x);
    return CheckInterval(a, 0.0, 1.0);
}

bool Intersect(const triangle & t, const point & p)
{
    vec v12(t.p2, t.p1);
    vec v13(t.p3, t.p1);
    vec v1p(p   , t.p1);
    vec v2p(p   , t.p2);
    vec v3p(p   , t.p3);
    const double totalArea = Cross(v1p, v2p).len() + Cross(v1p, v3p).len() + Cross(v2p, v3p).len(); //сумма площадей трех паралелограммов натянутых на вектора v1p, v2p и v3p
    const double trianglArea = Cross(v12, v13).len(); // удвоеная площадь треугольника
    return FuzzyEqual(totalArea, trianglArea); // totalArea >= trianglArea, равенство достигается только если точка находится внутри треугольника
    // решение может быть неустойчиво к малым смещениям точки вдоль нормали к плоскости треугольника, если проекция точки попадает в треугольник.
}

bool Intersect(const point & p, const triangle & t)
{
    return Intersect(t, p);
}

bool Intersect(const triangle & t, const section & s)
{
    vec v1(t.p1, s.p1);
    vec v2(s.p1, s.p2);
    vec v3(t.p1, t.p2);
    vec v4(t.p1, t.p3);
    const double delta = Triple(v2, v3, v4);
    if (delta == 0)
    {
        return Intersect(t, s.p1) || Intersect(section(t.p1, t.p2), s) || Intersect(section(t.p1, t.p3), s) || Intersect(section(t.p2, t.p3), s);
    }
    double tt = Triple(v1, v3, v4) / delta;
    if (!CheckInterval(tt, -1.0, 0.0))
        return false;
    double u = Triple(v2, v1, v4) / delta;
    if (!CheckInterval(u, 0.0, 1.0))
        return false;
    double v = Triple(v2, v3, v1) / delta;
    if (!CheckInterval(v, 0.0, 1.0))
        return false;
    if (!CheckInterval(v+u, 0.0, 1.0))
        return false;

    return true;
}

bool Intersect(const section & s, const triangle & t)
{
    return Intersect(t, s);
}

bool Intersect(const triangle & t1, const triangle & t2)
{
    if (Intersect(t1, section(t2.p1, t2.p2)))
        return true;
    if (Intersect(t1, section(t2.p1, t2.p3)))
        return true;
    if (Intersect(t1, section(t2.p2, t2.p3)))
        return true;

    if (Intersect(t2, section(t1.p1, t1.p2)))
        return true;
    if (Intersect(t2, section(t1.p1, t1.p3)))
        return true;
    if (Intersect(t2, section(t1.p2, t1.p3)))
        return true;

    return false;
}

section GetSection(const point & p1, const point & p2, const point &p3)
{
    assert(2 == ((p1 == p2) ? 0 : 1) + ((p1 == p3) ? 0 : 1) + ((p2 == p3) ? 0 : 1));
    return p1 == p2 ? section(p1, p3) : section(p1, p2);
}

bool Intersect(double x11, double y11, double z11, double x12, double y12, double z12, double x13, double y13, double z13,
               double x21, double y21, double z21, double x22, double y22, double z22, double x23, double y23, double z23)
{
    point p11(x11, y11, z11);
    point p12(x12, y12, z12);
    point p13(x13, y13, z13);
    point p21(x21, y21, z21);
    point p22(x22, y22, z22);
    point p23(x23, y23, z23);
    // находим топологии фигур 0 - точка, 2 - отрезок и 3 - треугольник
    int topology1 = ((p11 == p12) ? 0 : 1)
                  + ((p11 == p13) ? 0 : 1)
                  + ((p12 == p13) ? 0 : 1);

    int topology2 = ((p21 == p22) ? 0 : 1)
                  + ((p21 == p23) ? 0 : 1)
                  + ((p22 == p23) ? 0 : 1);

// разбор топологических случаев
    if (topology1 == 0 && topology2 == 0) return Intersect(p11                      , p21                      );
    if (topology1 == 0 && topology2 == 2) return Intersect(p11                      , GetSection(p21, p22, p23));
    if (topology1 == 0 && topology2 == 3) return Intersect(p11                      , triangle  (p21, p22, p23));
    if (topology1 == 2 && topology2 == 0) return Intersect(GetSection(p11, p12, p13), p21                      );
    if (topology1 == 2 && topology2 == 2) return Intersect(GetSection(p11, p12, p13), GetSection(p21, p22, p23));
    if (topology1 == 2 && topology2 == 3) return Intersect(GetSection(p11, p12, p13), triangle  (p21, p22, p23));
    if (topology1 == 3 && topology2 == 0) return Intersect(triangle  (p11, p12, p13), p21                      );
    if (topology1 == 3 && topology2 == 2) return Intersect(triangle  (p11, p12, p13), GetSection(p21, p22, p23));
    if (topology1 == 3 && topology2 == 3) return Intersect(triangle  (p11, p12, p13), triangle  (p21, p22, p23));

    assert(!"strange topologt");
    return false;
}

int main()
{
    // run test from file
    std::ifstream inp("test.txt");
    assert(inp.is_open());
    size_t n;
    inp >> n;
    std::cout << "Start test. Test in " << n << " cases." << std::endl;
    double x11, y11, z11, x12, y12, z12, x13, y13, z13, x21, y21, z21, x22, y22, z22, x23, y23, z23;
    for (size_t i = 0; i < n; ++i)
    {
         inp >> x11 >> y11 >> z11 >> x12 >> y12 >> z12 >> x13 >> y13 >> z13 >> x21 >> y21 >> z21 >> x22 >> y22 >> z22 >> x23 >> y23 >> z23;
         bool ret = Intersect(x11, y11, z11, x12, y12, z12, x13, y13, z13, x21, y21, z21, x22, y22, z22, x23, y23, z23);
         std::cout << "*   iter " << i + 1 << " result " << (ret ? "True" : "False") << std::endl;
    }

    cout << "End of test." << endl;
    return 0;
}


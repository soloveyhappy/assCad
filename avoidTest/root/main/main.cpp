#include "iostream"

#include "libavoid/geometry.h"
#include "libavoid/geomtypes.h"


int main()
{
    Avoid::Point p1(10.0, 10.0);
    Avoid::Point p2(40.0, 30.0);
    Avoid::Point p3(30.0, 50.0);
    Avoid::Point proj = Avoid::projection(p1, p2, p3);
    Avoid::Point p = p1 + p2;

    std::cout << p.x << " ; " << p.y;

    return 0;
}


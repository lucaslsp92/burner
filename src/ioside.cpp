#ifndef H_IOSIDE
#define H_IOSIDE

#include "Point.h"
#include "Particle.h"

int Particle::ioside (double DP, double geometry[])
{
    Point P (r, DP);

    double x = P.x;
    double y = P.y;
    double z = P.z;

    Point O;

    Point X (1,0,0);
    Point Y (0,1,0);
    Point Z (0,0,1);

    double H = 0.5;
    double L = 0.2; 
    double D = 10*DP;

    Point A (0, 0, -5*DP); 

    Region fluid = P.rectangleXY(L, H-DP, DP);

    if (fluid) 
        return 0;

    return 2;
}
#endif
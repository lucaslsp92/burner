#ifndef H_IOFLOW
#define H_IOFLOW

#include "Point.h"
#include "Particle.h"

int Particle::ioflow (double DP, double geometry[])
{
    Point P (r, DP);

    double x = P.x;
    double y = P.y;
    double z = P.z;

    Point O;
    Point X (1,0,0);
    Point Y (0,1,0);
    Point Z (0,0,1);

    //////////////////////////////////////////
    //// Lid-driven cavity flow 2D
    //// Variaveis (double, Point, etc..)
    double D = 0.15;        // cavity length
    double W = DP;         // cavity width

    if (z < 0.0 || z >= W)
    {
        return -1;
    }

    if (y > D)
    {
        return id+2;
    }

    /*if (x < 0.0 && y >= 0.0 && y <= D)
    {
        return id+2;
    }*/

    ///Return padrão (caso nenhuma operação seja feita)
    return id;

}

#endif
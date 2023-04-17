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
    //// Hydrostatic tank 3D
    ///Variaveis (double, Point, etc..)
    double H = 0.15;         // water column
    double L = 0.50;         // tank length
    double W = 1.0;          // tank width

    if (x > L*0.5-1*DP && x < L*0.5+1*DP && y > DP && y < 1.3*H-DP && z > 0.0 && z <= DP)
    {
        return id+1;
    }

    if (y >= 1.3*H)
    {
        return -1;
    }

    if (z < 0.0 || z >= DP)
    {
        return -1;
    }

    ///Return padrão (caso nenhuma operação seja feita)
    return id;

}

#endif
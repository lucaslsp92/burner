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
    double H = 11.5;         // water column
    double L = 10.0;         // tank length
    double W = 2.0;          // tank width

    if (y >= 1.1*H)
    {
        return -1;
    }

    /*if (z < 0.0 || z >= DP)
    {
        return -1;
    }*/

    ///Return padrão (caso nenhuma operação seja feita)
    return id;

}

#endif
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
    //// Dam break tank 2D
    double H = 0.07;         // water column
    double W = 0.32;         // water width
    double L = 0.7;          // tank length

    if (y >= H+3.0*DP)
    {
        return -1;
    }

    if(z <= 0.0 || z > DP)
    {
        return -1;
    }

    ///Return padrão (caso nenhuma operação seja feita)
    return id;

}

#endif
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
    //// Free-surface poiseuille flow
    //// Variaveis (double, Point, etc..)
    double H = 1.0;        // channel height
    double L = 2*H;        // channel lenght

    if(z < 0.0 || z >= DP)
    {
        return -1;
    }

    if(y >= 1.1*H)
    {
        return -1;
    }

    if(x < 0 /*&& y > 0 && y < H*/)
    {
        //return id+2;
        return -1;
    }

    if(x >= L /*&& y > 0*/)
    {
        //return 6;
        //return id+4;
        return -1;
    }

    ///Return padrão (caso nenhuma operação seja feita)
    return id;

}

#endif
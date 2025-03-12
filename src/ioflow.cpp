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
    double L = 0.05; 
    double D = 0.005;
    double hin = 0.005;
    double H = 0.06;

    if(x > 0.0 && y > 0.0 && x < L && y < L && z < 0.0)
        return id+2;

    /*if((x-L/2.0)*(x-L/2.0) + (y-L/2.0)*(y-L/2.0) > D/2.0*D/2.0 && z > H)
        return id+4;*/

    if((x-L/2.0)*(x-L/2.0) + (y-L/2.0)*(y-L/2.0) < D/2.0*D/2.0 && z > H)
        return id+6;

    ///Return padrão (caso nenhuma operação seja feita)
    return id;

}

#endif
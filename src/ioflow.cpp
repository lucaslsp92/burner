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
    //// Moonpool sloshing
    //// Variaveis (double, Point, etc..)
    double cresthillHeight = 0.028;
    double crestHillDistance = 9.0*cresthillHeight;
    double channelHeight = 3.035*cresthillHeight;

    /// Points

    ///Regiões (Region)

    ///Operações
    if(x<0.0 || x>crestHillDistance)
    {
        return -1;
    }

    if(z<0.0 || z>=DP)
    {
        return -1;
    }

    if(y>channelHeight)
    {
        return id+2;
    }

    ///Return padrão (caso nenhuma operação seja feita)
    return id;

}

#endif
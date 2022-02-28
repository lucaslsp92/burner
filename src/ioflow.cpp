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
    double waterDepth = 0.50;
    double tankHeight = waterDepth+0.3;

    /// Points

    ///Regi�es (Region)

    ///Opera��es
    if(y>tankHeight-waterDepth && (id==2 || id==3))
    {
        return -1;
    }

    if(z<0.0 || z>=DP)
    {
        return -1;
    }

    ///Return padr�o (caso nenhuma opera��o seja feita)
    return id;

}

#endif
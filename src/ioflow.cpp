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
    double cylinderLength = 0.1;
    double waterDepth = 9.0*cylinderLength;
    double tankLength = 30.0*cylinderLength;
    double tankHeight = waterDepth + 2.0*cylinderLength;

    /// Points

    ///Regiões (Region)

    ///Operações
    if(x<0.0 && y>0.0 && y<waterDepth && z>=0.0 && z<DP && (id==2 || id==3))
    {
        return id+4;
    }

    if(x>tankLength && z>=0.0 && z<DP && (id==2 || id==3))
    {
        return id+6;
    }

    /*if(x>tankLength)
    {
        return -1;
    }*/

    if(y>tankHeight && (id==2 || id==3))
    {
        return -1;
    }

    if(z<0.0 || z>=DP)
    {
        return -1;
    }

    ///Return padrão (caso nenhuma operação seja feita)
    return id;

}

#endif
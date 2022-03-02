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
    double cylinderDiameter = 1.0;
    double cylinderPosition = 9.5;
    double waterLength = 29.5;
    double waterHeight = 10.0;

    /// Points

    ///Regiões (Region)

    ///Operações
    if(x<0.0 || x>waterLength || y<0.0 || y>waterHeight)
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
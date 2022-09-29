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
    double centerX = 46.78;
    double centerY = 25.455;        /// DP = 0.005 m
    //double centerY = 25.455;         /// DP = 0.010 m
    double centerZ = 28.5388;
    double dumpValveRadius = 0.100;
    double dumpValveHeight = 0.120;
    double deckWidth = 0.800;
    double dumpValveOpening = 0.060;
    double inletRadius = 0.095;
    double inletHeight = 0.180;


    /// Points

    ///Regiões (Region)

    ///Operações

    ///Return padrão (caso nenhuma operação seja feita)
    if ((z-centerX)*(z-centerX)+(y-centerY)*(y-centerY) < (inletRadius)*(inletRadius) && (x-centerZ+inletHeight) < 0)
    {
        return id+2;
    }

    if (((y-centerY) < -deckWidth/2 || (y-centerY) > deckWidth/2 || (z-centerX) < -deckWidth/2 || (z-centerX) > deckWidth/2) || (x-centerZ) > dumpValveHeight+deckWidth)
    {
        return -1;
    }

    if((x-centerZ) > dumpValveHeight)
    {
        return -1;
    }

    return id;
}
#endif
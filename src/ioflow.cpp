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

    ///Variaveis (double, Point, etc..)
    double L = 10.0;        // pipe length
    double D = 0.078;       // pipe diameter
    double R = 1.5*D/2;     // elbow revolution radius
    double H = 0.01;        // inflow pipe length

    Point inflowBeginn (0,R+D/2+H,0);

    ///Regiões (Region)
    Region inflow = P.transformation(inflowBeginn,Y).cylinder(D/2, 3*DP);

    ///Operações
    if (inflow)
    {
        return id+2;
    }

    if (x >= L)
    {
        return -1;
    }

    ///Return padrão (caso nenhuma operação seja feita)

    return id;

}

#endif
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
    double tankLength = 1500.0;
    double tankHeight = 292.0;
    double waterHeight = 250.0;
    double beachLength = 200.0;

    ///Operações
    if(y>=tankHeight)
        return -1;

    if(z<0.0 || z>=DP)
        return -1;

    if(x<0.0 && y>129.0 && (id==2 || id==3))
        return id+2;

    if(x>850 && x<950 && y>0 && (id==2 || id==3))
        return id+4;

    ///Return padrão (caso nenhuma operação seja feita)
    return id;

}

#endif
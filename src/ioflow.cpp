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
/*
    //// Dam-break 2D ////
    double Wy = 1.0;
    //// Dam-break 2D ////
*/
/*
    //// Dam-break 3D with internal solid ////
    double Wz = 1.0;
    //// Dam-break 3D with internal solid ////
*/
    ///Regiões (Region)
    if (id==2 || id==3)
    {
/*
        //// Dam-break 2D ////
        if (y > Wy) return -1;    // Remove wall top
        //// Dam-break 2D ////
*/
/*
        //// Dam-break 3D with internal solid ////
        if (z > Wz) return -1;      // Remove wall top
        //// Dam-break 3D with internal solid ////
*/
    }
    ///Operações

    ///Return padrão (caso nenhuma operação seja feita)

    return id;

}

#endif

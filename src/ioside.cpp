#ifndef H_IOSIDE
#define H_IOSIDE

#include "Point.h"
#include "Particle.h"

int Particle::ioside (double DP, double geometry[])
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
    //// Lid-driven cavity flow 2D
    //// Variaveis (double, Point, etc..)
    double D = 0.15;        // cavity length
    double W = DP;         // cavity width

    ///Regiões (Region)
    Region fluid = P.rectangleXY(D, D, D);

    ///Operações
    if (fluid) 
        return 0;

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif
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
    //// Poiseuille flow 2D
    //// Variaveis (double, Point, etc..)
    double D = 0.1;         // pipe diameter
    double L = 0.5;         // pipe length

    ///Regiões (Region)
    Region fluid = P.rectangleXY(L, D, DP);

    ///Operações
    if (fluid) 
        return 0;

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif
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
    //// Dam break tank 2D
    //// Variaveis (double, Point, etc..)
    double L = 0.05; 
    double D = 0.005;
    double hin = 0.005;
    double H = 0.06;

    Point A(L/2.0, L/2.0, H);

    ///Regiões (Region)
    Region table = P.rectangleXY(L, L, L);
    Region inflow = P.transformation(A).cylinder(D/2.0, hin);

    ///Operações
    if(inflow)
      return -1;
    if (table)
       return -1;

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif
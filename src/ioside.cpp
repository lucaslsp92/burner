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

    //// Variaveis (double, Point, etc..)
    double H = 0.2;         // pipe diameter
    double L = 0.5;         // pipe length

    Point A(0.2075, 0.1, -5.0*DP);
    Point B(0.2250, 0.1, -5.0*DP);

    ///Regiões (Region)
    Region fluid = P.rectangleXY(L, H, DP);
    Region tank = P.rectangleXY(L, 1.1*H, 10.0*DP);
    Region cup = P.transformation(A).rectangleXY(0.085, 0.1175, 10.0*DP); 
    Region vazio = P.transformation(B).rectangleXY(0.05, 0.1, 10.0*DP);

    ///Operações
    if(cup && !vazio)
        return 4;
    if (fluid) 
        return 0;    
    if(tank)
        return -1;

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif
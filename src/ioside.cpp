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
    double H = 0.4;         // water column
    double W = 0.2;         // water length
    double L = 0.8;         // tank length
    double th = 0.5;        // tank height

    ///Regiões (Region)
    Region fluid = P.rectangleXY(W, H, DP);
    Region tank  = P.rectangleXY(L, th, 10*DP);

    ///Operações
    if (tank)
    {
        if (fluid) 
            return 0;

    	return -1;
    }

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif
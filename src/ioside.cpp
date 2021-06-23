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
    double H = 0.45;         // water column
    double W = 0.17;         // water length
    double L = 0.7;          // tank length
    double th = 0.60;        // tank height
    double Hb = 0.00;        // water height tank bottom
    double Wb = 0.0;         // water length tank bottom
    double tw = 0.1;         // tank width

    ///Regiões (Region)
    Region fluid = P.rectangleXY(W, H, tw);
    Region bottomFluid = P.rectangleXY(Wb, Hb, tw);
    Region tank  = P.rectangleXY(L, th, tw);

    ///Operações
    if (tank)
    {
        if (fluid || bottomFluid) 
            return 0;

    	return -1;
    }

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif
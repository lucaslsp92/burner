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
    //// Hydrostatic tank 3D
    ///Variaveis (double, Point, etc..)
    double H = 3.0;         // water column
    double L = 6.0;         // tank length
    double W = DP;          // tank width
    //double W = 2.5+DP;         // tank width  

    Point A (0, 0, 0); 

    ///Regi�es (Region)
    Region fluid = P.transformation(A).rectangleXY(L/2.0, H, W);
    Region tank  = P.rectangleXY(L, H*1.1, W);

    ///Opera��es
    if (tank)
    {
        if (fluid) 
            return 0;

    	return -1;
    }

    ///Return padr�o (constr�i a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif
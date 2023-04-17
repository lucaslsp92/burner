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
    double H = 0.15;         // water column
    double L = 0.50;         // tank length
    double W = 1.0;          // tank width

    Point R (L*0.5-2*DP, 0.0, 0.0);

    ///Regiões (Region)
    Region fluid = P.rectangleXY(L, H, DP);
    Region barrier = P.transformation(R).rectangleXY(4*DP, H*1.3, W);
    Region tankLeft  = P.rectangleXY(L*0.5-2*DP, H*1.3, W);
    Region tankRight = P.transformation(R).rectangleXY(L*0.5+2*DP, H*1.3, W);

    ///Operações
    if (tankLeft || tankRight)
    {
        if (barrier)
            return 4;
        else if (fluid && x < L*0.5-2*DP) 
            return 0;
        else if (fluid && x > L*0.5+2*DP)
            return 1; 

    	return -1;
    }

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif
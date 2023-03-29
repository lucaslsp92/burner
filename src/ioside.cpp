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
    //// Moonpool sloshing
    //// Variaveis (double, Point, etc..)
    double cylinderRadius = 0.08;
    double cylinderDepth = 3.0*cylinderRadius;
    double cylinderPosition = 114.0*cylinderRadius;
    double waterDepth = 8.0*cylinderRadius;
    double waterLength = 125.0*cylinderRadius;
    double tankLength = 125.0*cylinderRadius;
    double tankHeight = 10.0*cylinderRadius;

    /// Points
    Point tankCorner(0.0,0.0,-4*DP);
    Point cylinderCenter(cylinderPosition,waterDepth-cylinderDepth,-4*DP);

    ///Regiões (Region)
    Region water = P.rectangleXY(waterLength,waterDepth,DP); 
    Region tank = P.transformation(tankCorner).rectangleXY(tankLength,tankHeight,8*DP);
    Region cylinder = P.transformation(cylinderCenter).cylinder(cylinderRadius,8*DP);

    ///Operações
    if(tank)
    {
        if(cylinder)
            return 4;
        if(water)
            return 0;
        else
            return -1;
    }  

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif
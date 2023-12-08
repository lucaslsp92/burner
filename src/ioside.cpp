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
    double cylinderPosition = 20.0*cylinderRadius;
    double waterDepth = 16.0*cylinderRadius;
    double waterLength = 60.0*cylinderRadius;
    double tankLength = 60.0*cylinderRadius;
    double tankHeight = waterDepth + 2.0*cylinderRadius;

    /// Points
    Point tankCorner(0.0,0.0,-4*DP);
    Point cylinderCenter(cylinderPosition,waterDepth-cylinderDepth,-4*DP);

    ///Regi�es (Region)
    Region water = P.rectangleXY(waterLength,waterDepth,DP); 
    Region tank = P.transformation(tankCorner).rectangleXY(tankLength,tankHeight,8*DP);
    Region cylinder = P.transformation(cylinderCenter).cylinder(cylinderRadius,8*DP);

    ///Opera��es
    if(tank)
    {
        //if(cylinder)
            //return 4;
        if((x-cylinderPosition)*(x-cylinderPosition)+(y-waterDepth+cylinderDepth)*(y-waterDepth+cylinderDepth) < (cylinderRadius+DP)*(cylinderRadius+DP))
            return -1;
        if(water)
            return 0;
        else
            return -1;
    }  

    ///Return padr�o (constr�i a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif
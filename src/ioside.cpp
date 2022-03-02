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
    double cylinderDiameter = 1.0;
    double cylinderPosition = 9.5;
    double waterLength = 29.5;
    double waterHeight = 10.0;

    /// Points
    Point cylinderCenter(cylinderPosition,waterHeight/2,-4*DP);

    ///Regiões (Region)
    Region water = P.rectangleXY(waterLength,waterHeight,DP);
    Region cylinder = P.transformation(cylinderCenter).cylinder(cylinderDiameter/2,8*DP);

    ///Operações
    if(cylinder)
    {
        return 2;
    }
    if(water)
    {
        return 0;
    }

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif
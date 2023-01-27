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
    double waterDepth = 1.00;
    double waterLength = 2.2;
    double tankHeight = 1.3;
    double tankLength = 2.2;
    double draft = 0.11;
    double vesselLength = 0.40;
    double vesselHeight = 0.40;
    double moonpoolLength = 0.20;
    double recessLength = 0.10;
    double recessHeight = 0.05;
    double beachRatio = 2.5;
    double beachAngle = atan(1.0/beachRatio)*180/M_PI;
    double beachLength = 2.0*sqrt((waterDepth*waterDepth)+(beachRatio*waterDepth*beachRatio*waterDepth));

    /// Points
    Point tankVertex(-tankLength/2,-waterDepth,-4*DP);
    Point waterVertex(-waterLength/2,-waterDepth,0);
    Point leftVesselVertex(-vesselLength-moonpoolLength/2,-draft,-4*DP);
    Point rigthVesselVertex(moonpoolLength/2,-draft,-4*DP);
    Point recessVertex(0,-draft,-4*DP);
    Point leftBeachVertex(-waterLength/2,-waterDepth,-4*DP);
    Point rightBeachVertex(waterLength/2,-waterDepth,-4*DP);
    Point beachWaterVertex(-waterLength/2-beachRatio*waterDepth,-waterDepth,0);

    ///Regiões (Region)
    Region tank = P.transformation(tankVertex).rectangleXY(tankLength,tankHeight,8*DP);
    Region water = P.transformation(waterVertex).rectangleXY(waterLength,waterDepth,DP);
    Region leftVessel = P.transformation(leftVesselVertex).rectangleXY(vesselLength,vesselHeight,8*DP);
    Region rigthVessel = P.transformation(rigthVesselVertex).rectangleXY(vesselLength,vesselHeight,8*DP);
    Region recess = P.transformation(recessVertex).rectangleXY(recessLength,recessHeight,8*DP);
    Region leftBeach = P.transformation(leftBeachVertex,Z,-(90-beachAngle)).rectangleXY(beachLength,beachLength,8*DP);
    Region rightBeach = P.transformation(rightBeachVertex,Z,-beachAngle).rectangleXY(beachLength,beachLength,8*DP);
    Region beachWater = P.transformation(beachWaterVertex).rectangleXY(waterLength+2*beachRatio*waterDepth,waterDepth,DP);

    ///Operações
    if(leftVessel || rigthVessel || recess)
    {
        return 4;
    }
    if(water /*|| (beachWater && leftBeach) || (beachWater && rightBeach)*/)
    {
        return 0;
    }
    if(tank /*|| leftBeach || rightBeach*/)
    {
        return -1;
    }    

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif
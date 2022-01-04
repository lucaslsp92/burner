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
    double waterDepth = 1.0;
    double waterLength = 5.0;
    double tankHeight = 1.3;
    double tankLength = 8.0;
    double draft = 0.11;
    double vesselLength = 0.40;
    double vesselHeight = 0.30;
    double moonpoolLength = 0.20;
    double recessLength = 0.10;
    double recessHeight = 0.05;
    double beachLength = 1.5;
    double beachHypotenuse = 1.5*sqrt(beachLength*beachLength+tankHeight*tankHeight);
    double leftBeachAngle = atan(waterDepth/beachLength)*180/M_PI-180;
    double rightBeachAngle = 90-atan(waterDepth/beachLength)*180/M_PI;

    /// Points
    Point tankVertex(-tankLength/2,-waterDepth,-4*DP);
    Point waterVertex(-waterLength/2,-waterDepth,0);
    Point leftVesselVertex(-vesselLength-moonpoolLength/2,-draft,-4*DP);
    Point rigthVesselVertex(moonpoolLength/2,-draft,-4*DP);
    Point recessVertex(0,-draft,-4*DP);
    Point leftBeachVertex(-waterLength/2+beachLength,-waterDepth,-4*DP);
    Point rightBeachVertex(waterLength/2-beachLength,-waterDepth,-4*DP);

    ///Regiões (Region)
    Region tank = P.transformation(tankVertex).rectangleXY(tankLength,tankHeight,8*DP);
    Region water = P.transformation(waterVertex).rectangleXY(waterLength,waterDepth,DP);
    Region leftVessel = P.transformation(leftVesselVertex).rectangleXY(vesselLength,vesselHeight,8*DP);
    Region rigthVessel = P.transformation(rigthVesselVertex).rectangleXY(vesselLength,vesselHeight,8*DP);
    Region recess = P.transformation(recessVertex).rectangleXY(recessLength,recessHeight,8*DP);
    Region leftBeach = P.transformation(leftBeachVertex,Z,leftBeachAngle).rectangleXY(beachHypotenuse,beachHypotenuse,8*DP);
    Region rightBeach = P.transformation(rightBeachVertex,Z,rightBeachAngle).rectangleXY(beachHypotenuse,beachHypotenuse,8*DP);

    ///Operações
    if(leftVessel || rigthVessel || recess)
    {
        return 4;
    }
    if(water && !leftBeach && !rightBeach)
    {
        return 0;
    }
    if(tank && !leftBeach && !rightBeach)
    {
        return -1;
    }    

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif
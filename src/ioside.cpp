#ifndef H_IOSIDE
#define H_IOSIDE

#include "Point.h"
#include "Particle.h"
#include <algorithm>

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
    double centerX = 46.78;
    double centerY = 25.455;
    double centerZ = 28.5388;
    double dumpValveRadius = 0.100;
    double dumpValveHeight = 0.120;
    double deckWidth = 0.800;
    double dumpValveOpening = 0.060;
    double supportWidth = 0.05;
    double inletRadius = 0.095;
    double inletHeight = 0.180;

    /// Points
    Point DVC(centerZ,centerY,centerX);
    Point DCL(dumpValveHeight+centerZ,-deckWidth/2+centerY,deckWidth/2+centerX);
    Point DSK(dumpValveHeight+dumpValveOpening+centerZ,centerY,centerX);
    Point SPL(dumpValveHeight+centerZ,-supportWidth/2+centerY,-dumpValveRadius+centerX);
    Point SPR(dumpValveHeight+centerZ,-supportWidth/2+centerY,dumpValveRadius+3*DP+centerX);
    Point INC(centerZ-inletHeight,centerY,centerX);
    
    ///Regiões (Region)
    Region dumpValve = P.transformation(DVC,X).cylinder(dumpValveRadius,dumpValveHeight);
    Region deck = P.transformation(DCL,X).rectangleXY(deckWidth,deckWidth,deckWidth);
    Region disk = P.transformation(DSK,X).cylinder(dumpValveRadius,3*DP);
    Region supportLeft = P.transformation(SPL,X).rectangleXY(3*DP,supportWidth,dumpValveOpening+3*DP);
    Region supportRight = P.transformation(SPR,X).rectangleXY(3*DP,supportWidth,dumpValveOpening+3*DP);
    Region inlet = P.transformation(INC,X).cylinder(inletRadius,inletHeight);

    ///Operações
    if(/*dumpValve || inlet ||*/ (deck /*&& !(disk) && !(supportLeft) && !(supportRight)*/))
    {
        return -1;
    }

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif
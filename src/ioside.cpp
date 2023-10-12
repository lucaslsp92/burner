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
    double cylinderLength = 0.1;
    double cylinderWidth = 1.0*cylinderLength;
    double cylinderDepth = 3.0*cylinderLength + 0.5*cylinderLength;
    double cylinderPosition = 10.0*cylinderLength;
    double waterDepth = 9.0*cylinderLength;
    double waterLength = 30.0*cylinderLength;
    double tankLength = 30.0*cylinderLength;
    double tankHeight = waterDepth + 2.0*cylinderLength;

    int nFluid = 149;

    /// Points
    Point tankCorner(0.0,0.0,-4*DP);
    Point cylinderCenter(cylinderPosition,waterDepth-cylinderDepth,0.0);
    Point correction1Corner(0.0,waterDepth+0.0*DP,0.0);
    Point correction2Corner(0.0,waterDepth+1.0*DP,0.0);
    Point correction3Corner(0.0,waterDepth+2.0*DP,0.0);
    Point correction4Corner(0.0,waterDepth+3.0*DP,0.0);
    Point correction5Corner(0.0,waterDepth+4.0*DP,0.0);

    ///Regiões (Region)
    Region water = P.rectangleXY(waterLength,waterDepth,DP); 
    Region tank = P.transformation(tankCorner).rectangleXY(tankLength,tankHeight,8*DP);
    Region cylinder = P.transformation(cylinderCenter).rectangleXYCenter(cylinderWidth,cylinderLength,8*DP);
    Region correction1, correction2, correction3, correction4, correction5;
    if(tankLength/DP > nFluid)
    {
        correction1 = P.transformation(correction1Corner).rectangleXY(nFluid*DP,DP,DP);
        correction2 = P.transformation(correction2Corner).rectangleXY(0*DP,DP,DP);
        correction3 = P.transformation(correction3Corner).rectangleXY(0*DP,DP,DP);
        correction4 = P.transformation(correction4Corner).rectangleXY(0*DP,DP,DP);
        correction5 = P.transformation(correction5Corner).rectangleXY(0*DP,DP,DP);
    }
    else if(2.0*tankLength/DP > nFluid)
    {
       correction1 = P.transformation(correction1Corner).rectangleXY(int(tankLength/DP)*DP,DP,DP); 
       correction2 = P.transformation(correction2Corner).rectangleXY((nFluid-int(tankLength/DP))*DP,DP,DP); 
       correction3 = P.transformation(correction3Corner).rectangleXY(0*DP,DP,DP);
       correction4 = P.transformation(correction4Corner).rectangleXY(0*DP,DP,DP);
       correction5 = P.transformation(correction5Corner).rectangleXY(0*DP,DP,DP);
    }
    else if(3.0*tankLength/DP > nFluid)
    {
       correction1 = P.transformation(correction1Corner).rectangleXY(int(tankLength/DP)*DP,DP,DP); 
       correction2 = P.transformation(correction2Corner).rectangleXY(int(tankLength/DP)*DP,DP,DP); 
       correction3 = P.transformation(correction3Corner).rectangleXY((nFluid-int(2.0*tankLength/DP))*DP,DP,DP);
       correction4 = P.transformation(correction4Corner).rectangleXY(0*DP,DP,DP);
       correction5 = P.transformation(correction5Corner).rectangleXY(0*DP,DP,DP);
    }
    else if(4.0*tankLength/DP > nFluid)
    {
       correction1 = P.transformation(correction1Corner).rectangleXY(int(tankLength/DP)*DP,DP,DP); 
       correction2 = P.transformation(correction2Corner).rectangleXY(int(tankLength/DP)*DP,DP,DP); 
       correction3 = P.transformation(correction3Corner).rectangleXY(int(tankLength/DP)*DP,DP,DP);
       correction4 = P.transformation(correction4Corner).rectangleXY((nFluid-int(3.0*tankLength/DP))*DP,DP,DP);
       correction5 = P.transformation(correction5Corner).rectangleXY(0*DP,DP,DP);
    }
    else
    {
       correction1 = P.transformation(correction1Corner).rectangleXY(int(tankLength/DP)*DP,DP,DP); 
       correction2 = P.transformation(correction2Corner).rectangleXY(int(tankLength/DP)*DP,DP,DP); 
       correction3 = P.transformation(correction3Corner).rectangleXY(int(tankLength/DP)*DP,DP,DP);
       correction4 = P.transformation(correction4Corner).rectangleXY(int(tankLength/DP)*DP,DP,DP);
       correction5 = P.transformation(correction5Corner).rectangleXY((nFluid-int(4.0*tankLength/DP))*DP,DP,DP);
    }

    ///Operações
    if(tank)
    {
        if(cylinder)
            return 4;
        //if((x-cylinderPosition)*(x-cylinderPosition)+(y-waterDepth+cylinderDepth)*(y-waterDepth+cylinderDepth) < (cylinderRadius+DP)*(cylinderRadius+DP))
            //return -1;
        if(water)
            return 0;
        //if(correction1 || correction2 || correction3 || correction4 || correction5)
            //return 0;
        else
            return -1;
    }  

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif
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
    //// Pool weir fishway
    //// Variaveis (double, Point, etc..)
    double waterDepth = 0.23;
    double tankLength = 3.8;
    double tankHeight = 0.3;
    double weirHeight = 0.14;
    double weirWidth;
    double weirSpacing;
    double firstWeirPosition;
    double beachStart;
    double beachSlope = 1000.0;
    double beachLength = 1.0;
    double beachEnd;

    if(DP > 0.005)
        weirSpacing = 0.33;
    else
        weirSpacing = 0.335;

    if(6*DP > 0.01)
        weirWidth = 6*DP;
    else
        weirWidth = 0.01;    
    
    firstWeirPosition = 1.0+weirSpacing-weirWidth/2.0;
    beachStart = firstWeirPosition+6.0*weirSpacing;
    beachEnd = beachStart+beachLength;

    /// Points
    Point Ref(0, 0, 1);
    Point WP0(firstWeirPosition+0.0*weirSpacing, 0.0, -4*DP);
    Point WP1(firstWeirPosition+1.0*weirSpacing, 0.0, -4*DP);
    Point WP2(firstWeirPosition+2.0*weirSpacing, 0.0, -4*DP);
    Point WP3(firstWeirPosition+3.0*weirSpacing, 0.0, -4*DP);
    Point WP4(firstWeirPosition+4.0*weirSpacing, 0.0, -4*DP);
    Point WP5(firstWeirPosition+5.0*weirSpacing, 0.0, -4*DP);

    ///Regioes (Region)
    Region water = P.rectangleXY(tankLength,waterDepth,DP); 
    Region tank = P.rectangleXY(tankLength,tankHeight,8*DP);
    Region weir0 = P.transformation(WP0).rectangleXY(weirWidth,weirHeight,8*DP);
    Region weir1 = P.transformation(WP1).rectangleXY(weirWidth,weirHeight,8*DP);
    Region weir2 = P.transformation(WP2).rectangleXY(weirWidth,weirHeight,8*DP);
    Region weir3 = P.transformation(WP3).rectangleXY(weirWidth,weirHeight,8*DP);
    Region weir4 = P.transformation(WP4).rectangleXY(weirWidth,weirHeight,8*DP);
    Region weir5 = P.transformation(WP5).rectangleXY(weirWidth,weirHeight,8*DP);

    ///Operacoes
    if(tank)
    {
        if(water && !(weir0) && !(weir1) && !(weir2) && !(weir3) && !(weir4) && !(weir5) && !(-beachSlope*x + y + beachSlope*beachStart < 0.0 && y < weirHeight))
            return 0;
        else if(!(weir0) && !(weir1) && !(weir2) && !(weir3) && !(weir4) && !(weir5) && !(-beachSlope*x + y + beachSlope*beachStart < 0.0 && y < weirHeight))
            return -1;
    }  

    ///Return padrao (constroi a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif
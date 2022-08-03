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

    Point BRC(10.983,0,0);

    //////////////////////////////////////////
    //// Moonpool sloshing
    //// Variaveis (double, Point, etc..)
    double tankHeight = 27.856;
    double tankBreadth = 15.302;
    double tankFilling = 27.856 - 1.355;

    /// Points
    
    ///Regiões (Region)
    Region slopTankExternal = P.rectangleXY(tankBreadth,tankHeight,DP);
    Region belowLeftCurve = P.rectangleXY(0.765,20.229,DP);
    Region belowRightCurve = P.transformation(BRC).rectangleXY(4.319,3.858,DP);

    ///Operações
    if(slopTankExternal && y >= 25.0)
    {        
        if (((-0.780*x + 15.302*y -414.316952 > 0 && x > 0.000 && x < 15.302) ||                     /// TOP SLOPE
            ( 2.884*x +  1.811*y - 34.025530 < 0 && x > 0.765 && x <  2.576) ||                      /// LINEAR EQ. 1
            ( 2.454*x +  1.741*y - 31.889830 < 0 && x > 2.576 && x <  4.317) ||                      /// LINEAR EQ. 2
            ( 2.333*x +  1.752*y - 31.502025 < 0 && x > 4.317 && x <  6.069) ||                      /// LINEAR EQ. 3
            ( 2.122*x +  1.688*y - 29.587930 < 0 && x > 6.069 && x <  7.757) ||                      /// LINEAR EQ. 4
            ( 3.919*x +  3.226*y - 55.488285 < 0 && x > 7.757 && x < 10.983) ||                      /// LINEAR EQ. 5
            ( belowLeftCurve  && !((x-4.998)*(x-4.998) + (y-20.229)*(y-20.229) < 5.0*5.0)) ||        /// CIRCLE EQ. 1
            ( belowRightCurve && !((x-27.287)*(x-27.287) + (y-17.758)*(y-17.758) < 21.4*21.4))))     /// CIRCLE EQ. 2
        {
        }
        else if (y > tankFilling)
        {
            return -1;
        }
        else
        {
            return 0;
        }
    }

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif